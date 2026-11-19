"""Offline-first cache of public reference data (UniProt, PDBe/SIFTS, AlphaFold DB, ConSurf-DB).

Every getter answers from data/reference/ when it can. A miss fetches the file once,
stores the raw bytes and records url, date, sha256, size and HTTP status in
data/reference/manifest.json. BINDERSCOPE_OFFLINE=1 turns a miss into a clear error.

Only accessions, PDB ids, gene symbols, taxon ids and chain ids are ever put in a request.
Request builders validate each of them, so a binder sequence or a design id cannot end up
in a URL by accident.
"""

from __future__ import annotations

import contextlib
import csv
import datetime as _dt
import fcntl
import hashlib
import io
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CACHE_DIR = PACKAGE_ROOT / "data" / "reference"
USER_AGENT = "binderscope-refdata/0.1 (academic protein-design triage tool; public reference data only)"
OFFLINE_ENV = "BINDERSCOPE_OFFLINE"
CACHE_ENV = "BINDERSCOPE_REFDATA_DIR"

TIMEOUT_S = 60
RETRIES = 3
BACKOFF_S = 2.0
PACING_S = 0.4
PAE_MAX_BYTES = 4_000_000

ALLOWED_HOSTS = {
    "rest.uniprot.org",
    "www.ebi.ac.uk",
    "alphafold.ebi.ac.uk",
    "consurfdb.tau.ac.il",
    "files.rcsb.org",
    "rest.ensembl.org",
    "cdnjs.cloudflare.com",
}

_ACC_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
_PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
_CHAIN_RE = re.compile(r"^[A-Za-z0-9]{1,4}$")
_GENE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,19}$")

SIFTS_KINDS = ("uniprot", "uniprot_segments", "pfam", "interpro", "cath")

UNIPROT_SEARCH_FIELDS = (
    "accession,reviewed,id,gene_names,organism_name,organism_id,length,protein_name,fragment,sequence"
)


class RefDataError(RuntimeError):
    """Reference data could not be obtained."""


class RefDataOffline(RefDataError):
    """A cache miss happened while BINDERSCOPE_OFFLINE is set."""


class RefDataNotFound(RefDataError):
    """The remote server answered 404 for this record (recorded in the manifest)."""


def is_offline() -> bool:
    return os.environ.get(OFFLINE_ENV, "").strip().lower() not in ("", "0", "false", "no")


def cache_dir() -> Path:
    override = os.environ.get(CACHE_ENV)
    return Path(override) if override else DEFAULT_CACHE_DIR


def manifest_path() -> Path:
    return cache_dir() / "manifest.json"


# ---------------------------------------------------------------- validation

def check_accession(acc: str) -> str:
    acc = str(acc).strip().upper()
    if not _ACC_RE.match(acc):
        raise ValueError(f"not a UniProt accession: {acc!r}")
    return acc


def is_pdb_id(pdb_id: object) -> bool:
    """Whether this names a PDB entry, so a caller can skip a lookup instead of failing one.

    A predicted model (`af-q9nzq7-f1-model_v6`) has no PDB entry and therefore
    no SIFTS mapping. Asking for one raised four times per review and left four
    warnings on a review where nothing was wrong, which buries the warnings
    that mean something.
    """
    return bool(_PDB_RE.match(str(pdb_id).strip().lower()))


def check_pdb_id(pdb_id: str) -> str:
    pdb_id = str(pdb_id).strip().lower()
    if not _PDB_RE.match(pdb_id):
        raise ValueError(f"not a PDB id: {pdb_id!r}")
    return pdb_id


def check_chain(chain: str) -> str:
    chain = str(chain).strip()
    if not _CHAIN_RE.match(chain):
        raise ValueError(f"not a chain id: {chain!r}")
    return chain


def check_gene(gene: str) -> str:
    gene = str(gene).strip()
    if not _GENE_RE.match(gene):
        raise ValueError(f"not a gene symbol: {gene!r}")
    return gene


def check_taxon(taxon_id) -> int:
    try:
        value = int(taxon_id)
    except (TypeError, ValueError):
        raise ValueError(f"not a taxon id: {taxon_id!r}") from None
    if not 0 < value < 10**8:
        raise ValueError(f"not a taxon id: {taxon_id!r}")
    return value


# ------------------------------------------------------------------ manifest

def _read_manifest() -> dict:
    path = manifest_path()
    if not path.exists():
        return {"schema": 1, "entries": {}}
    with path.open() as fh:
        data = json.load(fh)
    data.setdefault("entries", {})
    return data


@contextlib.contextmanager
def _manifest_lock():
    cache_dir().mkdir(parents=True, exist_ok=True)
    with (cache_dir() / ".manifest.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def _record(rel: str, entry: dict) -> None:
    with _manifest_lock():
        manifest = _read_manifest()
        manifest["entries"][rel] = entry
        manifest["entries"] = dict(sorted(manifest["entries"].items()))
        tmp = manifest_path().with_suffix(".json.tmp")
        tmp.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
        tmp.replace(manifest_path())


def manifest() -> dict:
    """The parsed manifest (entries keyed by path relative to the cache directory)."""
    return _read_manifest()


# --------------------------------------------------------------------- fetch

_last_request: dict[str, float] = {}


def _today() -> str:
    return _dt.datetime.now().astimezone().strftime("%Y-%m-%d")


def _http_get(url: str) -> tuple[int, bytes, str]:
    """GET with a User-Agent, timeout, retries and backoff. Returns (status, body, content_type).

    404 is returned, not raised, so it can be recorded. Other failures raise RefDataError.
    """
    host = urllib.parse.urlsplit(url).hostname or ""
    if host not in ALLOWED_HOSTS:
        raise RefDataError(f"refusing to contact host outside the reference allow-list: {host}")
    last_error = None
    for attempt in range(RETRIES):
        wait = PACING_S - (time.monotonic() - _last_request.get(host, 0.0))
        if wait > 0:
            time.sleep(wait)
        _last_request[host] = time.monotonic()
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT_S) as resp:
                return resp.status, resp.read(), resp.headers.get("Content-Type", "")
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return 404, exc.read(), exc.headers.get("Content-Type", "")
            last_error = exc
            if exc.code not in (408, 425, 429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last_error = exc
        time.sleep(BACKOFF_S * (2 ** attempt))
    raise RefDataError(f"could not fetch {url}: {last_error}")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def cached_path(rel: str) -> Path:
    return cache_dir() / rel


def fetch_cached(rel: str, url: str, *, allow_missing: bool = False) -> Path | None:
    """Return the cached file for `rel`, fetching `url` on a miss.

    A remote 404 is recorded in the manifest (no file) so later calls know the record does
    not exist, even offline. With allow_missing a 404 returns None instead of raising.
    """
    path = cached_path(rel)
    if path.exists():
        return path
    entry = _read_manifest()["entries"].get(rel)
    if entry and entry.get("status") == 404:
        if allow_missing:
            return None
        raise RefDataNotFound(f"{url} is recorded as 404 in the manifest")
    if is_offline():
        raise RefDataOffline(
            f"{rel} is not in the reference cache and {OFFLINE_ENV} is set; "
            f"run scripts/fetch_reference.py with network access (source {url})"
        )
    status, body, ctype = _http_get(url)
    if status == 404:
        _record(rel, {"url": url, "retrieved": _today(), "sha256": None, "bytes": 0,
                      "status": 404, "file": None})
        if allow_missing:
            return None
        raise RefDataNotFound(f"{url} returned 404")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    tmp.write_bytes(body)
    tmp.replace(path)
    _record(rel, {"url": url, "retrieved": _today(), "sha256": _sha256(body),
                  "bytes": len(body), "status": status, "content_type": ctype, "file": rel})
    return path


def _json(path: Path):
    return json.loads(path.read_text())


# ------------------------------------------------------------------- UniProt

def _uniprot_url(acc: str, suffix: str) -> str:
    return f"https://rest.uniprot.org/uniprotkb/{acc}.{suffix}"


def get_uniprot(acc: str) -> dict:
    acc = check_accession(acc)
    return _json(fetch_cached(f"uniprot/{acc}.json", _uniprot_url(acc, "json")))


def get_uniprot_fasta(acc: str) -> str:
    acc = check_accession(acc)
    return fetch_cached(f"uniprot/{acc}.fasta", _uniprot_url(acc, "fasta")).read_text()


def fasta_sequence(fasta_text: str) -> str:
    return "".join(line.strip() for line in fasta_text.splitlines() if not line.startswith(">"))


def uniprot_summary(entry: dict) -> dict:
    """Accession, gene, organism, length and reviewed flag of a UniProt JSON entry."""
    genes = entry.get("genes") or []
    gene = genes[0].get("geneName", {}).get("value") if genes else None
    org = entry.get("organism", {})
    return {
        "accession": entry.get("primaryAccession"),
        "entry_name": entry.get("uniProtkbId"),
        "gene": gene,
        "organism": org.get("scientificName"),
        "taxon_id": org.get("taxonId"),
        "length": entry.get("sequence", {}).get("length"),
        "reviewed": str(entry.get("entryType", "")).startswith("UniProtKB reviewed"),
    }


def _parse_search_tsv(text: str) -> list[dict]:
    rows = []
    for rec in csv.DictReader(io.StringIO(text), delimiter="\t"):
        rows.append({
            "accession": rec["Entry"],
            "entry_name": rec["Entry Name"],
            "reviewed": rec["Reviewed"] == "reviewed",
            "gene_names": rec["Gene Names"].split(),
            "organism": rec["Organism"],
            "taxon_id": int(rec["Organism (ID)"]) if rec["Organism (ID)"] else None,
            "length": int(rec["Length"]),
            "protein_name": rec["Protein names"],
            "fragment": bool(rec["Fragment"].strip()),
            "sequence": rec["Sequence"],
        })
    return rows


def search_uniprot(gene: str | None = None, organism_id: int | None = None, *,
                   protein_name: str | None = None, size: int = 100) -> list[dict]:
    """Search UniProtKB by gene symbol and/or taxon (and optionally a protein-name word).

    The raw TSV answer is cached under uniprot/search/. Each row carries accession, entry
    name, reviewed flag, gene names, organism, taxon id, length, protein name, fragment flag
    and sequence.
    """
    if gene is None and organism_id is None:
        raise ValueError("search_uniprot needs a gene or an organism_id")
    parts, tag = [], []
    if gene is not None:
        gene = check_gene(gene)
        parts.append(f"gene_exact:{gene}")
        tag.append(f"gene-{gene}")
    if organism_id is not None:
        organism_id = check_taxon(organism_id)
        parts.append(f"organism_id:{organism_id}")
        tag.append(f"taxon-{organism_id}")
    if protein_name is not None:
        protein_name = check_gene(protein_name)
        parts.append(f"protein_name:{protein_name}")
        tag.append(f"name-{protein_name}")
    query = urllib.parse.quote(" AND ".join(parts), safe=":_")
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={query}"
           f"&fields={UNIPROT_SEARCH_FIELDS}&format=tsv&size={int(size)}")
    path = fetch_cached(f"uniprot/search/{'_'.join(tag)}.tsv", url)
    return _parse_search_tsv(path.read_text())


# ------------------------------------------------------------------ PDBe / SIFTS

def get_sifts(pdb_id: str, kind: str) -> dict:
    """PDBe/SIFTS residue-level mapping: uniprot, uniprot_segments, pfam, interpro or cath."""
    pdb_id = check_pdb_id(pdb_id)
    if kind not in SIFTS_KINDS:
        raise ValueError(f"kind must be one of {SIFTS_KINDS}, got {kind!r}")
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/{kind}/{pdb_id}"
    return _json(fetch_cached(f"pdbe/{pdb_id}/mappings_{kind}.json", url))[pdb_id]


def get_unobserved(pdb_id: str) -> dict:
    """Unobserved (unmodelled) residues of a PDB entry.

    Tries the PDBe unobserved_residues endpoint. When PDBe has no such record (it answers
    404 for some entries) falls back to the polymer_coverage endpoint, which lists the observed
    ranges per chain. Returns {"source": ..., "data": ...}.
    """
    pdb_id = check_pdb_id(pdb_id)
    base = "https://www.ebi.ac.uk/pdbe/api/pdb/entry"
    path = fetch_cached(f"pdbe/{pdb_id}/unobserved_residues.json",
                        f"{base}/unobserved_residues/{pdb_id}", allow_missing=True)
    if path is not None:
        return {"source": "unobserved_residues", "data": _json(path)[pdb_id]}
    path = fetch_cached(f"pdbe/{pdb_id}/polymer_coverage.json", f"{base}/polymer_coverage/{pdb_id}")
    return {"source": "polymer_coverage", "data": _json(path)[pdb_id]}


# ---------------------------------------------------------------- AlphaFold DB

def get_alphafold(acc: str, *, with_pae: bool = True) -> dict:
    """AlphaFold DB prediction for a UniProt accession.

    Returns {"metadata": dict, "model_path": Path, "pae": parsed json or None, "pae_path": Path or None}.
    The PAE json is returned only when the file is under PAE_MAX_BYTES; a larger file is
    still cached but pae and pae_path stay None.
    """
    acc = check_accession(acc)
    meta_path = fetch_cached(f"alphafold/{acc}/prediction.json",
                             f"https://alphafold.ebi.ac.uk/api/prediction/{acc}")
    metadata = _json(meta_path)[0]
    model_url = metadata["cifUrl"]
    model_name = model_url.rsplit("/", 1)[-1]
    model_path = fetch_cached(f"alphafold/{acc}/{model_name}", model_url)
    pae, pae_path = None, None
    if with_pae and metadata.get("paeDocUrl"):
        pae_name = metadata["paeDocUrl"].rsplit("/", 1)[-1]
        try:
            pae_path = fetch_cached(f"alphafold/{acc}/{pae_name}", metadata["paeDocUrl"],
                                    allow_missing=True)
        except RefDataOffline:
            pae_path = None
        if pae_path is not None:
            if pae_path.stat().st_size > PAE_MAX_BYTES:
                pae_path = None
            else:
                pae = _json(pae_path)
    return {"metadata": metadata, "model_path": model_path, "pae": pae, "pae_path": pae_path}


# --------------------------------------------------------------------- ConSurf-DB

_CONSURF_ROW = re.compile(
    r"^\s*(?P<pos>\d+)\s+(?P<seq>[A-Z-])\s+(?P<resname>[A-Z]{3}):(?P<resnum>-?\d+[A-Za-z]?):(?P<chain>\w+)\s+"
    r"(?P<score>-?\d+\.\d+)\s+(?P<grade>\d)\*?\s+"
)


def parse_consurf_grades(text: str) -> list[dict]:
    """Rows of a ConSurf-DB *_consurf_grades.txt file.

    insufficient is True where ConSurf marks the grade with * (too few aligned sequences or
    a wide confidence interval). Grades run 1 (variable) to 9 (conserved).
    """
    rows = []
    for line in text.splitlines():
        m = _CONSURF_ROW.match(line)
        if not m:
            continue
        rest = line[m.end():]
        grade_field = line[m.end("score"):m.end("grade") + 1]
        be = re.search(r"\s([be])\s", rest)
        rows.append({
            "pos": int(m["pos"]),
            "aa": m["seq"],
            "residue": m["resname"],
            "pdb_resnum": m["resnum"],
            "chain": m["chain"],
            "score": float(m["score"]),
            "grade": int(m["grade"]),
            "insufficient": "*" in grade_field,
            "buried": {"b": True, "e": False}.get(be.group(1)) if be else None,
        })
    return rows


def get_consurfdb(pdb_id: str, chain: str) -> list[dict] | None:
    """ConSurf-DB grades for a PDB chain, or None when unavailable (offline, 404, network error).

    Academic, non-commercial data: keep it an optional overlay.
    """
    pdb_id = check_pdb_id(pdb_id).upper()
    chain = check_chain(chain)
    name = f"{pdb_id}_{chain}_consurf_grades.txt"
    url = f"https://consurfdb.tau.ac.il/results/{pdb_id}/{chain}/{name}"
    try:
        path = fetch_cached(f"consurfdb/{name}", url, allow_missing=True)
    except RefDataError:
        return None
    return parse_consurf_grades(path.read_text()) if path else None


# ------------------------------------------------------------------- verification

def verify_cache() -> list[str]:
    """Problems found in the cache: missing files, size or sha256 mismatches. Empty when sound."""
    problems = []
    for rel, entry in _read_manifest()["entries"].items():
        if entry.get("status") == 404:
            continue
        path = cached_path(rel)
        if not path.exists():
            problems.append(f"{rel}: listed in manifest but missing")
            continue
        data = path.read_bytes()
        if len(data) != entry["bytes"]:
            problems.append(f"{rel}: size {len(data)} != manifest {entry['bytes']}")
        if _sha256(data) != entry["sha256"]:
            problems.append(f"{rel}: sha256 mismatch")
    return problems


# ------------------------------------------------------------------ ortholog panel

MIN_ORTHOLOG_IDENTITY = 0.70
PANEL_FILE = "panel.json"


def sequence_identity(query: str, reference: str) -> float:
    """Fraction of query residues identical to the reference in a global BLOSUM62 alignment."""
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"
    best = aligner.align(query, reference)[0]
    same = 0
    for (a0, a1), (b0, b1) in zip(*best.aligned):
        same += sum(1 for x, y in zip(query[a0:a1], reference[b0:b1]) if x == y)
    return same / len(query)


def rank_ortholog_candidates(rows: list[dict], reference_seq: str, paralog_seq: str | None = None, *,
                             min_identity: float = MIN_ORTHOLOG_IDENTITY, min_length: int | None = None,
                             max_length: int | None = None, gene: str | None = None) -> list[dict]:
    """Annotate UniProt search rows and order them best first.

    Eligible rows are non-fragment entries whose length lies between ``min_length`` and ``max_length`` (default 0.9 and 1.1 times the
    length of the reference) and whose sequence is at least ``min_identity`` identical to the reference and, when ``paralog_seq`` (the
    closest paralog of the reference) is given, closer to the reference than to that paralog: this separates the ortholog from the
    paralog where an entry has no gene name. Order: reviewed entries first, then length closest to the reference's, then identity to the
    reference, then accession so the choice is reproducible. ``gene``: the gene symbol of the target; ``gene_matches`` tells whether an
    entry carries it.
    """
    n = len(reference_seq)
    lo = int(0.9 * n) if min_length is None else min_length
    hi = int(1.1 * n) if max_length is None else max_length
    out = []
    for row in rows:
        if row["fragment"] or not lo <= row["length"] <= hi:
            continue
        to_ref = sequence_identity(row["sequence"], reference_seq)
        to_par = sequence_identity(row["sequence"], paralog_seq) if paralog_seq else None
        if to_ref < min_identity or (to_par is not None and to_ref <= to_par):
            continue
        item = {**{k: v for k, v in row.items() if k != "sequence"}, "identity_to_reference": round(to_ref, 4)}
        if to_par is not None:
            item["identity_to_paralog"] = round(to_par, 4)
        if gene:
            item["gene_matches"] = any(g.upper() == gene.upper() for g in row["gene_names"])
        out.append(item)
    out.sort(key=lambda r: (not r["reviewed"], abs(r["length"] - n), -r["identity_to_reference"], r["accession"]))
    return out


def get_panel() -> dict:
    """The specificity panel chosen by scripts/fetch_reference.py (data/reference/panel.json).

    {"target": {...}, "orthologs": [...], "paralogs": [...]}; each entry has accession, gene, organism,
    taxon_id, length, role ("default" or "optional") and the reason it was chosen.
    """
    path = cached_path(PANEL_FILE)
    if not path.exists():
        raise RefDataError(f"{path} is missing; run scripts/fetch_reference.py")
    return _json(path)
