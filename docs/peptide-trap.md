# Peptide Trap

The **Peptide Trap** tile retargets published peptide-binding scaffolds from
Sahtoe et al., [Design of amyloidogenic peptide traps](https://doi.org/10.1038/s41589-024-01578-5).
It accepts a sequence, matches compatible peptide windows, threads each window
onto a scaffold, and redesigns the binder interface with PyRosetta.

## Use in Studio

1. Open **All tools → Peptide Trap** and select a project.
2. Use **Load C37 paper example** or enter a target sequence / single-record FASTA.
3. Select templates and inspect the complex in the structure preview. The binder
   is green and the target peptide is orange. Matching placements use inclusive,
   one-based sequence positions.
4. Choose **Generate peptide traps** for the design stage alone, or **Open full
   pipeline** for the editable starter:

   **Peptide Trap → AlphaFold2 Initial Guess → Rosetta interface analysis → Results**

Completed design results also offer **Continue with AF2 → interface analysis**.
Rosetta has already designed the binder sequence, so ProteinMPNN is optional.
The AF2 starter uses one multimer model, three recycles, a whole-complex initial
guess, no binder/interface templates, and single-sequence input. These settings
can be edited on the canvas. Prediction from an initial guess is a consistency
check; it does not establish binding.

## Available templates and assumptions

This adapter assumes that the target fragment binds as an extended β-strand.
It does not predict the unbound target's conformation or require a structure for
the complete target sequence.

| Template | Peptide length | Buried/exposed mask | Source model |
| --- | --- | --- | --- |
| C37 | 8 | `10101011` | `C37.pdb` |
| C104.1 model | 10 | `1110101011` | `C104.1.pdb` |

The deposited `seqthread.py` defines `0` as buried and `1` as exposed. Buried
positions accept exactly `GFAMILVWYST`; internal proline is excluded, following
that script. A compatible mask is an input screen, not a binding prediction.

The upstream script associates the ten-position C104 mask with `C104.1xray`,
whose deposited coordinates contain only eight peptide residues. Studio uses
the complete ten-residue **C104.1 design model** with this mask instead. This
substitution is explicit in the template label and output provenance; it is not
a reconstruction of missing crystal residues.

The source archive also includes other designs, including hairpin models.
They are preserved for inspection, but only the two templates above are enabled
for retargeting. Initial scaffold generation, hairpin retargeting, experimental
characterization, and the paper's complete selection procedure are outside this
adapter's scope.

## Runtime and installation

```bash
python scripts/peptide-trap/install.py
```

The installer downloads the authors' [Zenodo deposit](https://doi.org/10.5281/zenodo.10391229),
verifies its SHA256, and extracts 19 PDB models and six source/instruction files
under the configured reference root:

```text
peptide-traps/sahtoe-10391229/
  data.tar.gz
  source_manifest.json
  source/
```

It then builds `mn-peptide-trap:sahtoe-10391229` from the already installed
`mn-protein-scoring-pyrosetta:latest` image. Use `--skip-build` to install only
the reference assets, or `--reference-dir PATH` to override their location.
The Dockerfile is `scripts/peptide-trap/Dockerfile`; it adds the Studio adapter
to the existing PyRosetta installation. The authors' archived source remains
unchanged. The runtime receives a staged copy of the verified native protocol.

Source archive SHA256:

```text
72e228dda952b82a24c1c9c4403bcc2d4340d6a25f4f06cb148954352d9b67c3
```

The deposit is CC BY 4.0. PyRosetta retains its separate license. Attribution and
adaptation details are in `scripts/peptide-trap/PROVENANCE.md`.

## Protocol and outputs

The adapter retains the deposited `design_after_seqthread.xml` interface
FastDesign/FastRelax movers, beta scoring, interface restrictions, layer design,
composition constraints, and fixed target identities. Rosetta can adjust both
backbones and their relative placement during design and relaxation. Repeat
counts are configurable. Chain A is the peptide; chain B is the binder.

The original analysis/filter block is omitted, including dependencies on
PSIPRED and DAlphaBall. The starter uses Studio's separate AF2 and Rosetta
analysis implementations, whose metrics should not be treated as the paper's
original selection procedure.

The design stage checks target-sequence preservation, chain roles, finite
scores, and backbone displacement from the threaded template. The default
maximum RMSD is 2 Å for each partner and the whole complex after separate rigid
alignment. Rejected structures remain in raw artifacts; accepted structures
become normalized `generation.backbone_sequence` candidates with binder
sequences and peptide-only target structures. Output provenance records the
template, target window, mask, source DOI, and adapter version.

This is a CPU workflow. The CPU allocation limits container usage; the Rosetta
protocol runs each placement/replicate sequentially. AF2 uses the selected GPU
when the downstream pipeline is run.

## Local integration check

On 26 October 2026, the published C37 sequence `DVRFNFRE` completed threading,
interface redesign and relaxation in the dedicated container, producing one
accepted candidate. The check used seed 42, one design repeat, one relaxation
repeat and one CPU core. The target sequence was preserved; the complex
backbone RMSD was approximately 0.31 Å. Its continuation completed one AF2
multimer model with one recycle, followed by Rosetta interface analysis, through
the live pipeline API. This verifies the software handoff on the published
example; it does not validate a new target or reproduce the paper's full study.

The example results and a saved full starter pipeline are in the local
**Peptide Trap example** project. The saved starter retains the normal two
design/relaxation repeats and three AF2 recycles.
