# Optional Chaperone interoperability

The Candidate collection node can use Chaperone's real `PipelineDataHandler` and `FilterBlock` through `mn_molecular_design_studio.chaperone_adapter`. Chaperone remains an optional execution dependency; the Studio app and its scheduler do not import it into the main Python environment.

Install the reproducible Python 3.12 Linux worker once with `scripts/setup_chaperone_worker.sh`. It verifies the locked Chaperone source revision, installs `requirements/chaperone-worker-py312-linux.lock` into a versioned environment under the Studio app home, and writes `chaperone-worker.json` there. The canvas then discovers the worker automatically; students do not configure Python paths. `CHAPERONE_REPO`, `CHAPERONE_SETUP_PYTHON`, and `MDS_CHAPERONE_WORKER_DIR` can be set by an administrator when the checkout or interpreter is in a different location. Environment overrides (`MDS_CHAPERONE_PYTHON`, `MDS_CHAPERONE_SRC`, `MDS_CHAPERONE_BLOCKS`) remain available for isolated testing. Studio calls the adapter through a JSON subprocess boundary, so the main app keeps its own Python runtime. The selected adapter, worker path and source revision, and Chaperone handler/schema/filter source files contribute to the graph cache key.

The adapter creates one Chaperone row per Studio candidate ID, carries the parent ID, model and chain provenance, and maps only explicitly known chain sequences. It stores the complete Studio candidate as a sidecar and adds the complete Chaperone output row as a receipt; unsupported fields therefore remain inspectable. Separate model predictions remain separate rows. An optional Chaperone filter expression is auxiliary evidence: all rows, including rejected and unresolved rows, return to Studio unchanged, and Studio's own acceptance and review rules remain authoritative.

The isolated interoperability tests use the real Chaperone package and local structural fixtures. The adapter calls Chaperone's `PipelineDataHandler`, `make_row`, `make_sequence_row`, and `get_sequences_from_row`; its optional metric check calls the native `FilterBlock`. Studio remains the scheduler and record authority. It does not run Chaperone's `ExecutionEngine` or make Chaperone a competing workflow runtime.

The optional compact-storage path uses Chaperone's `PipelineDataHandler.store()` to create CatCIF entries and `PipelineDataHandler.write_parquet()` for handler records. `read_catcif_cif_bytes()` materializes one referenced entry for downstream Studio consumers. `complex_analysis_handoff` recognizes those entries in structure history and materializes them in the analysis workspace without modifying the archive. The Studio sidecar carries fields outside Chaperone's native schema through Parquet restoration. PDB, mmCIF, JSONL Results and content-addressed Studio artifacts remain supported; no historical files or Database records are migrated.

CatCIF is optional and requires the configured Chaperone interpreter to include Chaperone's CatCIF dependency (`catcif_tools`). The main Studio Python environment is not modified. If the dependency or worker paths are unavailable, the optional conversion is unavailable rather than silently falling back to lossy conversion.

Focused interoperability and storage tests:

```sh
PYTHONPATH="$PWD:/path/to/chaperone/src:/path/to/chaperone/src/blocks" \
  /path/to/python3.12 -m pytest -q tests/test_chaperone_interop.py tests/test_chaperone_storage_interop.py
```

`requirements/chaperone-worker-py312-linux.lock` contains the exact resolved package versions used for the isolated interoperability tests. The Chaperone source checkout is separately pinned in its header; the setup script fails closed on a different checkout instead of silently changing schema behavior.

The storage fixture checks multiple CatCIF entries, separate model and parent identities, coordinate recovery, Parquet/sidecar restoration, a retrieved structure entering the existing Biotite analysis adapter, and unchanged original input files.
