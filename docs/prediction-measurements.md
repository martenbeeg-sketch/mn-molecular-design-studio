# Co-fold prediction measurements

The machine-readable catalog is served by the existing
`GET /api/v1/molecules/output-vocabulary` endpoint under
`prediction_measurements`. It is built from the installed Studio predictor
registry plus the native AF2-initial-guess, ESMFold2, BindCraft2 and RFD4
measurement routes. Each entry distinguishes an adapter-supported measurement
from a value guaranteed to exist on every run.

## Installed co-folding profiles

| Installed engine | Saved PAE matrix | Scalar pTM/iPTM | Notes |
| --- | --- | --- | --- |
| AlphaFold 3 JAX | Conditional | Conditional | Native summary confidence and per-model confidence sidecars |
| AlphaFold 3 PyTorch | Conditional | Conditional | Per-model confidence sidecars |
| OpenFold3 JAX | Conditional | Conditional | Native summary confidence |
| Boltz-2 | Conditional | Conditional | Full PAE is retained by the installed adapter |
| Protenix V1 | Conditional | Conditional | Summary/full-data sidecars |
| Protenix V2 | Conditional | Conditional | Summary/full-data sidecars |
| RF3 Foundry | Conditional | Conditional | Compressed or plain confidence sidecars |
| ESMFold | No | No interface score | Structure B-factors give pLDDT; multichain use is exploratory |
| ESMFold2 Full | Conditional | Conditional | Per-candidate confidence arrays and scalar rows |
| ESMFold2 Fast | Conditional | Conditional | Per-candidate confidence arrays and scalar rows |
| Chai-1 | Conditional | Conditional | Model-specific score sidecars |
| AtlasFold-M | Conditional | Conditional | Model-specific score sidecars |
| OpenFold3 OpenBind-0 | Conditional | Conditional | Confidence sidecars |
| OpenDDE | Conditional | Conditional | Summary/full-data sidecars |
| ColabFold | Conditional | Conditional | Per-rank score/PAE sidecars |
| BoltzGen Fold | Conditional | Conditional | Selected refold sample is matched to its native archive |
| Mol-1 | Conditional | Conditional | Requires saved atom-confidence output; this is enabled by default |

“Conditional” means that the selected run and its exact matching sidecar must
contain the value. Monomer iPTM sentinels are not interface scores. The
ipSAE & Structure Analysis tile reports saved structure Cα B-factors as
B-factors; it does not relabel them as pLDDT. Predictor-native pLDDT summaries
remain separate, including for Mol-1.

## Common PAE fields and conversions

When a saved PAE matrix matches the protein Cα inventory and the candidate has
explicit binder and target chain roles, Studio records the raw Å means for
binder→target, target→binder and their bidirectional mean. Each can also be
shown after division by 31 Å. The raw values, directions, matrix path, selected
chains and calculation method are retained in the candidate provenance. A
missing matrix, ambiguous role or axis mismatch leaves these fields absent.

The `/31` operation gives PAE a common numeric scale. It does not calibrate
predictors or turn one aggregation into another. BindCraft2 `i_pAE` is already
a binder→target directional mean on that scale; multiplying by 31 recovers its
raw Å scale. RFD4 `oracle_ppi.ipae_coi` is a native chain-of-interest summary
that averages forward/reverse chain-pair PAE and then partner chains. COI means
the chain of interest, normally the designed chain. The candidate's
`folding_engine` identifies its configured RFD4 oracle backend.

`ipSAE`, `iPTM`, pLDDT, contacts and RMSD remain separate measurements. ipSAE
is a distinct PAE-based Analysis calculation; contacts and RMSD require saved
coordinates and, for pose RMSD, a suitable reference and residue mapping.

## ipSAE & Structure Analysis

The tile submits each candidate as an independent job to the same CPU-limited
ipSAE scoring container, using the existing batch worker pool. It reads the
candidate's predicted complex, matching PAE matrix, source design structure
when available, and configured hotspots. It reports heavy-atom and residue
contact counts, inter-chain Cα clash pairs, per-chain Cα B-factor means,
hotspot contact coverage and distance, and binder-only fitted Cα RMSD after
sequence matching. The RMSD is a binder shape comparison; it is not
target-relative pose RMSD. Missing references, hotspots, or PAE are explicitly
reported as unavailable or not applicable.

When the original design complex is available, the tile also maps binder and
target residues by sequence and compares inter-chain residue-contact pairs
under the same configured heavy-atom cutoff. Results include the number of
design contacts retained, the fraction of design contact pairs retained, and
contact-set Jaccard overlap. This is a residue-contact comparison; it does not
assert preservation of hydrogen bonds or other specific chemical interactions.
The mapped residue coverage is recorded so alignment gaps are visible.
