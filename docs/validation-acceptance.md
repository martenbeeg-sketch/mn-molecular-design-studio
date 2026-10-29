# Validation & Acceptance

Validation & Acceptance is a standalone Studio tool tile for independently
re-predicting completed, filter-passing binder designs. It uses the existing
optimized-prediction workflow, candidate identity/provenance, scheduler, and
Results path. Its shared policy evaluator is also used by graph acceptance
filters; there is no second prediction runner or scheduler.

The tile reads binder format and target input type from the source candidate
and campaign settings. BindCraft2 `binder` maps to Miniprotein; scaffold
modalities such as VHH, scFv, Fab, and ARP stay distinct. FASTA target paths map
to FASTA / IDR. Ambiguous records require researcher confirmation, and mixed
formats or target types must be run separately.

The starting criteria mirror BindCraft2 source revision `4a56313` core and modality
settings. Only criteria with a compatible measurement in the selected engine
and run mode are active. The tile currently computes the same geometric
definitions as BindCraft2 for interface residues (binder residues with any
heavy atom within 4 Å of target) and interchain Cα clashes (≤2.5 Å). pTM and
iPTM are retained as engine-produced measurements; FASTA / IDR targets also
offer target-structure Cα B-factor pLDDT. These defaults are editable and are
recorded with every result.

BindCraft2 criteria without a compatible measurement stay inactive and are
shown as unavailable. In particular, a complex-only re-prediction does not
measure unbound binder pLDDT, i_pAE is not substituted with ipSAE, target-aligned
RMSD is not BindCraft2 Binder_RMSD, and scaffold-specific sequence/framework
filters need their scaffold annotations. Cyclic closure is informational unless
the researcher adds a threshold. BindCraft2 source i_pDAE is retained as a
source measurement and is the default ranking metric for BindCraft2 designs;
it does not gate acceptance.

Results keep native pass status, independent engine measurements, configured
acceptance, missing evidence, and rejection reasons separate. A missing or
invalid value cannot pass a mandatory rule.
