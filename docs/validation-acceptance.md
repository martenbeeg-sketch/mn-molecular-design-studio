# Validation & Acceptance

Validation & Acceptance is a standalone, saved-evidence-only Studio tile. It
applies configurable acceptance criteria to completed candidates' existing
model measurements and saved structures; it never launches a predictor. Its
shared policy evaluator is also used by graph acceptance filters.

The tile reads binder format and target input type from the source candidate
and campaign settings. BindCraft2 `binder` maps to Miniprotein; scaffold
modalities such as VHH, scFv, Fab, and ARP stay distinct. FASTA target paths map
to FASTA / IDR. Ambiguous records require researcher confirmation, and mixed
formats or target types must be run separately.

The starting criteria mirror BindCraft2 source revision `4a56313` core and modality
settings. Only criteria with a compatible measurement in the selected source
and run mode are active. The tile computes the same geometric
definitions as BindCraft2 for interface residues (binder residues with any
heavy atom within 4 Å of target) and interchain Cα clashes (≤2.5 Å). pTM and
iPTM are retained as source-model measurements. For BindCraft2 and RFD4 native
evidence, PAE can be displayed on the shared `PAE / 31 Å` numeric scale:
BindCraft2 `i_pAE` is a binder-to-target directional mean already normalized
by 31, while RFD4 `oracle_ppi.ipae_coi` is a bidirectional chain-of-interest
summary in Å. COI means the chain of interest selected by RFD4 (normally the
designed chain). The RFD4 candidate records the configured oracle engine; the
current Studio default is OpenFold3 JAX, while the upstream oracle keeps the
AF3-shaped metric contract. Dividing by 31 is an arithmetic scale conversion,
not calibration between predictors. Raw scores, direction and aggregation
remain visible and distinct.
FASTA / IDR targets do not
use target pLDDT or target-aligned pose RMSD. These defaults are editable and
recorded with every result.

BindCraft2 criteria without a compatible saved measurement stay inactive and
are shown as unavailable. In particular, a saved-evidence filter does not
calculate missing scores, ipSAE is not substituted for PAE, target-aligned
RMSD is not BindCraft2 Binder_RMSD, and scaffold-specific sequence/framework
filters need their scaffold annotations. Cyclic closure is informational unless
the researcher adds a threshold. BindCraft2 source i_pDAE is retained as a
source measurement and is the default ranking metric for BindCraft2 designs;
it does not gate acceptance.

The ipSAE Analysis tile also derives binder→target, target→binder, and
bidirectional contact PAE from the selected prediction's saved PAE matrix and
coordinates. It records each raw Å direction and the corresponding `/31` scale
conversion separately. The same tile reports a Studio recomputation of the
BindCraft2 `i_pDAE` definition: contact-masked anchored interface TM scores,
using the configurable Cα cutoff (default 8 Å) and the maximum across binder
and target anchors. This is not the BindCraft2 model's native score; it uses
the selected predictor's PAE, and no score is emitted if its matrix cannot be
mapped exactly to the saved protein Cα inventory. The computed score is
available as an explicit Validation & Acceptance criterion after ipSAE
Analysis. The renamed **ipSAE & Structure Analysis** tile also measures saved
complex heavy-atom and residue contacts, Cα clashes, hotspot contacts, and
binder-only fitted Cα RMSD when the corresponding reference data exists. When
the source design complex is available, it reports sequence-mapped preservation
of the designed binder–target residue-contact pairs. That result is not a
chemical interaction or hydrogen-bond preservation score. Its contact fraction
and Jaccard overlap can be chosen as separate acceptance measurements; missing
design-complex evidence remains unavailable. Structure Cα B-factors are
reported explicitly as B-factors and are not treated as pLDDT.

Results keep native pass status, independent engine measurements, configured
acceptance, missing evidence, and rejection reasons separate. A missing or
invalid value cannot pass a mandatory rule.
