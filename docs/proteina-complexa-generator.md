# Proteina-Complexa Generator

The **Proteina-Complexa** tile runs the complete installed native design pipeline. The separate **Proteina-Complexa Generator** tile runs only backbone and sequence generation, using the same target, ligand, motif, and sampling form. Its API workflow is `design.proteina_complexa_generate`.

Generator mode forces Single Pass and disables the native reward model. It does not run AF2/RF3 ranking, sequence redesign, or evaluation and needs only the selected generator and autoencoder checkpoints. Generation results have no confidence score or filter pass.

## Canvas use

- Start with the **Proteina-Complexa · Generator → MPNN → RF3** template, or add a Generator block. Click it to open the same task form used by the tile. A connected Target block is optional.
- Connect its output to ProteinMPNN / LigandMPNN, a compatible prediction node, or Results. Sequence redesign is optional because the generator already produces a sequence.
- To reuse a completed standalone generation, select its project and use **Generated project designs**, or click **Continue on pipeline canvas** in its results. Existing successful-design validation batches still require explicit filter passes.
- The default sequence node preserves ligand coordinates, bonds, and verified fixed motif residues. It rejects motif redesign if the native residue mapping is incomplete or does not preserve the motif residue identities.
- RF3 receives the exact ligand atom graph. Boltz-2 accepts the original standard CCD ligand code; custom ligand fragments require RF3. AlphaFold2 cannot predict ligand complexes. ESMFold predicts the binder alone. Foundry MPNN requires CIF input; use the shared ProteinMPNN / LigandMPNN node for these PDB generation outputs.

## Verification, 17 October 2026

- Full pytest: **318 passed**, one existing Starlette/httpx deprecation warning.
- All three task configs pass `complexa validate generate` in `mn-proteina-complexa:latest`, with only generator/autoencoder weights mounted. Their conditional input features also initialize successfully on CPU.
- The installed native motif matcher and sidecar exporter also pass a CPU self-match of the six prepared motif residues, preserving residue identities.
- Regression coverage includes generation-only scheduling, motif ordering through sorted JSON, project ownership and source eligibility, the shared tile/node forms, API preflight without a separate Target block, exact ligand bonds, fixed motif chain remapping, redesigned sequences in RF3 inputs, and unsupported ligand prediction routes.
- Both running local APIs expose the new workflow; the live Streamlit catalog and generator form were checked.
- CPU logs: `/tmp/mds-complexa-generator-preflight/{protein_binder,ligand_binder,motif_ligand}`. No new GPU jobs were submitted for the extraction.

Earlier full native pipeline smokes completed for protein (`20261017-160605-d9df8b98`), ligand (`20261017-161040-df0b7573`), and motif (`20261017-161854-a9edfea8`), with zero filter-passing designs. A motif input ordering mismatch was found and corrected during extraction: explicit residue order now travels with atom selections through sorted JSON. The earlier motif smoke predates this correction and does not validate the corrected motif mapping. The extracted generator and its downstream combinations have not had a fresh end-to-end GPU smoke.

The scheduler and paused V7 campaign were not changed.
