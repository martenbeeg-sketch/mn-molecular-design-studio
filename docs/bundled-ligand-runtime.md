# Bundled ligand runtime

Studio includes the ligand runtime code it uses at
`mn_molecular_design_studio/engine/workflow_engine/vendor/ligand_runtime`. This
runtime is part of the Studio source and package. Ligand workflow execution must
not look up a sibling `ovo-ligand` checkout or require `MN_LIGAND_ROOT`.

The vendored snapshot is based on `ovo-ligand` commit
`f87c9d3bcd9c7e81aef1161d341db8651949585f`. It contains the shared job and
artifact modules, ligand chemistry and MD services, workflow implementations,
manifests, and the OpenMM MM/GBSA helper used by the bundled workflows. The
separate ligand app UI is not included. Workflow changes needed by Studio should
be made in this bundled copy; upstream syncs should be deliberate and tested.

The Studio adapters mount this directory read-only in the matching engine
containers. Their `PYTHONPATH` points to the same mounted directory, which
contains the bundled `mn_ligand` Python package. The complex-preparation adapter
also uses this directory as the source root for its host-side helper and nested
MODELLER/PDBFixer jobs.

This removes the source-checkout dependency. Scientific execution still needs
the configured container images and external runtime assets such as MODELLER,
force-field files, and model/reference data.
