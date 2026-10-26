ARG BASE_IMAGE=mn-protein-scoring-pyrosetta:latest
FROM ${BASE_IMAGE}
LABEL org.opencontainers.image.source="https://doi.org/10.5281/zenodo.10391229" \
      org.opencontainers.image.version="sahtoe-10391229-studio-v1"
COPY mn_molecular_design_studio/engine/workflow_engine/tools/peptide_trap_runtime.py /opt/studio/peptide_trap_runtime.py
COPY scripts/peptide-trap/PROVENANCE.md /opt/studio/PEPTIDE_TRAP_PROVENANCE.md
RUN chmod a+r /opt/studio/peptide_trap_runtime.py /opt/studio/PEPTIDE_TRAP_PROVENANCE.md && python -c "import pyrosetta; from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects; from pyrosetta.rosetta.protocols.relax import FastRelax"
ENV PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ENTRYPOINT ["python", "/opt/studio/peptide_trap_runtime.py"]
