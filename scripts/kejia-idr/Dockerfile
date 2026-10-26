# syntax=docker/dockerfile:1
ARG BASE_IMAGE=mn-protein-scoring-pyrosetta:latest
FROM ${BASE_IMAGE}
LABEL org.opencontainers.image.source="https://github.com/drhicks/Kejia_peptide_binders" \
      org.opencontainers.image.revision="7f65c18cc57b3e72d290c885c540276c7d9c6ab8"
ADD --checksum=sha256:2bbd4064994d46ff3ea51c1f347f0ab9781c06e797cbf94232051ec803cd4b9e https://raw.githubusercontent.com/drhicks/Kejia_peptide_binders/7f65c18cc57b3e72d290c885c540276c7d9c6ab8/threading/thread_peptide_sequence_new.py /opt/kejia/threading/thread_peptide_sequence_new.py
ADD --checksum=sha256:a2859b869854d877a4a04590f51374540fed551915ae9bf1d338020df13f9a41 https://raw.githubusercontent.com/drhicks/Kejia_peptide_binders/7f65c18cc57b3e72d290c885c540276c7d9c6ab8/threading/blosum62.py /opt/kejia/threading/blosum62.py
ADD https://raw.githubusercontent.com/drhicks/Kejia_peptide_binders/7f65c18cc57b3e72d290c885c540276c7d9c6ab8/LICENSE /opt/kejia/LICENSE
COPY mn_molecular_design_studio/engine/workflow_engine/tools/kejia_idr_threading.py /opt/studio/kejia_idr_threading.py
RUN chmod -R a+rX /opt/kejia /opt/studio && python -c "import pandas, pyrosetta; from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol"
ENV PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ENTRYPOINT ["python", "/opt/studio/kejia_idr_threading.py"]
