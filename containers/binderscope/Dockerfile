FROM python:3.12-slim-bookworm

COPY mn_molecular_design_studio/vendor/binderscope /opt/binderscope
RUN pip install --no-cache-dir numpy==2.5.5 scipy==1.18.1 pandas==3.1.1 \
    scikit-learn==1.9.1 biopython==1.88 gemmi==0.7.5 tmtools==0.3.0 \
    Jinja2==3.1.6 openpyxl==3.1.5 jsonschema==4.26.0
COPY mn_molecular_design_studio/engine/workflow_engine/tools/binderscope.py /opt/studio/binderscope_runtime.py
COPY mn_molecular_design_studio/engine/workflow_engine/tools/binderscope_context.py /opt/studio/binderscope_context.py
ENV PYTHONPATH=/opt/binderscope/binderscope/src PYTHONDONTWRITEBYTECODE=1
LABEL org.opencontainers.image.source="https://github.com/schnickn/binderscope" \
      org.opencontainers.image.revision="c61ef9ca3dfc64dbde018ae7329d7538d8c583f1"
WORKDIR /work
