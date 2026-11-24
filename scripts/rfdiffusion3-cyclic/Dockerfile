FROM mn-foundry:cu128

ARG FOUNDRY_COMMIT
LABEL org.opencontainers.image.source="https://github.com/RosettaCommons/foundry" \
      org.opencontainers.image.revision="${FOUNDRY_COMMIT}"

COPY src/foundry /opt/studio-foundry/foundry
COPY models/rfd3/src/rfd3 /opt/studio-foundry/rfd3
COPY models/rfd3/configs /opt/studio-foundry/rfd3/configs
ENV PYTHONPATH=/opt/studio-foundry

RUN python -c "from rfd3.inference.input_parsing import DesignInputSpecification; assert 'cyclic_chains' in DesignInputSpecification.model_fields" && \
    rfd3 design --help >/dev/null
