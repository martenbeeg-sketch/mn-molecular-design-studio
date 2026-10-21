ARG BASE_IMAGE=mn-bindcraft2:v1.0.4-cu13
FROM ${BASE_IMAGE}

ARG BINDCRAFT2_SOURCE_REVISION=d502f19a187a1bc8c05e833ca9cbab7bc1934fe7

LABEL org.opencontainers.image.title="BindCraft2 with pairwise PAE export" \
      org.opencontainers.image.version="v1.0.4-pae-v1" \
      org.opencontainers.image.revision="$BINDCRAFT2_SOURCE_REVISION" \
      org.opencontainers.image.description="Studio-owned overlay that writes residue-ordered PAE JSON sidecars beside BindCraft2 CIF outputs"

COPY bindcraft2-pae.patch /tmp/bindcraft2-pae.patch
COPY apply_bindcraft2_pae_patch.sh /tmp/apply_bindcraft2_pae_patch.sh
RUN bash /tmp/apply_bindcraft2_pae_patch.sh /opt/bindcraft /tmp/bindcraft2-pae.patch \
    && python3 -m compileall -q /opt/bindcraft/bindcraft/af2.py /opt/bindcraft/bindcraft/MPNN_stage.py /opt/bindcraft/bindcraft/campaign.py \
    && rm /tmp/bindcraft2-pae.patch /tmp/apply_bindcraft2_pae_patch.sh

WORKDIR /work
