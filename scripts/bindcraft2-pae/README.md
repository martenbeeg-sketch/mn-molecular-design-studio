# BindCraft2 PAE output overlay

This patch keeps the complete directional AlphaFold PAE matrix from BindCraft2 validation predictions and writes it beside each refolded and ranked complex as `<structure-stem>.pae.json`. The matrix and per-residue pLDDT vector follow the matching CIF's chain and residue order; `chain_order` records the source and output chain mapping. Existing symmetric PAE metrics, filtering, and BindCraft2 ranking are unchanged.

The patch is maintained here in the Design Studio repository, outside the BindCraft2 checkout. To apply it to a source checkout for inspection or a future image build:

```bash
bash scripts/bindcraft2-pae/apply_bindcraft2_pae_patch.sh /path/to/BindCraft2
```

If an upstream update changes the patch context, application stops before editing files. Review and adapt the patch against that version, then rebuild the derived image.

Build the Studio image overlay from the existing CUDA image with:

```bash
bash scripts/build_bindcraft2_pae_image.sh
```

The default output image is `mn-bindcraft2:c0ee638-cu13-pae-v1`. Set `BINDCRAFT2_BASE_IMAGE` and/or `BINDCRAFT2_PAE_IMAGE` to override those tags.
