# RFdiffusion β-pairing sources

`source/make_interface_tensor.py` is based on the MIT-licensed script distributed in the paper's [`interface_tensors.zip`](https://files.ipd.uw.edu/pub/strand_pairing_binders_2025/interface_tensors.zip). The paper's ALK2 example PDB is included for the input-page example.

Studio applies three correctness fixes to the upstream tensor script:

- include the final valid binder β-strand register;
- assign the full requested binder secondary-structure segment, including at the C terminus;
- map one-based target PDB residue IDs to the combined binder-plus-target tensor without a one-residue offset.

The input workflow copies one selected PDB chain, strips unrelated molecules, and renumbers the chain contiguously. It then maps the user's original residue selection to this job-local numbering before creating tensors. RFdiffusion receives the resulting full interface tensors through its standard `ScaffoldedSampler` SS/adjacency input.

The downloaded RFdiffusion and AF2-initial-guess repositories and the Zenodo release archive are kept in the sibling `tools-to-implement` directory for source review. Studio uses the already-installed `mn-rfdiffusion:latest` image and its existing AF2 Initial Guess workflow; AF2 weights and the original DeepMind-licensed AF2 source are not copied into this package.
