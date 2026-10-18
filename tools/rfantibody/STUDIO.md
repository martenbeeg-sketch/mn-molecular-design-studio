# RFantibody in Molecular Design Studio

Studio vendors the MIT-licensed [RosettaCommons/RFantibody](https://github.com/RosettaCommons/RFantibody) source at commit `3e5a5bb1cb98befde9fc049ae3ff58cabb491f67`.

Build the ready-to-run CUDA 11.8 image from the repository root with:

```bash
scripts/build_rfantibody_image.sh
```

The image installs the pinned Python environment and downloads the RFdiffusion antibody, ProteinMPNN, and RF2 checkpoints from the locations documented by upstream. The optional AbMPNN checkpoint is supplied at runtime from Studio's configured `reference_files/mpnn/abmpnn.pt` file, so it is not copied into the image.

RFantibody's upstream `LICENSE` covers the source code. The model weights have their own terms; the optional AbMPNN checkpoint is distributed under CC BY 4.0. Review the upstream model and dataset attribution before redistributing the image or weights.
