# Generator examples

Offline example coordinates and task configurations used by the standalone tool
forms and pipeline node forms. Loading an example makes a private input copy;
these originals are never edited by a job.

## Sources

- `1CRN.pdb`, `1IEP.pdb`: [RCSB PDB](https://www.rcsb.org/). Used for the public
  examples on Tamarind's PXDesign and Genie3 pages. Downloads and SHA256 hashes
  are recorded in `download-provenance.json`.
- `pxdesign/`: [PXDesign examples](https://github.com/bytedance/PXDesign/tree/main/examples),
  installed source revision `f788441`. The PD-L1 input and quick-start YAML are
  included; MSA files and model weights are not packaged.
- `rfdiffusion/`: [RFdiffusion examples](https://github.com/RosettaCommons/RFdiffusion/tree/main/examples).
  Input PDBs, native example scripts, TIM barrel scaffolds and target fold
  conditioning tensors are included. The installed image does not identify its
  exact code revision. Public example 5TPN motif selections follow Tamarind.
- `genie3/`: [Genie3 design data](https://github.com/aqlaboratory/genie3/tree/main/data/design),
  installed source revision `d77ae5a`. BinderBench problem files and target
  coordinates, and all 30 MotifBench problem/motif pairs are included. Large
  target MSAs are omitted because generation does not require them.
- `protpardelle/`: [Protpardelle-1C examples](https://github.com/ProteinDesignLab/protpardelle-1c/tree/main/examples),
  installed source revision `ee37840`. Native sampling YAML and matching motif
  structures include MotifBench, BindCraft targets and nanobody examples.

Upstream licenses are retained in `source-licenses/`. Structures retain their
original attribution/header records. Checkpoints remain in the configured
common reference directory, outside the application package.

The public Tamarind Protpardelle-1C URL did not expose a tool form during review;
that tool uses the same Studio visual layout with upstream task configurations.
