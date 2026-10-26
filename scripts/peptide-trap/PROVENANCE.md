# Peptide-trap source provenance

Danny D. Sahtoe and colleagues, *Design of amyloidogenic peptide traps*,
Nature Chemical Biology (2024), DOI 10.1038/s41589-024-01578-5.

Authors' source/model deposit: https://doi.org/10.5281/zenodo.10391229
Archive: data.tar.gz, SHA256
72e228dda952b82a24c1c9c4403bcc2d4340d6a25f4f06cb148954352d9b67c3

The deposit declares Creative Commons Attribution 4.0 International:
https://creativecommons.org/licenses/by/4.0/
The separately installed PyRosetta distribution retains its own license.

Studio preserves the deposited source and models. Its adapter parameterizes
sequence/register matching and uses the authors' interface FastDesign/relax
settings with fixed target identities. External PSIPRED and DAlphaBall dependent
analysis is omitted from the runtime protocol; Studio's independent prediction
and analysis steps are separate. Input and output backbone geometry is checked.
This adaptation does not recreate the paper's initial scaffold generation,
experimental characterization, or entire selection procedure.

The deposited seqthread.py associates a ten-position burial mask with the
C104.1xray model, whose peptide coordinates have only eight residues. Studio
uses the complete ten-residue C104.1 design model with that mask instead.
C37 retains the directly matching eight-residue source model and mask.
