# Kejia IDR binder workflow

The **Kejia IDR Binder** tile runs template threading and flexible peptide docking
in `mn-kejia-idr:7f65c18`. Its PDB candidates can continue through the existing
Studio ProteinMPNN and AlphaFold2 Initial Guess workflows. The **Open full
pipeline** button creates an editable canvas draft with these three steps.
Existing threaded results also have a continuation button.

## Installed runtime and templates

Run `bash scripts/kejia-idr/install.sh` from the repository to install the pinned
upstream template collection under the configured reference directory and build
the dedicated image. An optional first argument overrides the reference root.
The image extends Studio's existing PyRosetta scoring image, which must already
be installed and licensed for the intended use. It includes the upstream
threading code and its local BLOSUM module. No new model weights are required by
the threading step; downstream ProteinMPNN and AlphaFold2 use their existing
Studio images and weight installations.

Source: [Kejia peptide binders, revision 7f65c18](https://github.com/drhicks/Kejia_peptide_binders/tree/7f65c18cc57b3e72d290c885c540276c7d9c6ab8).
The upstream source and template collection retain their original license.
The PyRosetta license remains separate. Source checksums are enforced when
building the threading image; jobs record the source revision and image tag.

## Inputs and handoff

Supply one sequence fragment or one FASTA record containing 4–40 standard
amino acids. Choose an uploaded PDB/mmCIF with exactly two protein chains, or an
installed template library. Assign the binder and peptide chains for uploads;
the job stages them as A and B with contiguous residue numbering. Original
uploaded files remain unchanged.

Walle and Walle2 are the authors' recommended starting libraries. Mini and
Strand primarily contain strand-pairing geometries. Library jobs use a
deterministic, sorted selection up to the requested maximum template count.
The template limit controls compute cost; it is not a template-quality ranking.

Native threading evaluates overlapping sequence/template registers. If the target
is longer than the template peptide, an output represents a subsequence of the
requested target. Every candidate records that range and the actual threaded
target sequence. Studio validates the output sequence before writing candidates.
Only binder chains are redesigned downstream; the threaded target is retained.

The dedicated adapter preserves native binary silent files and scorefiles,
exports each accepted pose as PDB, and writes a JSON handoff. It also enforces
the documented `max_farep` threshold before serializing a candidate; the pinned
upstream script does not consistently skip candidates exceeding that threshold.
Jobs with no surviving registers complete with an empty candidate collection,
allowing the pipeline to skip downstream work.

## Pipeline scope

The starter contains **Kejia threading → ProteinMPNN → AlphaFold2 Initial Guess
→ Results**. Its defaults follow the upstream initial pass: five sequences per
backbone, temperature 0.1, and five AF2 multimer models using a whole-complex
initial guess. It uses the Studio MPNN/AF2 implementations and candidate format,
so it should not be described as an exact reproduction of the paper's complete
pipeline. Predicted confidence scores do not establish experimental binding.

Iterative MPNN/AF2 cycles and additional filters can be added on the canvas.
The upstream sequence clustering, averaging of five AF2 predictions, optional
sequence-only confirmation, and modified RFdiffusion refinement are not added
automatically. The existing vanilla RFdiffusion image does not establish
compatibility with that modified branch.

## Example

**Load Tamarind example inputs** loads `FNCREWCWN` and the upstream
`templates/walle/conf_2.pdb`, matching the [Tamarind IDR page](https://app.tamarind.bio/tools/idr).
This template has a 223-residue binder and a 12-residue peptide. A reduced
integration check can use one flexible-docking repeat, one downstream sequence,
and one AF2 model; these settings check execution and handoffs, not design yield.

Verified locally on 2026-10-26 with that reduced example:

| Step | Completed output | Run ID |
| --- | --- | --- |
| Kejia threading, one repeat per register | Four threaded complexes | `20261026-150724-50003413` |
| Existing ProteinMPNN, one selected candidate | One redesigned binder | `20261026-151034-ceda44f9` |
| Existing AF2 multimer initial guess, one model and one recycle | One predicted complex | `20261026-151143-0394ff21` |

All stages retained the exact nine-residue target sequence and 223-residue
binder length. The AF2 result had binder pLDDT 95.82 and ipTM 0.338; its interface
confidence is low. This is an execution/handoff check, not a validated binder.
The three completed jobs are grouped in the **Kejia IDR example** Studio project.
The live API lists and validates the new workflow. The 72 focused checks cover
input validation, form submission and canvas handoff, candidate integrity,
pipeline contracts, and cluster capability reporting. The full canvas scheduler
was covered by contract tests; the reduced example's model stages were invoked
sequentially through the existing workflow functions.
