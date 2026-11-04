# Protein Library and source structures

Protein Library contains reusable source proteins and protein-containing systems. A protein's sequence/composition is its identity; an experimental structure, a supplied template, and a prediction are distinct representations. A target is a role in a run or campaign, rather than a separate copy of a protein.

Studio retains declared PDB/mmCIF input files in the Database's content-addressed file store and records the originating run, input field, owner and project. These source representations are marked portable for the existing molecular Database synchronization and file-copy mechanism. Each owner/project has its own reference collection, so project access rules remain in effect. Sequence-only IDR targets are reusable Library records without an invented structural template. The ideal RFD4 sequence-carrier CIF is not an experimental reference structure.

Input collection runs when Studio tags a submitted job, and again through its completion registration hook. A failed prediction does not invalidate its supplied source protein. Explicit manual/disabled collection policies still apply. Analysis inputs and validation candidates are not automatically promoted to reusable source designs. Acceptance, evaluation and generated-result collection remain separate operations.

The packaged 8JYR and 7QF4 examples carry source-reference evidence only when their full submitted protein-chain sequences uniquely match the example, including copy counts. A changed residue invalidates that match. This records the PDB as a sequence source, and never claims its coordinates were supplied as a prediction template or assigns unspecified target/binder roles. Original PDB coordinates can have missing residues, glycans or a different composition: their own canonical record is linked as the sequence source instead of attaching unequal coordinates to the submitted sequence identity.

Database defaults include these reference collections. The details pane links a matched sequence system to its source PDB, and distinguishes source reference coordinates from the representation associated with a selected Evaluation. The Complex Analysis landing table also distinguishes a supplied target structure, an IDR sequence-only target, and a PDB-derived input system.

Saved Database records render while a background job-index refresh runs. An identical observation refresh preserves its evaluation projection; changed scientific inputs or scores still invalidate and rebuild that projection.
