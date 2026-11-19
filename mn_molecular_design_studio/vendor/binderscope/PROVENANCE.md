# Integrated BinderScope analysis source

Source: https://github.com/schnickn/binderscope
Revision: c61ef9ca3dfc64dbde018ae7329d7538d8c583f1
Licence: MIT; see LICENSE, including the separate BindCraft licence boundary.

The analysis library, its profile, validation schemas/scripts and offline review
assets are included here as source, not loaded from another checkout. Studio's
adapter and scheduler integration live outside this directory. The upstream
launcher/server and app installation are not used. BindCraft itself is not
included. The review assets retain their bundled 3Dmol licence.

The CPU runtime builds from these repository files. Building needs the ordinary
Python dependency packages; running analysis requires neither a BinderScope
checkout, a BinderScope server, GitHub access nor a network connection.

Native BindCraft campaigns use the upstream readers and pipeline. Saved SQL
designs use a separate, explicitly labelled cohort adapter. They cannot establish
native campaign funnel, rejected/near-miss or stopping statistics. Author chain
and residue identities are retained alongside temporary canonical coordinates.
