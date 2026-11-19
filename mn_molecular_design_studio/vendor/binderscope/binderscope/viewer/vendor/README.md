# Vendored libraries

The viewer loads these files from disk (file://), so nothing is fetched at view time.

## 3Dmol.js 2.5.5

| item | value |
| --- | --- |
| file | `3Dmol-min.js` |
| version | 2.5.5 |
| source | https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.5/3Dmol-min.js |
| retrieved | 2026-10-01 |
| sha256 | `f7cc78921ae72e7623e89cdd111434f58c2efddd2ffda1cd212644b406fb8016` |
| size | 537792 bytes |
| licence | BSD-3-Clause, see `3Dmol-LICENSE.txt` |

The licence text comes from https://raw.githubusercontent.com/3dmol/3Dmol.js/2.5.5/LICENSE
(retrieved 2026-10-01, sha256 `4c6eaaed856f3f28a3b1a98e74f4a8a71618de7d51ea4155c29f6f793bcef861`). It also covers the code 3Dmol.js incorporates from
GLmol, Three.js and jQuery. Keep it next to the library in any bundle that ships the script.

No fonts are vendored. To re-fetch: `.venv/bin/python scripts/fetch_reference.py --force-vendor`.
