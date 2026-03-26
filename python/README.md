# Python Components

Python tools in this directory provide orchestration around the C/C++ VSEARCH extension and Java RDP classifier.

## Commands

- `python/get_rdp_classifier.py`: Download latest RDP release and pretrained assets from SourceForge.
- `python/vsearch_plus/rdp_tav_taxonomy.py`: Compile/run Java native paired-end NB taxonomy assignment.
- Root wrappers are also available: `get_rdp_classifier.py` and `rdp_tav_taxonomy.py`.

Use `uv run` for execution when available.
