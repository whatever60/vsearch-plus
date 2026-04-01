# Data Directory

This directory is reserved for local runtime data that should not be bundled in the source distribution.

Expected contents:

- `data/rdp_classifier/`: downloaded pretrained data and `manifest.json`
- `data/manifests/`: small tracked notes about external data layout
- `../extern/java/rdp_classifier/`: downloaded RDP classifier code and jars

The first real `rdp-classifier` run downloads the required RDP assets automatically.
