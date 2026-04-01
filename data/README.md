# Data Directory

This directory is reserved for local runtime data that should not be bundled in the source distribution.

Expected contents:

- `data/rdp_classifier/`: downloaded RDP runtime, pretrained data, and `manifest.json`
- `data/manifests/`: small tracked notes about external data layout

Typical setup:

```bash
./scripts/vsearch-plus-rdp-download --output-root data/rdp_classifier
```
