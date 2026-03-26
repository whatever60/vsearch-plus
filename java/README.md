# Java Components

This project reuses the official RDP Classifier Java runtime for taxonomy assignment.

- RDP executable jar is downloaded under `data/third_party/rdp_classifier/`.
- Taxonomy inference is orchestrated by Python tools in `python/vsearch_plus/`.
- No custom Java sources are currently required for the TAV taxonomy pipeline.

## Why Java Is Still First-Class Here

RDP's core classifier implementation is Java, and we call that runtime directly for inference.
The Python layer adds paired-end TAV semantics on top of the stock single-sequence classifier.
