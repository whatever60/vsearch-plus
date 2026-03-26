# Java Components

This project reuses the official RDP Classifier Java runtime for taxonomy assignment.

- RDP executable jar is downloaded under `data/third_party/rdp_classifier/`.
- Taxonomy inference is orchestrated by Python tools in `python/vsearch_plus/`.
- Native paired-end NB extension sources live under `java/src/org/vsearchplus/rdp/`.

## Why Java Is Still First-Class Here

RDP's core classifier implementation is Java, and we call that runtime directly for inference.
The local Java extension keeps stock model internals and extends the query unit from single sequence to paired TAV.
