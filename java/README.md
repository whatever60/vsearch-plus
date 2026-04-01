# Java Components

This project reuses the official RDP Classifier Java runtime for taxonomy assignment.

- RDP executable jar is downloaded under `data/rdp_classifier/`.
- Taxonomy inference is launched through `./scripts/vsearch-plus-rdp-tav`.
- Native paired-end NB extension sources live under `java/src/main/java/org/vsearchplus/rdp/`.

## Why Java Is Still First-Class Here

RDP's core classifier implementation is Java, and we call that runtime directly for inference.
The local Java extension keeps stock model internals and extends the query unit from single sequence to paired TAV.
The Java bootstrap launcher also replaces the old Python taxonomy wrapper by resolving manifest paths, compiling the paired extension, and launching the classifier JVM.
