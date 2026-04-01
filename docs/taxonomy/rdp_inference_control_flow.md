# RDP Inference Control Flow for Native TAV Classification

## Goal

Document the stock RDP inference path and where paired-end TAV logic extends it.

## Stock RDP path (single-sequence)

1. `ClassifierMain.main()` dispatches to classify command path.
2. `Main.main()` parses options and builds `MultiClassifier`.
3. `ClassifierFactory` loads training resources / properties.
4. `Classifier.classify(...)` computes:
   - orientation handling
   - k-mer word feature extraction
   - genus posterior scoring
   - bootstrap confidence over query words
5. `ClassificationResultFormatter` renders output.

## Native paired TAV extension path

Entrypoint: `./scripts/rdp-classifier`

1. Java bootstrap launcher resolves RDP jar/model paths from `data/rdp_classifier/manifest.json`.
2. Java bootstrap launcher compiles and invokes the paired classifier classes in `java/src/main/java/org/vsearchplus/rdp`.
   - Launcher accepts stock-style classify aliases for core options (`-o/-f/-c/-w/-s/-t/-g`) plus paired input args.
   - Launcher supports stock-style `-b/--bootstrap_outfile` and `-h/--hier_outfile`.
   - Metadata/biom stock options are rejected explicitly in paired mode.
3. Java `PairedClassifierMain` reads paired records from:
   - split inputs (`--input`, `--input2`) or
   - interleaved input (`--interleaved`).
4. Java `PairedNaiveBayesClassifier.classifyPair(...)` performs native paired NB:
   - orientation normalization on each anchor using stock `TrainingInfo.isSeqReversed`
   - word-probability row construction for R1 and R2 using stock model parameters
   - full paired posterior scoring over combined rows
   - bootstrap over combined R1+R2 word pool
   - lineage confidence accumulation and stock result object creation
5. Output is written via stock `ClassificationResultFormatter`.

## Reuse versus new code

Reused directly:

- `ClassifierFactory`
- `TrainingInfo` and trained probability structures
- `ClassifierSequence`
- `ClassificationResult` and `RankAssignment`
- `ClassificationResultFormatter`

New extension code:

- `PairedNaiveBayesClassifier` (paired likelihood/bootstrap core)
- `PairedClassifierMain` (paired IO and CLI bridge)
- `RdpTavTaxonomyMain` (Java launcher that resolves assets and compiles/runs the paired extension)

## Output semantics

- One TAV-level output stream (`--output`).
- Optional short-sequence list (`--shortseq-outfile`).
- No per-anchor output merge artifacts.
