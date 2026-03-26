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

Entrypoint: `python3 rdp_tav_taxonomy.py`

1. Python wrapper resolves RDP jar/model paths from manifest.
2. Python wrapper compiles and invokes Java extension classes in `java/src/org/vsearchplus/rdp`.
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
- Python wrapper that builds/runs Java extension in-project

## Output semantics

- One TAV-level output stream (`--output`).
- Optional short-sequence list (`--shortseq-outfile`).
- No per-anchor output merge artifacts.
