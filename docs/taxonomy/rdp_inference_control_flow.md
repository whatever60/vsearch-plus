# RDP Inference Control Flow for TAV Taxonomy

## Goal

Document the stock RDP classifier inference flow and how `vsearch_plus` reuses it while adding paired-end TAV semantics.

## Stock RDP Java Flow (Single-Sequence Backbone)

### 1. CLI entry dispatch

- Entry jar: `data/third_party/rdp_classifier/extracted/.../dist/classifier.jar`
- Main class: `edu/msu/cme/rdp/classifier/cli/ClassifierMain.java`
- `ClassifierMain.main()` dispatches `classify` to `edu/msu/cme/rdp/multicompare/Main.main()`.

### 2. Classify command parse and setup

- File: `edu/msu/cme/rdp/multicompare/Main.java`
- `Main.main()` parses `-o`, `-f`, `-g` or `-t`, `-w`, and sample files.
- It constructs `MultiClassifier` and calls `multiCompare(...)`.

### 3. Model and training data load

- File: `edu/msu/cme/rdp/classifier/utils/ClassifierFactory.java`
- `ClassifierFactory.getFactory(gene)` loads training resources from bundled `/data/classifier/...`.
- `ClassifierFactory.setDataProp(prop, false)` switches to external pretrained model files via property file.

### 4. Per-sequence inference loop

- File: `edu/msu/cme/rdp/multicompare/MultiClassifier.java`
- `multiCompare(...)` reads sequences from `MCSample.getNextSeq()` and runs:
  - `classifier.classify(new ClassifierSequence(seq), min_bootstrap_words)`
- File: `edu/msu/cme/rdp/classifier/Classifier.java`
- `Classifier.classify(...)` performs:
  - orientation check/reverse-complement when needed
  - word-index feature extraction
  - posterior scoring across taxa
  - bootstrap confidence estimation (100 trials)
  - rank assignment generation

### 5. Output formatting

- File: `edu/msu/cme/rdp/classifier/io/ClassificationResultFormatter.java`
- `getOutput(..., allRank/fixRank/filterbyconf/db/biom)` emits tab-delimited assignments.

## `vsearch_plus` TAV Extension Flow

### 1. Python entrypoints

- `rdp_tav_taxonomy.py` (root convenience wrapper)
- `python/rdp_tav_taxonomy.py` (Python CLI wrapper)
- `python/vsearch_plus/rdp_tav_taxonomy.py` (implementation)

### 2. Input normalization to paired TAV records

- `build_records_from_pair_fastas(...)`: reads `tav_left.fa` + `tav_right.fa`
- `build_records_from_catalog(...)`: reads `tav_denoised.tsv` style catalogs
- Each record is normalized to `(tav_id, abundance, left_seq_id, right_seq_id, left_seq, right_seq)`.

### 3. Reuse stock classifier unchanged

- `run_rdp(...)` executes stock jar command twice:
  - once on left anchors
  - once on right anchors
- This intentionally reuses stock Java inference exactly for each end.
- Pretrained models are reused via `-t .../rRNAClassifier.properties` discovered from `manifest.json`.

### 4. Pair-aware aggregation layer

- `parse_allrank_tsv(...)` parses per-end RDP allrank output.
- `aggregate_rank(...)` and `aggregate_lineage(...)` combine R1/R2 rank assignments with:
  - `--pair-filter both`: require confident agreement from both ends
  - `--pair-filter any`: allow one confident end; resolve conflicts by higher confidence
- Prefix-consistent truncation is enforced: once a rank is unresolved, deeper ranks are blank.

### 5. Output products

- `write_paired_table(...)` writes `PREFIX.paired.tsv` with side-by-side and paired fields.
- `write_rank_counts(...)` writes `PREFIX.rank_counts.tsv` with abundance-weighted rank totals.

## Bottom Line

The expensive taxonomy classifier remains the stock RDP Java implementation. The extension is intentionally a thin paired-end orchestration + aggregation layer on top.
