# TAV Extension Implementation Log

## 2026-03-24 - Initial integration into VSEARCH source tree

### Scope implemented in this pass
- Added a new VSEARCH extension module in `vsearch/src/tav_extension.cc` and `vsearch/src/tav_extension.h`.
- Wired paired-end-native TAV logic into existing command flow without replacing existing single-end behavior.
- Reused VSEARCH FASTX parsing (`fastx_open`, `fastx_next`), reverse-complement (`reverse_complement`), and FASTA writing (`fasta_print_general`) primitives.

### Implemented command routing
- `--fastx_uniques`:
  - Existing single-end path unchanged.
  - If `--reverse` is provided, route to paired-end TAV dereplication (`tav_fastx_uniques`).
- `--cluster_unoise`:
  - Existing VSEARCH path unchanged.
  - If `--reverse` is provided, route to paired-end TAV denoising (`tav_cluster_unoise`).
- `--uchime3_denovo`:
  - Existing VSEARCH path unchanged.
  - If `--reverse` is provided, route to paired-end TAV chimera detection (`tav_uchime3_denovo`).
- `--usearch_global`:
  - Existing VSEARCH path unchanged.
  - If `--reverse` is provided, route to paired-end TAV search (`tav_usearch_global`).

### Option validity updates
- Added support for `--reverse` and `--fastaout_rev` with `--fastx_uniques` in paired mode.
- Added support for `--reverse` with `--usearch_global` in paired mode.
- Added support for `--tabbedout` with `--cluster_unoise` and `--uchime3_denovo` for TAV TSV outputs.

### TAV data format used
- Catalog TSV header:
  - `tav_id\tabundance\theader\tleft_anchor\tright_anchor`

### Algorithmic notes (v1)
- Pairing is ordered and R2 is reverse-complemented at ingestion.
- Anchor length uses `--fastq_trunclen` when specified (>0), otherwise full overlap of available R1/R2 lengths.
- Pair distance is scalar sum of per-anchor mismatch/length differences.
- UNOISE-like assignment rule uses:
  - `abund(Q) / abund(C) <= 2^(-1 - alpha * D_pair)` with `alpha = --unoise_alpha`.
- UCHIME3-like paired mode evaluates middle/left/right breakpoint classes and compares best one-parent vs best two-parent explanations.
- Paired UCHIME3 parent preselection includes a stock-like 32bp window voting stage over combined left/right match vectors when possible.
- Paired breakpoint/class selection now uses stock-like vote-derived `h` scoring (with score tie-break) per candidate parent pair.
- Global assignment in paired mode uses per-end and total identity checks based on `--id`.

### Known limitations in this pass
- This pass is scalar and correctness-first; no SIMD/multithreading optimization.
- Paired table/sample parsing is intentionally simple (`sample=` header tag if present, otherwise `sample_unknown`, or `--sample` override).
- Paired UCHIME3 classification currently uses a simplified paired analog (selected two-parent model fully reconstructs query pair, while neither selected parent alone is perfect) rather than full stock UCHIME3 internal scoring.

### Next steps queued
- Build and fix compile/runtime issues.
- Add a reproducible smoke test with synthetic paired FASTQ inputs.
- Run end-to-end on user-provided S3 data (subsample if needed).
- Inspect outputs and calibrate thresholds/options for biological plausibility.

## 2026-03-24 - Build, smoke tests, and real-data subsample

### Build status
- VSEARCH builds successfully with the new TAV extension code integrated into `bin/vsearch`.

### Synthetic validation status
- Built a paired synthetic dataset and executed all requested command stages end-to-end:
  - `--fastq_filter` (paired)
  - `--fastx_uniques` (paired TAV mode via `--reverse`)
  - `--cluster_unoise` (paired TAV mode via `--reverse`)
  - `--uchime3_denovo` (paired TAV mode via `--reverse`)
  - `--usearch_global` (paired TAV mode via `--reverse`)
- Output behavior was coherent:
  - deterministic TAV IDs
  - exact-pair derep behaved as expected
  - assignment table matched known synthetic sample composition

### Real-data run details
- Input downloaded from S3:
  - `merged_isolate.fq.gz`
  - `merged_isolate_R2.fq.gz`
- Subsampled to first 10,000 read pairs for fast iteration.
- Initial filter settings (`--fastq_trunclen 120`, `--fastq_minlen 120`) retained zero reads.
  - Diagnostic: all reads in this subsample are length 110.
  - Decision: adjusted to `--fastq_trunclen 100`, `--fastq_minlen 100`, `--fastq_maxee 5.0`.

### Real-data output summary (10k read-pair subsample)
- `tav_catalog.tsv`: 711 lines (710 TAV records + header)
- `tav_denoised.tsv`: 190 lines (189 centroids + header)
- `tav_nonchim.tsv`: 83 lines (82 nonchimeras + header)
- `tav_chim.tsv`: 108 lines (107 chimeras + header)
- `query_to_tav.tsv`: 9854 lines (9853 assigned queries + header)
- `tav_table.tsv`: 83 lines (82 TAV rows + header)

### Biological/algorithmic observations from real-data subsample
- One dominant TAV is highly abundant after denoising (expected in isolate-heavy samples).
- Chimera reporting produced a mix of `LEFT_BREAK`, `RIGHT_BREAK`, and `MIDDLE_BREAK` classes, indicating all three paired-end-native chimera modes were exercised.
- Assignment rate was high on this subsample (`9853 / 9973` filtered queries assigned), suggesting representative coverage is reasonable for current thresholds.

### Key implementation nuance discovered
- Existing VSEARCH `--fastq_filter` already applies pair-aware keep/discard semantics when `--reverse` is provided (pair retained only if both ends pass). This phase-1 requirement is therefore satisfied by extending usage and routing around that existing behavior rather than reimplementing filter internals.
