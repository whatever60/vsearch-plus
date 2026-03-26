# Paired `usearch_global` Parity (Stock VSEARCH vs Extension)

This document captures parity status between stock `vsearch --usearch_global` and the paired-end extension path in this workspace.

## Scope and intent

- Command: `vsearch --usearch_global` with paired input via `--reverse`.
- Design target: preserve stock search backbone semantics, then extend only what is necessary to treat one read pair as one query unit.

## Current paired-mode contract

### Required paired inputs

- Query input: `--usearch_global` (R1) + `--reverse` (R2).
- Database input:
  - one interleaved paired FASTA/FASTQ file via `--db` (left record then right record), or
  - split paired FASTA/FASTQ files via `--db` (left) + `--db2` (right).
- Catalog TSV input is intentionally removed for paired `usearch_global` parity.

### Supported paired outputs

- `--alnout`
- `--userout`
- `--uc`
- `--blast6out`
- `--matched`, `--notmatched`
- `--matched` + `--matched2` (split paired outputs)
- `--notmatched` + `--notmatched2` (split paired outputs)
- `--dbmatched`, `--dbnotmatched`
- `--dbmatched` + `--dbmatched2` (split paired outputs)
- `--dbnotmatched` + `--dbnotmatched2` (split paired outputs)
- `--otutabout`, `--biomout`, `--mothur_shared_out`

### Not supported in paired mode (current)

- `--fastapairs`
- `--qsegout`
- `--tsegout`
- `--samout`
- `--lcaout`
- Custom `--userfields`

## Step-by-step parity map

### Step 1: Command dispatch and mode guards

Stock implementation: `cmd_usearch_global()` routes non-paired runs to `usearch_global()`. It requires `--db`, requires `--id` in `[0,1]`, and accepts the full stock output set (including SAM/segment/alignment-rich outputs).

Paired extension implementation: `cmd_usearch_global()` routes `--usearch_global` plus `--reverse` to `tav_usearch_global()`. It enforces:
- `--strand plus` only
- no `--fastapairs`, `--qsegout`, `--tsegout`, `--samout`, `--lcaout`
- no custom `--userfields`
- split-output pairing constraints (`--matched2` requires `--matched`, etc.)
- at least one paired-supported output
- `--db` and valid `--id`

Extension mechanism: paired mode is a dedicated execution path with explicit CLI constraints.

### Step 2: Database ingestion and target construction

Stock implementation: `search_prep()` loads DB from UDB or FASTA/FASTQ, applies optional DB masking (`dust`/`hardmask`), then builds `dbindex` for k-mer lookup (`dbindex_prepare`, `dbindex_addallsequences`).

Paired extension implementation: `tav_usearch_global()` loads paired targets from:
- interleaved `--db` (left/right alternating), or
- split `--db` (left) + `--db2` (right)
Catalog TSV input is rejected. Each DB pair is converted to one `TavRecord`:
- anchor length `anchor_len = min(left_len, right_len)` (or bounded by `--fastq_trunclen` via `get_anchor_len`)
- min/max length filtering with `--filter any|both`
- DB masking on both ends (`opt_dbmask`)
- header base normalized by removing trailing `/1` or `/2`
- abundance taken from left record (warning on mismatch)

Extension mechanism: one paired target record replaces one single-sequence target.

### Step 3: Query ingestion and paired normalization

Stock implementation: each query sequence is read once (or twice with reverse strand search when enabled), optionally query-masked, then searched.

Paired extension implementation: forward (`--usearch_global`) and reverse (`--reverse`) queries are read synchronously. For each pair:
- right read is reverse-complemented
- both ends are truncated to `anchor_len`
- query masking (`opt_qmask`) is applied on both ends
- abundance is `max(fwd_abundance, 1)`
- query unit becomes one paired `TavRecord`

Extension mechanism: query right-end orientation is normalized before matching, and search operates on paired units.

### Step 4: K-mer candidate scoring and top-hit preselection

Stock implementation: `search_onequery()` extracts unique query k-mers (`unique_count`), accumulates DB k-mer hit counts through `dbindex`, and keeps top candidates via min-heap (`search_topscores`).

Paired extension implementation:
- builds paired postings once (`left_postings`, `right_postings`) from DB k-mers
- for each query, computes `db_kmer_scores[i] = left_overlap_i + right_overlap_i`
- keeps candidates with `score >= min(opt_minwordmatches, qk_left + qk_right)`
- uses the same heap comparator behavior (higher k-mer count first, then shorter target length, then earlier index)

Extension mechanism: stock k-mer backbone is preserved with additive left+right evidence.

### Step 5: Unaligned filtering and delayed-alignment staging

Stock implementation: candidate loop applies `search_acceptable_unaligned()`, pushes survivors into delayed alignment batches (`MAXDELAYED`), and stops by `maxaccepts/maxrejects` limits.

Paired extension implementation: mirrors this control flow with:
- `paired_unaligned_filters_pass()` (paired analog of stock unaligned filters:
  `maxqsize`, `mintsize`, `minsizeratio`, `maxsizeratio`, `minqt`, `maxqt`, `minsl`, `maxsl`, `idprefix`, `idsuffix`, `self`, `selfid`)
- pending batch processed at `MAXDELAYED` or loop end
- effective limits:
  - `maxaccepts_effective`
  - `maxrejects_effective`
  - `maxhits_considered = maxaccepts + maxrejects - 1`

Extension mechanism: delayed-batch semantics are kept, but filters evaluate paired totals/per-end equivalents.

### Step 6: Alignment, aligned filters, and accept/weak/reject assignment

Stock implementation: delayed hits are globally aligned (`search16` SIMD batch, with linear-memory fallback), then `align_trim()` computes terminal/internal stats, `search_acceptable_aligned()` applies aligned thresholds, and hits are classified:
- accepted if `id >= 100 * opt_id`
- weak if aligned-filters pass but `id < 100 * opt_id`
- rejected otherwise

Paired extension implementation: each pending hit is aligned per end with `align_one_end_stock_style()` (linear-memory aligner), then aggregated:
- mismatches/gaps/columns/indels summed across ends
- identity metrics `id0..id4` computed and selected by `opt_iddef`
- `mid` and coverage metrics computed on combined totals
- aligned filter family applied in paired form:
  `weak_id`, `maxsubs`, `maxgaps`, `mincols`, `leftjust`, `rightjust`, `query_cov`, `target_cov`, `maxid`, `mid`, `maxdiffs`
- final accepted/weak/rejected decision uses the same `opt_id` cutoff semantics as stock

Extension mechanism: accept/weak/reject logic is stock-equivalent after two-end stat aggregation.

### Step 7: Hit retention, ordering, and report window

Stock implementation: `search_joinhits()` keeps only accepted/weak hits and sorts by stock comparator:
1. accepted before weak
2. aligned before non-aligned (if present)
3. higher identity
4. lower target index
Reporting applies `--maxhits`, `--top_hits_only`, and `--uc_allhits`.

Paired extension implementation: same retention and ordering semantics over `TavHit`, then:
- `toreport = min(opt_maxhits, report_hits.size())`
- `--top_hits_only` truncates below top identity
- UC output respects `--uc_allhits`

Extension mechanism: stock ranking semantics preserved on paired hit objects.

### Step 8: Per-query output emission and match accounting

Stock implementation: writes configured outputs (`alnout`, `userout`, `blast6`, `uc`, etc.); query is matched if any kept hit exists; `dbmatched[target]` increments for every accepted/weak hit by query abundance (`--sizein`) or by 1.

Paired extension implementation:
- matched query if `report_hits` non-empty; unmatched otherwise
- `dbmatched_abundance[target]` increments for every accepted/weak hit with stock-like `--sizein` weighting
- top hit only drives OTU assignment
- paired output renderers:
  - `--userout`: fixed schema (`query`, `target`, `id_left`, `id_right`, `id_total`, `d_left`, `d_right`, `d_total`)
  - `--blast6out`, `--uc`: stock-shaped paired adaptations
  - `--alnout`: paired summary plus separate R1/R2 alignments
  - `--matched`/`--notmatched`: interleaved paired FASTA or split with `--matched2`/`--notmatched2`
- if unmatched and table output enabled:
  - default: unassigned row behavior (`nullptr` target)
  - with `--unknown_name`: assign unmatched pairs to that explicit target name

Extension mechanism: stock accounting rules are preserved; formatting is paired-specific.

### Step 9: End-of-run table/database outputs

Stock implementation: adds zero-count rows for DB targets not observed in matches, writes OTU/BIOM/mothur tables, and writes `dbmatched`/`dbnotmatched` FASTA outputs from accumulated counts.

Paired extension implementation:
- adds all paired DB targets to OTU table with zero-count completion
- optionally adds `--unknown_name` row
- writes OTU/BIOM/mothur tables
- writes paired DB outputs as interleaved or split FASTA:
  - matched targets use accumulated `dbmatched_abundance`
  - non-matched targets use original DB abundance

Extension mechanism: terminal aggregation semantics are stock-aligned with paired output encoding.

## Remaining paired-only deltas

- Paired `--userout` currently emits a fixed paired summary schema:
  - `query`, `target`, `id_left`, `id_right`, `id_total`, `d_left`, `d_right`, `d_total`
- Alignment-rich single-sequence outputs (`fastapairs`, `qsegout`, `tsegout`, `samout`, `lcaout`) are intentionally unsupported in paired mode.

## Example paired command

```bash
vsearch --usearch_global query_r1.fa \
  --reverse query_r2.fa \
  --db tav_db_left.fa \
  --db2 tav_db_right.fa \
  --id 0.97 \
  --alnout paired.aln \
  --userout paired_userout.tsv \
  --uc paired.uc \
  --blast6out paired.b6 \
  --otutabout table.tsv
```
