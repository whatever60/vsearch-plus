# UCHIME3-Denovo Parity Notes: Stock vs Native Paired

This document captures parity status between stock `vsearch --uchime3_denovo` and the native paired implementation in this workspace.

## Scope and intent

- Stock `uchime3_denovo` evaluates one denoised sequence stream.
- Paired native implementation evaluates synchronized denoised sequence pairs `(R1, R2)`.
- Design target: preserve stock command intent (abundance-aware de novo chimera detection), while extending to paired data with explicit breakpoint logic.

## Input and output semantics

### Stock mode: single-end `uchime3_denovo`

- Input: denoised FASTA file to `--uchime3_denovo`.
- Outputs: single-end chimera/nonchimera/report outputs.

### Paired mode: `uchime3_denovo` with paired input

- Input: two synchronized sequence files:
  - `--uchime3_denovo` for left reads
  - second positional input for right reads (`R2`)
- Alternative input mode: `--interleaved` with one interleaved paired FASTA/FASTQ stream
- Input preprocessing reuses paired loader behavior:
  - synchronized paired FASTX loading and record-count guards
  - qmask/hardmask preprocessing on both ends
- Current paired outputs:
  - `--tabbedout`: paired-extension report (`parent_a`, `parent_b`, `breakpoint_class`, scores, delta)
  - `--uchimeout`: stock-shaped UCHIME tabular output adapted to paired scoring
  - `--uchimealns`: paired chimera alignment summaries (left/right context)
  - `--borderline`: accepted for CLI parity (typically empty with uchime3 logic)
  - `--nonchimeras` + `--nonchimeras2`: paired FASTA outputs (split `R1` and `R2`)
  - `--chimeras` + `--chimeras2`: paired FASTA outputs (split `R1` and `R2`)
  - `--nonchimeras_tsv`: paired TSV catalog rows for non-chimeras
  - `--chimeras_tsv`: paired TSV catalog rows for chimeras

### Clarifications: catalog input and masking

- Catalog input note:
  - Stock `uchime3_denovo` does not accept TSV catalog input.
  - Paired `uchime3_denovo` also only accepts FASTA/FASTQ pairs (split positional `R1/R2` input or `--interleaved`).
  - `--chimeras_tsv` / `--nonchimeras_tsv` are paired output formats, not input formats.

- `qmask` / `hardmask` note:
  - `--qmask dust` runs DUST low-complexity masking on each end before candidate scoring/classification.
  - In this implementation, DUST scans 64 nt windows (stepping by 32 nt), computes a low-complexity score from repeated 3-mer composition, and masks windows whose score exceeds an internal threshold (`dust_level = 20`).
  - With soft masking (`--qmask dust` and no `--hardmask`), masked positions are lowercased.
  - With hard masking (`--hardmask` active), masked/lowercase positions are converted to `N`.
  - Masked positions are then excluded from downstream k-mer matching behavior (`MASK_SOFT` path treats lowercase as masked; hardmasked `N` is also excluded as ambiguous).

## Step-by-step parity map

### Step 1: Command dispatch and hard guards

Stock implementation: `cmd_chimera()` routes plain `--uchime3_denovo` (without paired input) to `chimera()`. Parameter guards enforce `abskew >= 1`, `xn > 1`, and `dn > 0`.

Paired native implementation: `cmd_chimera()` detects paired input (`R2` as second positional input or `--interleaved`) and routes to `uchime3_denovo_paired()`. It enforces the same `abskew/xn/dn` guards, enforces split FASTA pairing (`--chimeras` with `--chimeras2`, `--nonchimeras` with `--nonchimeras2`), and requires at least one output among paired FASTA, TSV catalog, `--tabbedout`, `--uchimeout`, `--uchimealns`, or `--borderline`.

Native mechanism: paired mode is an explicit alternate execution path, not a flag inside the stock `chimera()` loop.

### Step 2: Input loading, masking, and deterministic processing order

Stock implementation: `chimera()` reads the denoised FASTA from `--uchime3_denovo`, optionally applies dust/hardmask (`opt_qmask` path), and then calls `db_sortbyabundance()`. Query processing order is abundance descending.

Paired native implementation: `uchime3_denovo_paired()` loads paired FASTX input directly (split positional `R1/R2` or interleaved via `--interleaved`). The loader enforces synchronized record counts, applies the same qmask/hardmask operations to both ends, preserves separate R1 and R2 headers, requires the same first whitespace-delimited token on both ends, and then sorts records by:
1. abundance descending
2. header ascending
3. left sequence ascending
4. right sequence ascending

Native mechanism: each processing unit is a native `record_paired_s` pair `(R1,R2,abundance,header)` rather than a single sequence.

### Step 3: Parent search space construction and growth rule

Stock implementation: in denovo mode, `dbindex_prepare()` initializes an initially empty parent index. `chimera()` sets `opt_self=1`, `opt_selfid=1`, and `opt_maxsizeratio = 1/abskew` so candidate parents must be sufficiently more abundant. After each query is classified as non-chimeric (`status < suspicious`), `dbindex_addsequence(seqno, opt_qmask)` inserts that sequence for future queries.

Paired native implementation: `uchime3_denovo_paired()` now reuses the shared paired dbindex layer directly:
- `dbindex_prepare_paired(&records_paired, 1, opt_qmask)` initializes a paired k-mer index over the in-memory paired record store.
- the native engine keeps a logical paired target mapping through `target_seqnos_r1_paired` and `target_seqnos_r2_paired`.
- after each non-chimeric query, `dbindex_addsequence_paired(current_seqno, opt_qmask)` is called to incrementally expose that pair as future parent material.

Native mechanism: both modes use a monotonic non-chimeric parent index with incremental growth; paired mode keeps end-specific semantics inside the shared paired dbindex layer rather than by managing two separate stock dbindex insertions.

### Step 4: Candidate enumeration for the current query

Stock implementation: each query is split into `parts=4` (`partition_query()` for uchime/uchime2/uchime3). `search_onequery()` plus `search_joinhits()` collects accepted hits from each part, deduplicates targets, and builds `cand_list`. Full-query global alignments are then computed for every candidate (`search16`, with linear-memory fallback when SIMD overflows).

Paired native implementation:
- `search_onequery()` / `search_joinhits()` are not called directly in paired uchime mode.
- Query kmers are looked up through stock dbindex match lists (`dbindex_getmatchlist`/`dbindex_getmatchcount`/`dbindex_getmapping`).
- R1 k-mers only contribute from even mapped sequence slots (left-end parents), and R2 k-mers only contribute from odd mapped sequence slots (right-end parents), then both contributions are summed to a pair-level parent score.
- Candidate parents are scored by summed left/right k-mer overlap.
- The abundance gate (`p.abundance >= abskew * q.abundance`) is enforced.
- Top candidates are retained with stock-like bounded heap behavior (`maxaccepts + maxrejects` style cap).
- Full-query alignments are then computed for each retained candidate to build match vectors used by parent selection.

Native mechanism: this is now a paired analog of stock candidate discovery + full-query alignment materialization, rather than direct iteration over all eligible parents.

### Step 5: One-parent baseline (`best_one`)

Stock implementation: there is no explicit pre-parent-search `best_one` variable in `chimera.cc`. The closest stock counterpart is `QT` in `eval_parents()`, where:
- `QA` = identity(query, parent A)
- `QB` = identity(query, parent B)
- `QT = max(QA, QB)`
This one-parent baseline exists only after two parents have already been selected/evaluated. If two parents are not found (`find_best_parents() == 0`), stock returns `Status::no_parents` and emits no-parent outputs (placeholder parent fields in `uchimeout`).

Paired native implementation: `best_one_score` is now computed from the selected two-parent evaluation itself as:
`best_one = max(match_QA, match_QB)` (count scale), mirroring stock's `QT = max(QA,QB)` baseline semantics within the same evaluation stage.
For no-parent cases, `best_one` falls back to full query-length count for reporting continuity.

Native mechanism: the paired path now aligns with stock by deriving one-parent baseline from the evaluated parent pair, rather than running a separate pre-model one-parent optimization sweep.

### Step 6: Two-parent preselection (stock-like 32 bp voting)

Stock implementation detail (`find_matches()` + `find_best_parents()`):
- Build `match[candidate][qpos]` from the full-query alignment CIGAR:
  - `1` if aligned query/parent symbols overlap at `qpos`
  - `0` otherwise
- For each candidate, compute `smooth[candidate][qpos]` = number of matches in the trailing 32 bp window ending at `qpos`.
- For each `qpos`, compute `maxsmooth[qpos]` = best window score among candidates not yet selected.
- "Win" definition: a candidate wins position `qpos` if `smooth[candidate][qpos] == maxsmooth[qpos]` and `maxsmooth[qpos] > 0`.
- Parent 1 = candidate with maximum total wins.
- Before picking parent 2, for every `qpos` where parent 1 was tied for best, zero all candidates' `match` entries inside that winning 32 bp window (stock's way to remove already-explained windows).
- Recompute smoothing/wins on the masked matrix and pick parent 2.

Paired native implementation detail:
- Use concatenated query axis:
  `[left query positions in forward order][right query positions in forward order]`.
- Build candidate match vectors from full-query alignment CIGARs (paired analog of stock `find_matches()`).
- Run the same 32 bp smoothing + per-position winner counting + winner-window wipeout to choose parent 1 then parent 2.
- Parent 1/2 winner selection now calls the same low-level helper (`select_best_two_parents_from_match_matrix`) from both stock `find_best_parents()` and paired Step 6.
- No fallback all-pairs sweep is performed when no two parents are found; behavior is stock-like no-parent handling.

Low-level relationship in current code:

- Stock Step 6 call stack:
  ```text
  find_matches(...)
    -> find_best_parents(...)
    -> select_best_two_parents_from_match_matrix(...)
  ```
- Paired Step 6 call stack:
  ```text
  compute_match_vector(...)
    -> parse_cigar_operations_from_string(...)
    -> parse_cigar_string(...)
    -> select_best_parent_pair(...)
    -> select_best_two_parents_from_match_matrix(...)
  ```

Native mechanism: paired parent preselection now uses alignment-derived match vectors and shared stock winner-window selection logic rather than direct character-equality heuristics or all-pairs fallback.

### Step 7: Diff coding and `h` optimization for two-parent models

Stock implementation: `eval_parents()` builds aligned strings for query and two parents, marks ignored columns (gap columns, neighbors of gaps, and ambiguous symbols), then emits `diffs` symbols (`A`, `B`, `N`, `?`, or space). It scans breakpoint position `i` across the alignment and tracks left/right counts:
- `Y` (yes) votes
- `N` (no) votes
- `A` (abstain/other)
Here, `Y/N/A` are vote roles (support, contradiction, abstain), not nucleotide letters.
It computes:
`left_h = left_y / (xn*(left_n + dn) + left_a)`
`right_h = right_y / (xn*(right_n + dn) + right_a)`
`h = left_h * right_h`
and also checks the reverse orientation (swap yes/no) when anti-chimera orientation is stronger. Best `h` determines the model.

Paired native implementation:
- For each selected parent pair, both ends are globally aligned and converted into unified gapped strings (`Q`, `A`, `B`) using stock-like insertion normalization.
- Ambiguous and gap-adjacent columns are marked ignored exactly in stock style.
- `diffs` are emitted with stock symbols (`A/B/N/?/space`), with reverse-orientation swap handling (`A<->B`) when anti-chimera orientation is stronger.
- `h` optimization scans one concatenated paired axis with right-end scan direction reversed (left: start->end, right: end->start), yielding one `best_i`.
- Paired breakpoint class is then derived from `best_i` location:
  - in left end -> `LEFT_BREAK`
  - at end boundary -> `MIDDLE_BREAK`
  - in right end -> `RIGHT_BREAK`

Low-level relationship in current code:

- Paired Step 7 call stack:
  ```text
  evaluate_parent_pair(...)
    -> build_end_alignments(...)
    -> parse_cigar_operations_from_string(...)
    -> parse_cigar_string(...)
    -> stock-style diff/vote/model scoring
  ```

Native mechanism: the paired path now uses stock-style diff coding and one-pass breakpoint scan; paired class labels are metadata derived from the chosen breakpoint position.

### Step 8: Final chimeric decision rule

Stock implementation (`uchime2/3` path in `eval_parents()`): after best model selection, define:
- `cols`: number of non-ignored alignment columns
- `match_QM`: count of columns where query matches the selected two-parent model
- `QA`: `100 * match_QA / cols` (query vs parent A identity)
- `QB`: `100 * match_QB / cols` (query vs parent B identity)
- `QT = max(QA, QB)` (best one-parent identity)
- `QM = 100 * match_QM / cols` (two-parent model identity)
`uchime3_denovo` marks chimeric when the model is perfect and best single parent is imperfect:
`match_QM == cols` and `QT < 100.0`.

Paired native implementation: after selecting/evaluating two parents on the concatenated paired axis, it uses the same uchime3 condition in paired form:
- `match_QM`: non-ignored columns explained by the selected two-parent model
- `cols`: total non-ignored columns
- `QT = max(QA, QB)`: best one-parent identity
Classify as chimeric iff:
- `match_QM == cols`
- `QT < 100.0`

Native mechanism: parity is kept in the decision rule itself; paired breakpoint class (`LEFT_BREAK`/`MIDDLE_BREAK`/`RIGHT_BREAK`) is retained as annotation in paired reports, not as an extra classification gate.

### Step 9: Parent-index update after classification

Stock implementation: all non-chimeric outcomes (`status < suspicious`) are added to the denovo parent dbindex, enabling incremental parent search for later, lower-abundance queries.

Paired native implementation: all non-chimeric pairs are added to the paired denovo parent index via `dbindex_addsequence_paired(...)`; chimeric pairs are never added.

Native mechanism: identical growth policy, but parent units are paired records.

### Step 10: Output materialization

Stock implementation file/schema summary:
- `--chimeras`, `--nonchimeras`: single-end FASTA (one record per query sequence).
- `--uchimeout`: stock UCHIME tabular schema.
  - default form includes columns equivalent to:
    `h, query, parentA, parentB, top, QM, QA, QB, AB, QT, LY, LN, LA, RY, RN, RA, divdiff, YN`
  - with `--uchimeout5`, the `top` column is omitted (stock compatibility mode).
- `--uchimealns`: stock textual alignment report with aligned query/parents plus `Diffs`, `Votes`, and `Model` lines.
- `--borderline`: FASTA output for borderline calls (normally unused in uchime2/3 decision path).

Paired native implementation comparison against stock:
- Same option name, mostly same schema family:
  - `--uchimeout` keeps stock-like column semantics (including `--uchimeout5` behavior), but values come from paired scoring/counts.
- Same option name, paired-adapted body format:
  - `--uchimealns` now includes paired alignment blocks with stock-style `Diffs/Votes/Model` lines split by side (`*_LEFT`, `*_RIGHT`), while keeping paired context labels.
- Same option name, different FASTA layout:
  - `--chimeras`/`--chimeras2` and `--nonchimeras`/`--nonchimeras2` are split paired FASTA outputs (`R1` file and `R2` file).
- Same option name, effectively empty in paired uchime3:
  - `--borderline` accepted for CLI parity but not populated by current paired decision path.
- Extra paired-only outputs not present in stock:
  - `--tabbedout` paired report with 7 columns:
    `query_tav, parent_a, parent_b, breakpoint_class, best_one_score, best_two_score, delta`
  - `--chimeras_tsv` / `--nonchimeras_tsv` paired catalog TSV with 5 columns:
    `tav_id, abundance, header, left_anchor, right_anchor`

Native mechanism: paired outputs carry more metrics than stock because the paired algorithm has extra latent variables (two ends + breakpoint class) that do not exist in stock's single-axis model.

## Operational constraints in paired mode

- Paired native mode is entered when paired input is provided (`R2` as second positional input or `--interleaved`).
- At least one output must be specified among:
  - `--nonchimeras` + `--nonchimeras2`
  - `--chimeras` + `--chimeras2`
  - `--tabbedout`
  - `--uchimeout`
  - `--uchimealns`
  - `--borderline`
  - `--chimeras_tsv`
  - `--nonchimeras_tsv`

## Parity status summary

- Current status: paired `uchime3_denovo` follows paired FASTX ingestion routing (split positional input or `--interleaved`) and removes catalog-mode input.
- Current status: paired mode now mirrors key stock denovo behavior for parent-pool growth (nonchimera-only), index-backed candidate discovery, 32bp winner-vote parent preselection, and stock-style one-pass `h` breakpoint optimization.
- Current algorithm status: paired chimera core now follows stock diff/vote/model mechanics on a concatenated paired axis; paired breakpoint classes are retained as derived annotations.
- Main open parity work items:
  - calibrate/validate paired classification behavior against broader real datasets
  - decide whether to further tighten candidate discovery to an even closer `search_onequery()`/`search_joinhits()` analogue
  - continue harmonizing paired `--uchimealns` body formatting with stock text layout where practical

## Example paired command

```bash
./bin/vsearch \
  --uchime3_denovo denoised_r1.fasta denoised_r2.fasta \
  --tabbedout chim_report.tsv \
  --uchimeout chim_uchime.tsv \
  --uchimealns chim_alns.txt \
  --nonchimeras_tsv nonchim_pairs.tsv \
  --chimeras_tsv chim_pairs.tsv \
  --nonchimeras nonchim_pairs.fa \
  --nonchimeras2 nonchim_pairs_r2.fa \
  --chimeras chim_pairs.fa \
  --chimeras2 chim_pairs_r2.fa
```
