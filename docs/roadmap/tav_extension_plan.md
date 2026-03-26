# TAV-native extension plan for UNOISE3, UCHIME3-denovo, and usearch_global

## 1. Scope and explicit constraints

This document captures the design for a **paired-end-native** extension of the ASV pipeline around the concept of a **Terminal Anchor Variant (TAV)**.

### 1.1 Problem statement

For long amplicons where paired-end overlap is insufficient for read merging, define the fundamental denoised variant unit by the ordered pair of terminal anchors rather than by one fully observed merged amplicon sequence.

### 1.2 Working definition

A **Terminal Anchor Variant (TAV)** is an ordered pair of denoised terminal anchor sequences:

- `left_anchor`: the forward-orientation anchor derived from R1
- `right_anchor`: the forward-orientation anchor derived from R2 after reverse-complementing into amplicon orientation

A TAV is therefore:

```text
TAV = (left_anchor, right_anchor)
```

with optional metadata such as abundance, internal numeric ID, taxonomy, etc.

### 1.3 Requirements from the project brief

The implementation should follow these rules.

1. **Paired-end reads are first-class citizens.**
   The algorithms should not be written as a single-end method plus an afterthought pairing layer.

2. **This should be a new project, not a fork.**
   The project structure and code style may mimic VSEARCH, and useful low-level functions may be reused, but the new command implementations should live in **new source files**.

3. **Core algorithm complexity should not be artificially reduced.**
   CLI, logging, and secondary I/O may be simplified, but the algorithmic extension itself should be faithful to the spirit of the original methods and should accept that paired-end-native logic increases complexity.

4. **The three commands of interest are the new core deliverables.**
   - TAV-native UNOISE3-like denoising
   - TAV-native UCHIME3-denovo-like chimera detection
   - TAV-native usearch_global-like read assignment and count-table generation

5. **Output model changes are expected.**
   A TAV catalog is not one FASTA of merged sequences. The primary sequence outputs become:
   - one FASTA for left-anchor representatives
   - one FASTA for right-anchor representatives
   - one catalog/manifest mapping IDs to paired anchors

6. **Initial implementation can ignore speed-focused engineering.**
   SIMD reuse is optional for the first pass. A clean, correct scalar implementation is acceptable.

### 1.4 Non-goals for v1

These are explicitly out of scope for the first implementation pass.

- Full CLI parity with VSEARCH
- Full logging parity with VSEARCH
- UDB support
- Bit-exact reproduction of every VSEARCH output format
- Production-level multithreading
- SIMD acceleration
- Exact reuse of every VSEARCH internal struct

The first version should prioritize **correct algorithmic structure**.

---

## 2. Design principles

### 2.1 Preserve the original command meanings as much as possible

The paired-end extension should preserve the command-level semantics:

- **UNOISE3-like step:** denoise observations into representative variants
- **UCHIME3-denovo-like step:** remove chimeras using abundance-aware de novo logic
- **usearch_global-like step:** assign reads to representatives and build tables

The change is that the atomic sequence object is no longer a single observed sequence, but an **ordered pair of observed terminal anchors**.

### 2.2 Reuse single-end primitives, not single-end command structure

The safest reusable pieces are the low-level operations that still make sense per anchor:

- FASTA / FASTQ pair parsing
- reverse-complement logic
- header/sample parsing
- global alignment
- edit distance / mismatch counting
- dereplication helpers
- table writing helpers

The parts that should become TAV-native are the command-level algorithms:

- candidate generation
- distance aggregation across anchors
- chimera model structure
- centroid / parent selection
- representative output model
- count-matrix row semantics

### 2.3 Keep the pair ordered and oriented

Everywhere in the system:

- `left_anchor` always means the amplicon-left observed anchor
- `right_anchor` always means the amplicon-right observed anchor in forward amplicon orientation

That means R2 is reverse-complemented during ingestion, not later.

### 2.4 Use a catalog as a first-class artifact

Even if two FASTA files are emitted, a TAV is fundamentally one paired object. Therefore every command should also be able to read/write a catalog such as:

```text
TAV_ID    left_header    right_header    left_seq    right_seq    abundance
```

This catalog becomes the canonical representation. The paired FASTA files are convenient secondary views.

---

## 3. Core data model

### 3.1 Observation-level record

Each input read pair is converted into one normalized observation:

```text
Observation = (
    query_id,
    sample_id,
    left_anchor,
    right_anchor,
    abundance = 1
)
```

### 3.2 Dereplicated pair record

After exact paired dereplication:

```text
PairUnique = (
    left_anchor,
    right_anchor,
    abundance,
    member_count,
    representative_header
)
```

### 3.3 Representative TAV record

After denoising:

```text
TAV = (
    tav_id,
    left_anchor,
    right_anchor,
    abundance,
    source = denoised,
    status = provisional | nonchimera
)
```

### 3.4 Assignment record

For the search/counting step:

```text
Assignment = (
    query_id,
    tav_id,
    left_score,
    right_score,
    total_score,
    left_distance,
    right_distance,
    total_distance
)
```

---

## 4. Repository structure and important files

The intention is to create new algorithm files while reusing existing utility-like code where helpful.

## 4.1 Proposed new files

```text
src/
  tav_types.h
  tav_pairio.h
  tav_pairio.cc
  tav_catalog.h
  tav_catalog.cc
  tav_distance.h
  tav_distance.cc
  tav_align.h
  tav_align.cc
  tav_unoise3.h
  tav_unoise3.cc
  tav_chimera3.h
  tav_chimera3.cc
  tav_global.h
  tav_global.cc
  tav_table.h
  tav_table.cc
  tav_results.h
  tav_results.cc
```

### 4.1.1 What each file should own

- `tav_types.*`
  - shared structs
  - observation / TAV / assignment types
  - enum types for statuses and match classes

- `tav_pairio.*`
  - synchronized paired FASTQ / FASTA reading
  - orientation normalization
  - anchor extraction
  - sample extraction from headers

- `tav_catalog.*`
  - canonical catalog reading/writing
  - left/right FASTA export and import
  - ID allocation

- `tav_distance.*`
  - pair-level distance and score aggregation
  - scalar edit distance / mismatch counting wrappers
  - prefix/suffix scoring helpers for chimera logic

- `tav_align.*`
  - alignment wrappers
  - optional reuse of existing global aligner later
  - no SIMD requirement in v1

- `tav_unoise3.*`
  - exact pair dereplication
  - denoising loop
  - centroid selection / sequence-to-TAV mapping output

- `tav_chimera3.*`
  - abundance-sorted de novo chimera detection over TAVs
  - parent-pair evaluation
  - breakpoint-case evaluation

- `tav_global.*`
  - paired query to paired representative assignment
  - best-hit selection
  - optional per-query intermediate match output

- `tav_table.*`
  - count matrix generation
  - OTU-style table output
  - BIOM export

- `tav_results.*`
  - text outputs for mapping/match/alignment summaries

## 4.2 Existing code that is worth reusing conceptually

Even though the new project is not a fork, these existing VSEARCH areas are the most useful templates:

- pair parsing and pair processing patterns
- FASTA/FASTQ IO and header/sample parsing conventions
- global alignment utilities
- output table concepts (`otutabout`, BIOM)
- existing split of responsibilities among clustering, chimera detection, and search

## 4.3 Minimal touch points if integrating into a VSEARCH-like shell later

If this were eventually wired into a VSEARCH-like dispatcher, the minimal top-level additions would be:

- new command declarations in the central command registry
- new option parsing for:
  - paired FASTQ inputs
  - paired FASTA catalog inputs
  - anchor length
  - pair-aware thresholds
- calls from the dispatcher into:
  - `tav_unoise3()`
  - `tav_chimera3_denovo()`
  - `tav_usearch_global()`

The core logic should still remain in the new files above.

---

## 5. Overall execution model for the TAV pipeline

The whole paired-end-native pipeline should be:

```text
paired FASTQ
  -> normalize / orient / extract anchors
  -> exact paired dereplication
  -> TAV-UNOISE3 denoising
  -> TAV-UCHIME3-denovo chimera removal
  -> paired query vs TAV assignment
  -> count matrix / BIOM
```

### 5.1 Shared preprocessing contract

All three commands should use the exact same preprocessing rules:

1. read synchronized pair
2. quality-filter each end as desired
3. orient R2 into forward amplicon space
4. cut fixed-length anchors
5. discard pairs if either end cannot produce the required anchor
6. preserve sample labels

That ensures that the representatives and the queries live in the same sequence space.

---

## 6. Symbol glossary used in the command plans

To avoid ambiguity, use the following symbols consistently.

### 6.1 Pair-level objects

- `Q` = query TAV or query paired observation
- `Q_L` = query left anchor
- `Q_R` = query right anchor

- `C` = candidate centroid / representative / TAV
- `C_L` = candidate left anchor
- `C_R` = candidate right anchor

- `A`, `B` = candidate parent TAVs for chimera detection
- `A_L`, `A_R` = left/right anchors of parent A
- `B_L`, `B_R` = left/right anchors of parent B

### 6.2 Scores and distances

- `d(X, Y)` = distance between two single anchors
- `s(X, Y)` = score between two single anchors
- `D_pair(Q, C)` = paired distance between query `Q` and candidate `C`
- `S_pair(Q, C)` = paired score between query `Q` and candidate `C`

### 6.3 Prefix/suffix helper scores for chimera logic

- `s_pref(X, Y, i)` = score of prefix-to-prefix comparison up to split position `i`
- `s_suf(X, Y, i)` = score of suffix-to-suffix comparison after split position `i`

### 6.4 Abundance terms

- `abund(X)` = abundance of object `X`
- `alpha` = UNOISE alpha parameter
- `gamma` = de novo chimera abundance skew threshold (`abskew`)

---

## 7. Command-specific plan: UNOISE3 -> TAV-UNOISE3

## 7.1 Command goal

Denoise paired anchor observations into representative TAVs.

## 7.2 High-level command name

Suggested internal command name:

```text
tav_unoise3
```

This name is deliberately descriptive even if the eventual public CLI chooses another spelling.

## 7.3 Side-by-side: original vs extension

| Aspect | Original single-sequence logic | TAV-native extension |
|---|---|---|
| Atomic object | one observed sequence | ordered pair `(left_anchor, right_anchor)` |
| Exact dereplication key | sequence | exact anchor pair |
| Representative | one denoised sequence | one denoised paired-anchor representative |
| Distance evaluation | one sequence vs one centroid | left-anchor distance + right-anchor distance |
| Abundance rule | one sequence abundance ratio | pair abundance ratio |
| Output representative | one FASTA | left FASTA + right FASTA + catalog |
| Optional map output | query -> centroid | query pair -> TAV |

## 7.4 Original VSEARCH logic to preserve conceptually

Conceptually preserve these properties of `cluster_unoise`:

1. abundance-aware denoising
2. clustering depends jointly on abundance ratio and sequence distance
3. chimera removal is **not** part of this command
4. output can include representative sequences and optional mapping-like outputs

## 7.5 TAV-native observation model

Each normalized paired observation becomes:

```text
Q = (Q_L, Q_R)
```

where:

- `Q_L` is the left observed anchor
- `Q_R` is the right observed anchor

After exact paired dereplication, the denoising unit is the unique pair `(Q_L, Q_R)` with abundance `abund(Q)`.

## 7.6 Pair-level distance definition

Define the pair distance as:

```text
D_pair(Q, C) = d(Q_L, C_L) + d(Q_R, C_R)
```

where:

- `d(Q_L, C_L)` = single-anchor distance between the query left anchor and candidate left anchor
- `d(Q_R, C_R)` = single-anchor distance between the query right anchor and candidate right anchor

For v1, `d()` may be implemented as one of the following, in order of simplicity:

1. Hamming distance if anchor lengths are fixed and indels are disallowed
2. edit distance
3. global-alignment-derived mismatch count

### 7.6.1 Optional balanced-distance guard

To avoid pathological assignments where almost all discrepancy is concentrated in one anchor, add an optional per-end guard:

```text
max(d(Q_L, C_L), d(Q_R, C_R)) <= d_end_max
```

This is not a simplification; it is a paired-end-specific consistency constraint.

## 7.7 Pair-level abundance rule

Preserve the UNOISE-style abundance-vs-distance structure:

```text
abund(Q) / abund(C) <= 2^(-1 - alpha * D_pair(Q, C))
```

where:

- `abund(Q)` = abundance of query pair `Q`
- `abund(C)` = abundance of candidate centroid `C`
- `alpha` = denoising alpha parameter
- `D_pair(Q, C)` = total paired distance defined above

This is the direct paired-end generalization of the original one-sequence rule.

## 7.8 Denoising loop

### 7.8.1 Input order

Sort dereplicated exact pairs by decreasing abundance.

### 7.8.2 Main loop

For each dereplicated pair `Q` in abundance order:

1. compare `Q` to existing TAV centroids `C`
2. compute `D_pair(Q, C)`
3. check paired abundance-vs-distance acceptance rule
4. if at least one centroid is acceptable:
   - assign `Q` to the best acceptable centroid
5. otherwise:
   - create a new TAV centroid from `Q`

### 7.8.3 Best acceptable centroid selection

A conservative default tie-break order is:

1. smallest `D_pair(Q, C)`
2. then smallest `max(d(Q_L, C_L), d(Q_R, C_R))`
3. then highest `abund(C)`
4. then lexical ID stability

## 7.9 Output artifacts

### 7.9.1 Mandatory

- `tav_left.fa`
- `tav_right.fa`
- `tav_catalog.tsv`

### 7.9.2 Optional intermediate mapping output

- `query_to_tav.tsv`

Suggested fields:

```text
query_id    tav_id    d_left    d_right    d_total    abundance_query    abundance_tav
```

### 7.9.3 Optional alignment output

- `query_to_tav_left.aln`
- `query_to_tav_right.aln`

## 7.10 Important new files for this command

- `tav_unoise3.cc`
- `tav_distance.cc`
- `tav_pairio.cc`
- `tav_catalog.cc`
- optionally `tav_align.cc`

## 7.11 Minimal reuse from existing single-end logic

The minimal reusable conceptual pieces are:

- scalar per-anchor alignment/distance routine
- FASTA/FASTQ parsing and header handling
- abundance extraction and sample parsing conventions

What should **not** be preserved structurally:

- one-record-in, one-record-out object model
- representative-as-single-sequence assumption

---

## 8. Command-specific plan: UCHIME3-denovo -> TAV-UCHIME3-denovo

## 8.1 Command goal

Remove chimeric TAVs from abundance-sorted denoised TAV representatives.

## 8.2 High-level command name

Suggested internal command name:

```text
tav_uchime3_denovo
```

## 8.3 Side-by-side: original vs extension

| Aspect | Original single-sequence logic | TAV-native extension |
|---|---|---|
| Query object | one denoised sequence | one denoised TAV `(Q_L, Q_R)` |
| Parent object | one denoised sequence | one parent TAV `(A_L, A_R)` or `(B_L, B_R)` |
| Parent model | one-parent or two-parent segmented explanation of one sequence | one-parent or two-parent segmented explanation across two observed anchors plus one unobserved middle region |
| Breakpoint classes | segmented within one full sequence | breakpoint in left anchor, unobserved middle, or right anchor |
| Abundance rule | parents more abundant than chimera | same, but applied to parent TAVs |
| Output | chimeras/nonchimeras FASTA | left FASTA + right FASTA + catalog + chimera report |

## 8.4 Query and parent notation

Query TAV:

```text
Q = (Q_L, Q_R)
```

Candidate parent TAVs:

```text
A = (A_L, A_R)
B = (B_L, B_R)
```

Abundance terms:

```text
abund(Q), abund(A), abund(B)
```

Abundance skew threshold:

```text
gamma
```

Parent eligibility requires:

```text
abund(A) >= gamma * abund(Q)
abund(B) >= gamma * abund(Q)
```

## 8.5 Single-parent null model

Define the best one-parent explanation score as:

```text
H_1(P) = s(Q_L, P_L) + s(Q_R, P_R)
```

where `P` is a candidate parent TAV.

The best single-parent explanation is:

```text
Best_1(Q) = max over candidate parents P of H_1(P)
```

Meaning:

- compare the whole left anchor of the query to the whole left anchor of one parent
- compare the whole right anchor of the query to the whole right anchor of that same parent
- sum the two scores

## 8.6 Three paired-end-native chimera cases

### 8.6.1 Case M: breakpoint in the unobserved middle

Here the left observed anchor is inherited from parent A and the right observed anchor is inherited from parent B.

Score:

```text
H_mid(A, B) = s(Q_L, A_L) + s(Q_R, B_R)
```

Meaning:

- the entire left anchor of `Q` is explained by `A`
- the entire right anchor of `Q` is explained by `B`
- the crossover happens somewhere in the unobserved interior between anchors

### 8.6.2 Case L: breakpoint inside the left anchor

Let `i` be a split position inside `Q_L`.

Score:

```text
H_L(A, B, i) = s_pref(Q_L, A_L, i) + s_suf(Q_L, B_L, i) + s(Q_R, B_R)
```

Meaning:

- the prefix of the query left anchor up to position `i` is explained by parent A
- the suffix of the query left anchor after position `i` is explained by parent B
- the entire right anchor is explained by parent B

Best score for case L:

```text
H_L*(A, B) = max over all valid split positions i of H_L(A, B, i)
```

### 8.6.3 Case R: breakpoint inside the right anchor

Let `j` be a split position inside `Q_R`.

Score:

```text
H_R(A, B, j) = s(Q_L, A_L) + s_pref(Q_R, A_R, j) + s_suf(Q_R, B_R, j)
```

Meaning:

- the entire left anchor is explained by parent A
- the prefix of the right anchor is explained by parent A
- the suffix of the right anchor is explained by parent B

Best score for case R:

```text
H_R*(A, B) = max over all valid split positions j of H_R(A, B, j)
```

## 8.7 Best two-parent explanation

For one ordered parent pair `(A, B)`:

```text
H_2(A, B) = max( H_mid(A, B), H_L*(A, B), H_R*(A, B) )
```

The best two-parent explanation for `Q` is:

```text
Best_2(Q) = max over all eligible ordered parent pairs (A, B), A != B, of H_2(A, B)
```

## 8.8 Chimera evidence measure

Define:

```text
Delta(Q) = Best_2(Q) - Best_1(Q)
```

Interpretation:

- small `Delta(Q)` means one parent explains the query almost as well as any two-parent model
- large `Delta(Q)` means the query is substantially better explained as a chimera

## 8.9 Classification output should record breakpoint class

Because paired-end data yield three biologically distinct evidence classes, classification should include the winning class:

- `LEFT_BREAK`
- `MIDDLE_BREAK`
- `RIGHT_BREAK`

This is not a cosmetic detail; it captures the strength and type of evidence.

## 8.10 Output artifacts

### 8.10.1 Mandatory

- `nonchim_left.fa`
- `nonchim_right.fa`
- `nonchim_catalog.tsv`

### 8.10.2 Optional

- `chim_left.fa`
- `chim_right.fa`
- `chim_catalog.tsv`
- `tav_chimera.tsv`

Suggested `tav_chimera.tsv` fields:

```text
query_tav    parent_A    parent_B    breakpoint_class    best_one_score    best_two_score    delta
```

## 8.11 Important new files for this command

- `tav_chimera3.cc`
- `tav_distance.cc`
- `tav_align.cc`
- `tav_catalog.cc`

## 8.12 Minimal reuse from original logic

Preserve conceptually:

- abundance-sorted de novo processing
- parent must be more abundant than query
- best one-parent vs best two-parent competition

Replace structurally:

- one contiguous sequence segmentation logic
- single-file alignment visualization assumptions
- full-sequence-only parent representation

---

## 9. Command-specific plan: usearch_global -> TAV-usearch_global

## 9.1 Command goal

Assign paired query reads to paired representative TAVs and generate count matrices / BIOM output.

## 9.2 High-level command name

Suggested internal command name:

```text
tav_usearch_global
```

## 9.3 Side-by-side: original vs extension

| Aspect | Original single-sequence logic | TAV-native extension |
|---|---|---|
| Query unit | one FASTA/FASTQ sequence | one normalized read pair -> one `(Q_L, Q_R)` query pair |
| Database unit | one target sequence or UDB target | one TAV representative `(C_L, C_R)` |
| Alignment/scoring | one global sequence-vs-sequence comparison | one left-anchor comparison + one right-anchor comparison |
| Acceptance threshold | one-sequence identity / distance threshold | aggregated pair threshold plus optional per-end thresholds |
| Table row key | target sequence / OTU ID | TAV ID |
| Optional pairwise output | one alignment / hit record | paired hit record, optionally split left/right alignments |

## 9.4 Query construction

Each synchronized paired query is converted into:

```text
Q = (Q_L, Q_R)
```

using the same normalization pipeline used during representative generation.

Representative catalog entry:

```text
C = (C_L, C_R)
```

## 9.5 Pair-level scoring and distance

Define the aggregated pair score as:

```text
S_pair(Q, C) = s(Q_L, C_L) + s(Q_R, C_R)
```

and aggregated pair distance as:

```text
D_pair(Q, C) = d(Q_L, C_L) + d(Q_R, C_R)
```

### 9.5.1 Aggregated identity

If identity-based acceptance is preferred, define:

```text
ID_pair(Q, C) = (matches_left + matches_right) / (aligned_left + aligned_right)
```

where:

- `matches_left` = number of matches in the left-anchor alignment
- `aligned_left` = aligned columns in the left-anchor alignment
- `matches_right` = number of matches in the right-anchor alignment
- `aligned_right` = aligned columns in the right-anchor alignment

### 9.5.2 Recommended acceptance rule

Use both a total and per-end threshold:

```text
ID_pair(Q, C) >= id_total_min
ID_left(Q, C) >= id_end_min
ID_right(Q, C) >= id_end_min
```

This avoids pathological assignments where one end is strong and the other is weak.

## 9.6 Best-hit assignment

For each query pair `Q`:

1. evaluate all candidate TAVs `C`
2. keep those passing the chosen total/per-end thresholds
3. select the best hit by:
   1. highest `S_pair(Q, C)` or highest `ID_pair(Q, C)`
   2. then smallest `D_pair(Q, C)`
   3. then highest representative abundance
   4. then stable TAV ID order

## 9.7 Count matrix semantics

If the best accepted hit is `C*`, increment:

```text
count[sample(Q), tav_id(C*)] += abundance(Q)
```

where `sample(Q)` is parsed from the query header.

This preserves the conceptual role of `otutabout`, but with the row key now representing a paired-anchor TAV rather than one merged sequence.

## 9.8 Output artifacts

### 9.8.1 Mandatory

- `tav_table.tsv`
- optionally `tav_table.biom`

### 9.8.2 Optional

- `query_to_tav.tsv`
- `query_to_tav_left.aln`
- `query_to_tav_right.aln`

Suggested `query_to_tav.tsv` fields:

```text
query_id    tav_id    id_left    id_right    id_total    d_left    d_right    d_total
```

## 9.9 Important new files for this command

- `tav_global.cc`
- `tav_distance.cc`
- `tav_table.cc`
- `tav_results.cc`
- `tav_pairio.cc`
- `tav_catalog.cc`

## 9.10 Minimal reuse from original logic

Preserve conceptually:

- best-hit assignment
- optional OTU/BIOM style count matrices
- optional alignment and match reports

Replace structurally:

- one query record -> one target record assumption
- table row semantics tied to a single target sequence
- one-alignment-per-hit display assumption

---

## 10. Minimal modifications needed vs new code

To avoid overcomplicating the plan, the right decomposition is:

## 10.1 What should be completely new

These should be implemented as new code, not by patching old command logic line-by-line.

- paired observation object model
- paired representative catalog
- pair-level distance aggregation
- TAV-native denoising loop
- TAV-native chimera case logic
- TAV-native global assignment loop
- paired sequence output writers

## 10.2 What can be reused or wrapped

These can be reused as standalone utilities or copied/adapted as needed.

- per-anchor global alignment
- reverse-complement helpers
- FASTA/FASTQ parsing
- sample extraction conventions
- BIOM / OTU table formatting ideas

## 10.3 What the “minimal modifications” principle means in practice

The command-specific plans above are designed so that each new command is built by wrapping the same single-anchor primitive twice, then adding a pair-level decision layer.

That is the minimal-change mental model:

- **single-end primitive stays single-end**
- **pair-level command logic becomes new**

This avoids rewriting every low-level detail while still making paired-end native to the command semantics.

---

## 11. Suggested implementation roadmap

## 11.1 Phase 0: shared infrastructure

Implement first:

- `tav_types.*`
- `tav_pairio.*`
- `tav_catalog.*`
- `tav_distance.*`
- optional `tav_align.*`

Deliverables:

- read synchronized pairs
- normalize orientation
- cut anchors
- write paired FASTA + catalog
- scalar distance/score between anchor pairs

## 11.2 Phase 1: exact paired dereplication

Implement exact dereplication of `(left_anchor, right_anchor)` pairs.

Deliverables:

- abundance-sorted unique pair list
- representative headers
- test fixtures showing that exact-pair dereplication preserves pair structure

## 11.3 Phase 2: TAV-usearch_global

Build the easiest full end-to-end command first.

Why first:

- easiest to validate manually
- uses the same paired query normalization needed everywhere
- produces immediately useful count tables
- establishes the canonical TAV catalog format

## 11.4 Phase 3: TAV-UNOISE3

Add denoising over exact paired dereplication.

Validation targets:

- pair-level distance behaves as intended
- abundance-vs-distance rule is monotone and interpretable
- outputs stable paired representatives

## 11.5 Phase 4: TAV-UCHIME3-denovo

Add paired-end-native chimera detection.

Validation targets:

- synthetic chimeras with breakpoint in left anchor
- synthetic chimeras with breakpoint in right anchor
- synthetic chimeras with breakpoint in unobserved middle
- nonchimeric controls where the best one-parent model remains competitive

## 11.6 Phase 5: polish and unification

Add:

- cleaner command entrypoints
- simplified CLI
- optional alignment outputs
- BIOM export
- regression tests

---

## 12. Testing plan

## 12.1 Unit tests

- reverse-complement + orientation normalization
- anchor extraction
- exact pair dereplication
- pair distance aggregation
- prefix/suffix split scoring on one anchor
- catalog round-trip tests

## 12.2 Property tests

- identical pairs have zero pair distance
- pair distance is symmetric
- middle-break chimera score is invariant to hidden middle length because it uses only anchors
- left-break and right-break models reduce to whole-anchor scores when split is at extremes

## 12.3 Synthetic biological tests

- nonchimeric parent pair controls
- left-anchor crossover controls
- right-anchor crossover controls
- unobserved-middle crossover controls
- unrelated-end pairing controls

## 12.4 End-to-end tests

- paired FASTQ -> TAV-UNOISE3 -> TAV-UCHIME3-denovo -> TAV-usearch_global -> table

---

## 13. Deliverables summary

## 13.1 Core commands

- `tav_unoise3`
- `tav_uchime3_denovo`
- `tav_usearch_global`

## 13.2 Core outputs

For representatives:

- left-anchor FASTA
- right-anchor FASTA
- TAV catalog TSV

For denoising/chimera/search diagnostics:

- query-to-TAV mapping TSV
- chimera explanation TSV
- optional split left/right alignments

For abundance output:

- count matrix TSV
- optional BIOM

---

## 14. Final recommendation

The cleanest way to keep the project faithful to the original methods while making paired-end data first-class is:

1. define TAV as an ordered pair of oriented anchors
2. build all three commands around that object directly
3. reuse single-anchor scoring/alignment primitives
4. create new pair-native command files rather than patching existing single-sequence command files
5. use a catalog as the canonical representative artifact

This gives a coherent paired-end-native system without pretending that non-overlapping read pairs are equivalent to merged reads.

---

## 15. Source context used while drafting this plan

This plan was drafted against the current public VSEARCH command/model context:

- `--cluster_unoise` performs denoising but leaves de novo chimera removal to `--uchime3_denovo`
- `--uchime3_denovo` is abundance-aware and uses higher default `abskew` than older denovo modes
- `--usearch_global` supports FASTA/FASTQ query input and can emit `otutabout` / BIOM-style outputs
- `otutabout` represents a sample-by-target abundance matrix

The concrete command and behavior assumptions in this document should therefore be read as a paired-end-native extension of those public semantics, not as a proposal to preserve every internal implementation detail.
