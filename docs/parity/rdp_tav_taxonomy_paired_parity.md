# RDP Taxonomy Paired-End Parity (TAV)

## Scope

This document compares:

- stock RDP classifier behavior (single-sequence taxonomy inference)
- `vsearch_plus` paired-end TAV taxonomy extension

Goal: keep the stock classifier backbone unchanged while adding first-class paired-end TAV handling.

## Command Surface

### Stock baseline

Stock invocation classifies one sequence stream per run:

```bash
java -jar classifier.jar classify -f allrank -o out.tsv input.fa
```

### Paired TAV extension

`vsearch_plus` adds a paired command wrapper:

```bash
python3 rdp_tav_taxonomy.py \
  --left-fasta tav_left.fa \
  --right-fasta tav_right.fa \
  --output-prefix tav_taxonomy \
  --pair-filter both \
  --min-conf 0.8
```

or catalog mode:

```bash
python3 rdp_tav_taxonomy.py \
  --tav-catalog tav_denoised.tsv \
  --output-prefix tav_taxonomy_catalog
```

## Step-by-Step Backbone Parity

### 1. Classifier algorithm parity: preserved

- Each end is classified by the stock RDP Java jar.
- No algorithm changes are made to `Classifier.classify(...)`.
- Orientation logic, word extraction, bootstrap scoring, and confidence computation remain stock behavior.

### 2. Model parity: preserved and explicit

- Extension uses pretrained model properties from downloaded RDP assets:
  - `data/.../classifier/16srrna/rRNAClassifier.properties`
- This is wired through stock `-t` training-property input.

### 3. Per-end output parity: preserved

- Raw per-end outputs are still stock `allrank` tables:
  - `PREFIX.left.allrank.tsv`
  - `PREFIX.right.allrank.tsv`
- Format is unchanged from stock output conventions.

### 4. Paired aggregation: extension-only layer

- New logic combines per-end rank assignments into one TAV-level lineage.
- `--pair-filter both`:
  - require both ends confident and agreeing at each rank
- `--pair-filter any`:
  - allow one confident end; if both confident and conflict, choose higher-confidence end

### 5. Prefix truncation rule: explicit extension behavior

- Once one rank cannot be assigned in paired output, downstream ranks are blank.
- This avoids creating lineage paths that skip unresolved intermediate ranks.

### 6. Count table extension

- New `PREFIX.rank_counts.tsv` summarizes paired assignments by rank and abundance.

## Real Data Example (from `data/real_out2`)

### Stock-style per-end records (unchanged)

From `tav_taxonomy.left.allrank.tsv`:

- `TAV000001;size=8242` -> genus `Pseudarthrobacter` (0.98)
- `TAV000002;size=310` -> genus `Pseudomonas` (0.93)

### Paired output behavior (`--pair-filter both`, `--min-conf 0.8`)

From `tav_taxonomy.paired.tsv`:

- `TAV000001`
  - left genus: `Pseudarthrobacter` (0.98)
  - right genus: `Pseudarthrobacter` (0.55)
  - paired genus: blank (right confidence below 0.8)
  - paired lineage stops at family `Micrococcaceae`

- `TAV000002`
  - left genus: `Pseudomonas` (0.93)
  - right genus: `Pseudomonas` (0.93)
  - paired genus: `Pseudomonas` (0.93, source `both`)

- `TAV000003`
  - left genus: `Collimonas` (0.70)
  - right genus: `Noviherbaspirillum` (0.07)
  - paired lineage stops at order due to insufficient/conflicting support

## Intentional Disparities

These are intentional extension additions, not parity gaps:

1. Paired TAV input handling (`left+right` or catalog)
2. Pair aggregation policy (`--pair-filter any|both`)
3. Paired consolidated table and rank count outputs

## Conclusion

Parity is achieved for the inference backbone: stock RDP Java classification is reused directly. New behavior is isolated to pair-aware aggregation and reporting.
