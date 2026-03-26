# TAV Taxonomy Assignment (RDP)

This module performs taxonomy assignment for TAV sequences as native paired-end objects.

## Quickstart

1. Download or refresh RDP assets:

```bash
python3 get_rdp_classifier.py --output-root data/third_party/rdp_classifier
```

2. Run native paired classification from split paired inputs:

```bash
python3 rdp_tav_taxonomy.py \
  --input data/real_out2/tav_left.fa \
  --input2 data/real_out2/tav_right.fa \
  --output data/real_out2/tav_taxonomy_native.tsv \
  --format allrank \
  --conf 0.8 \
  --min-words 5 \
  --java-opt=-Xmx2g
```

3. Optional interleaved mode:

```bash
python3 rdp_tav_taxonomy.py \
  --input interleaved_tav.fa \
  --interleaved \
  --output tav_taxonomy_native.tsv
```

## Output

The command writes one stock-compatible classification stream where each line is one TAV record.

- `--output`: primary taxonomy output (`allrank`, `fixrank`, `filterbyconf`, or `db`)
- `--shortseq-outfile`: optional IDs of short/unclassifiable paired records

No per-anchor merge outputs are produced.
