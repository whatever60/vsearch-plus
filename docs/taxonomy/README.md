# TAV Taxonomy Assignment (RDP)

This directory documents and supports taxonomy assignment for paired-end TAV sequences.

## Quickstart

1. Download or refresh RDP assets:

```bash
python3 get_rdp_classifier.py --output-root data/third_party/rdp_classifier
```

2. Run paired TAV taxonomy on real TAV outputs:

```bash
python3 rdp_tav_taxonomy.py \
  --left-fasta data/real_out2/tav_left.fa \
  --right-fasta data/real_out2/tav_right.fa \
  --output-prefix data/real_out2/tav_taxonomy \
  --pair-filter both \
  --min-conf 0.8 \
  --java-opt=-Xmx2g
```

3. Optional catalog-driven mode:

```bash
python3 rdp_tav_taxonomy.py \
  --tav-catalog data/real_out2/tav_denoised.tsv \
  --output-prefix data/real_out2/tav_taxonomy_catalog \
  --pair-filter any \
  --min-conf 0.8 \
  --java-opt=-Xmx2g
```

## Outputs

Given `--output-prefix PREFIX`, the command writes:

- `PREFIX.left.allrank.tsv`: raw RDP output for left anchors
- `PREFIX.right.allrank.tsv`: raw RDP output for right anchors
- `PREFIX.paired.tsv`: paired TAV-level taxonomy table after aggregation
- `PREFIX.rank_counts.tsv`: abundance-weighted taxon counts by rank
