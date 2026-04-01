# TAV Taxonomy Assignment (RDP)

This module performs taxonomy assignment for TAV sequences as native paired-end objects.

## Quickstart

1. Run native paired classification from split paired inputs:

```bash
./scripts/rdp-classifier \
  --input data/real_out2/tav_left.fa \
  --input2 data/real_out2/tav_right.fa \
  --output data/real_out2/tav_taxonomy_native.tsv \
  --format allrank \
  --conf 0.8 \
  --min-words 5 \
  --java-opt=-Xmx2g
```

The first real `rdp-classifier` run downloads the required RDP assets automatically.

2. Optional interleaved mode:

```bash
./scripts/rdp-classifier \
  --input interleaved_tav.fa \
  --interleaved \
  --output tav_taxonomy_native.tsv
```

## Output

The command writes one stock-compatible classification stream where each line is one TAV record.

- Primary output option (stock aliases supported):
  - `--output` / `--outputFile` / `-o`
- Supported stock-style config aliases:
  - `--format` / `-f` (`allrank`, `fixrank`, `filterbyconf`, `db`)
  - `--conf` / `-c`
  - `--min-words` / `--minWords` / `-w`
  - `--shortseq-outfile` / `--shortseq_outfile` / `-s`
  - `--bootstrap_outfile` / `-b` (bootstrap summary output)
  - `--hier_outfile` / `-h` (hierarchy count output)
  - `--train-prop` / `--train_propfile` / `-t`
  - `--gene` / `-g`
  - `--queryFile` / `-q` (accepted, ignored)

Not currently supported in paired mode:

- `--metadata` / `-d`
- `--biomFile` / `-m`
- `--format biom`

No per-anchor merge outputs are produced.
