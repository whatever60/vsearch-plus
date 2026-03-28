# Native Paired `fastx_uniques`: SOP and Implementation Checklist

This document mirrors the paired-native workflow we used for `cluster_unoise`, `usearch_global`, and `uchime3_denovo`, but for `vsearch --fastx_uniques`.

It is based on:

- `docs/parity/fastx_uniques_paired_extension.md`
- `docs/callstacks/fastx_uniques.txt`
- `src/derep.cc`
- `src/derep_paired.cc`

## Current Status

- Paired `--fastx_uniques` now routes to `derep_paired()` in `src/derep_paired.cc`.
- Single-end `--fastx_uniques` still routes to stock `derep(...)` in `src/derep.cc`.
- The new paired public surface lives in `src/derep_paired.h`.
- the now-removed legacy `tav_extension.cc` path is no longer on the paired `fastx_uniques` CLI path and is no longer part of the build.
- `src/derep_paired.cc` now presents a single stock-shaped paired entrypoint, `derep_paired(...)`; the temporary helper lambdas were refolded into that function.
- Native paired smoke tests passed for:
  - split input (`R1 R2`)
  - interleaved input (`--interleaved`)
  - split/interleaved output agreement for `--tabbedout`, `--fastaout`, and `--fastaout_rev`

## Purpose

The goal is to make paired `fastx_uniques` look like a stock-owned derep command extension, not like a lingering custom sidecar in the now-removed legacy `tav_extension.cc` path.

The working rule stays the same:

```text
copy stock first
rename with _paired
replace one read with two reads
reuse stock primitives
justify every divergence
```

## SOP

### 1. Freeze stock ownership and callstack first

For `fastx_uniques`, the stock owner is `src/derep.cc`, not `src/unique.cc`.

That means the paired native port belongs beside `derep`, not beside low-level uniqueness helpers.

### 2. Define paired semantics before porting

Paired `fastx_uniques` turns one derep unit into one synchronized read pair:

- anchor from R1
- anchor from R2
- exact key is `(left_anchor, right_anchor)`
- abundance is aggregated only when both anchors match exactly

Shared paired rules:

- split mode reads R1 and R2 in lockstep
- interleaved mode consumes `(R1,R2)` tuples
- anchor length is `min(len(R1), len(R2))`, also bounded by `--fastq_trunclen`
- R2 stays in native read orientation in external outputs

### 3. Port by stock module boundary

Target module mapping:

- `src/derep.cc` -> `src/derep_paired.cc`
- `src/derep.h` -> `src/derep_paired.h`

Do not move this port into `src/unique.cc`, because `unique.cc` is not the stock command owner.

### 4. Keep paired naming stock-shaped

Main mapping:

- stock `derep(...)`
- paired `derep_paired(...)`

For this command the paired wrapper can stay narrow, because the stock command owner is already a single top-level routine rather than a large thread/search engine.

### 5. Reuse stock primitives aggressively

The native paired port should reuse stock helpers wherever the behavior is still the same:

- `fastx_open`
- `fastx_next`
- `fastx_get_*`
- `fopen_output`
- `fasta_print_general`
- stock parser behavior for `--notrunclabels`
- stock abundance parsing semantics (`--sizein` aware FASTX abundance)

### 6. Allow only truly necessary paired-only helpers

For paired `fastx_uniques`, the necessary differences are narrow:

- paired input synchronization
- shared anchor-length calculation across both ends
- exact paired key formation
- split paired FASTA outputs
- paired catalog schema

Anything beyond that should be questioned.

### 7. Retire the legacy extension path cleanly

Completion criteria for this command:

- paired CLI routing no longer calls `tav_fastx_uniques(...)`
- the build no longer compiles or links the now-removed legacy `tav_extension.cc` path
- the paired command lives in derep-owned files
- split and interleaved paired smoke tests agree

## Concrete Implementation Checklist

1. Document the stock and paired callstacks in `docs/callstacks/fastx_uniques.txt`.
2. Create `src/derep_paired.h` and declare `derep_paired(...)` there.
3. Create `src/derep_paired.cc` beside `src/derep.cc`.
4. Route paired `--fastx_uniques` input in `src/vsearch.cc` to `derep_paired(...)`.
5. Keep single-end routing on stock `derep(...)` unchanged.
6. Wire `src/derep_paired.cc` and `src/derep_paired.h` into the build.
7. Remove `tav_extension.cc` from the build.
8. Verify the final link line no longer includes `tav_extension.o`.
9. Run split paired smoke tests.
10. Run interleaved paired smoke tests.
11. Confirm split/interleaved output agreement.
12. Update parity and roadmap docs so they describe the native derep-owned path rather than the retired tav path.

## Definition Of Done

This command is done when all of these are true:

- paired `fastx_uniques` is owned by `src/derep_paired.cc`
- stock single-end `fastx_uniques` remains on `src/derep.cc`
- paired split and interleaved input both work
- `--tabbedout`, `--fastaout`, and `--fastaout_rev` work in paired mode
- the link line no longer includes `tav_extension.o`
- docs describe the native path accurately
