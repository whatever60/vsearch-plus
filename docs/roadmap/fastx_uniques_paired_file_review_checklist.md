# Native Paired `fastx_uniques`: File-by-File Review Checklist

Reviewed files/modules:

- `src/derep.cc`
- `src/derep_paired.cc`
- `src/derep.h`
- `src/derep_paired.h`
- `src/vsearch.cc`
- `src/tav_extension.cc` (legacy paired path)

Status legend:

- `Good`: stock-owned or paired-owned in the right place, with the right responsibility.
- `Needs closer stock parity`: the paired behavior is correct in spirit, but still diverges more than necessary from stock structure.
- `Should be deleted/refolded`: legacy helper or path that should no longer own paired `fastx_uniques`.

## File-Level Notes

- Stock `fastx_uniques` is owned by `src/derep.cc`.
- Native paired `fastx_uniques` is now owned by `src/derep_paired.cc`.
- CLI routing in `src/vsearch.cc` now chooses between stock `derep(...)` and native `fastx_uniques_paired(...)` based on whether paired input is present.
- `src/tav_extension.cc` no longer owns the paired command path and is no longer part of the build.

## `src/derep.cc`

| Function | Role | Status | Review note |
| :-- | :-- | :-- | :-- |
| `derep(...)` | Stock owner for single-end `fastx_uniques` | Good | Remains the stock single-end implementation and the structural reference for the paired port. |

## `src/derep_paired.h`

| Declaration | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `fastx_uniques_paired(...)` | paired public analog of stock `derep(...)` ownership | Good | Correct public home for the native paired command surface. |

## `src/derep_paired.cc`

| Function / object | Stock counterpart or role | Status | Review note |
| :-- | :-- | :-- | :-- |
| `fastx_uniques_paired(...)` | derep-owned paired command implementation | Good | Correct native owner for paired `fastx_uniques`; handles split/interleaved paired input, paired key aggregation, and paired FASTA/catalog outputs. |
| `get_anchor_len_paired` lambda | none in stock single-end path | Good | Necessary paired-only helper because paired anchors depend on both ends and `--fastq_trunclen`. |
| `append_unique_record` lambda | inner derep update path | Good | Necessary local paired aggregation logic; keeps the implementation compact instead of inventing a separate legacy helper layer. |

## `src/derep.h`

| Declaration | Role | Status | Review note |
| :-- | :-- | :-- | :-- |
| `derep(...)` | stock single-end public surface | Good | Kept unchanged, which is what we want. |

## `src/vsearch.cc`

| Branch | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `--fastx_uniques` single-end route | stock `derep(...)` | Good | Unchanged. |
| `--fastx_uniques` paired route | native `fastx_uniques_paired(...)` | Good | Correct native dispatch; the paired command no longer routes to `tav_fastx_uniques(...)`. |

## `src/tav_extension.cc` (legacy paired `fastx_uniques` section)

| Function / object | Target home | Status | Review note |
| :-- | :-- | :-- | :-- |
| `tav_fastx_uniques(...)` | `src/derep_paired.cc:fastx_uniques_paired(...)` | Should be deleted/refolded | The native paired command no longer routes here and the file is no longer part of the build. |

## Checklist Goal

This checklist is complete when all of these are true:

- paired `fastx_uniques` no longer routes to `tav_fastx_uniques(...)`
- native paired ownership sits under `src/derep_paired.cc`
- build/link no longer includes `tav_extension.o`
- split and interleaved paired smoke outputs agree
