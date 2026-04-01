# Native Paired `fastx_uniques`: File-by-File Review Checklist

Status legend:

- `Good`: already close to stock shape, with only explicit single-read to paired-read differences.
- `Needs closer stock parity`: there is a clear stock counterpart, but the paired implementation still diverges more than necessary.
- `Should be deleted/refolded`: should be inlined, merged into a stock-shaped function, or moved to a more appropriate shared module.

## `cpp/src/derep.cc` / `cpp/src/derep.h`

| Function / declaration | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `derep(...)` | `derep(...)` | Good | Stock owner for single-end `fastx_uniques`; this is the explicit function the paired port mirrors. |

## `cpp/src/derep_paired.cc` / `cpp/src/derep_paired.h`

| Function / declaration | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `derep_paired(...)` | `derep(...)` | Good | Correct public and implementation-level paired analog of stock `derep(...)`. The temporary helper lambdas were refolded so this module now presents a single stock-analog entrypoint. |

## `cpp/src/vsearch.cc`

| Function / branch | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `--fastx_uniques` single-end route | `derep(...)` | Good | Unchanged stock dispatch. |
| `--fastx_uniques` paired route | `derep(...)` | Good | Correct native dispatch to the stock-shaped paired analog `derep_paired(...)`. |
