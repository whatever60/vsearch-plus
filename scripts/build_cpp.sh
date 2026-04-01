#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="$repo_root/cpp"
build_dir="$repo_root/build/cpp"

mkdir -p "$build_dir"

find "$source_dir/src" -name '*.o' -delete
find "$source_dir/src" -name '*.a' -delete
find "$source_dir/src" -name '.deps' -type d -prune -exec rm -rf {} +
rm -f "$source_dir/Makefile" "$source_dir/src/Makefile" "$source_dir/man/Makefile"

if [[ ! -x "$source_dir/configure" ]]; then
  (
    cd "$source_dir"
    ./autogen.sh
  )
fi

(
  cd "$build_dir"
  CFLAGS="${CFLAGS:--O2}" CXXFLAGS="${CXXFLAGS:--O2}" "$source_dir/configure" "$@"
  make ARFLAGS="${ARFLAGS:-cr}"
)
