#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ $# -ge 1 && -n "$1" ]]; then
  prefix="$1"
elif [[ -n "${PREFIX:-}" ]]; then
  prefix="$PREFIX"
else
  echo "Usage: ./scripts/stage_install.sh <prefix>" >&2
  echo "Or set PREFIX in the environment." >&2
  exit 1
fi

prefix="$(cd "$(dirname "$prefix")" && pwd)/$(basename "$prefix")"

"$repo_root/scripts/build_cpp.sh"

install -d "$prefix/bin"
install -d "$prefix/share/vsearch-plus/scripts"
install -d "$prefix/share/vsearch-plus/java/src/main"

install -m 755 "$repo_root/build/cpp/bin/vsearch" \
  "$prefix/bin/vsearch-plus"
install -m 755 "$repo_root/scripts/rdp-classifier" \
  "$prefix/bin/rdp-classifier"
install -m 755 "$repo_root/scripts/get_rdp_classifier.py" \
  "$prefix/share/vsearch-plus/scripts/get_rdp_classifier.py"

rm -rf "$prefix/share/vsearch-plus/java/src/main/java"
cp -a "$repo_root/java/src/main/java" \
  "$prefix/share/vsearch-plus/java/src/main/java"
