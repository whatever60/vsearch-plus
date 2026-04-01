#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
launcher_src="$repo_root/java/src/main/java/org/vsearchplus/rdp/RdpTavTaxonomyMain.java"
launcher_dir="$repo_root/build/java/launcher_classes"
javac_bin="${JAVAC_BIN:-javac}"

mkdir -p "$launcher_dir"
"$javac_bin" -d "$launcher_dir" "$launcher_src"
