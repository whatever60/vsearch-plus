#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

conda env create -f "$repo_root/environment.yml" || conda env update -f "$repo_root/environment.yml"
echo "Activate with: conda activate vsearch-plus"
