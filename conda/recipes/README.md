# Conda Recipes

This directory contains the local recipe scaffolding for the Bioconda packaging work around:

- the compiled `vsearch-plus` executable
- the public `rdp-classifier` launcher
- the internal RDP downloader helper used by `rdp-classifier`

Scaffolding directories are already present for:

- `conda/recipes/vsearch-plus-cpp/`
- `conda/recipes/vsearch-plus-java/`
- `conda/recipes/vsearch-plus-python/`
- `conda/recipes/vsearch-plus/`

The aggregate recipe now lives under `conda/recipes/vsearch-plus/`.
The split package directories remain as placeholders in case the packaging strategy changes later.
