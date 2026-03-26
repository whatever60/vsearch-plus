#!/usr/bin/env python3
"""Compatibility wrapper for python/rdp_tav_taxonomy.py."""

import pathlib
import sys

repo_root = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(repo_root / "python"))

from vsearch_plus.rdp_tav_taxonomy import main


main()
