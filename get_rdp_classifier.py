#!/usr/bin/env python3
"""Compatibility wrapper for python/get_rdp_classifier.py."""

import pathlib
import sys

repo_root = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(repo_root / "python"))

from vsearch_plus.rdp_download import main


main()
