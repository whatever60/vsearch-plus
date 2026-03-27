#!/usr/bin/env python3

"""Build compiled docs in a logical order and convert to DOCX."""

from datetime import datetime, timezone
from pathlib import Path
import re
import subprocess


repo_root = Path(__file__).resolve().parents[2]
docs_dir = repo_root / "docs"
readme_path = docs_dir / "README.md"
compiled_md_path = docs_dir / "compiled_docs.md"
compiled_docx_path = docs_dir / "compiled_docs.docx"

readme_text = readme_path.read_text(encoding="utf-8")

ordered_rel_paths = re.findall(r"- `([^`]+\.md)`", readme_text)
ordered_rel_paths = [path for path in ordered_rel_paths if path.startswith("parity/")]

seen = set()
ordered_unique_rel_paths = []
for rel_path in ordered_rel_paths:
  if rel_path not in seen:
    ordered_unique_rel_paths.append(rel_path)
    seen.add(rel_path)

sections = [
  "# VSEARCH Plus Documentation",
  "",
  f"_Generated on {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} from Markdown files under `docs/parity/`._",
  "",
]

for i, rel_path in enumerate(ordered_unique_rel_paths):
  source_path = docs_dir / rel_path
  source_text = source_path.read_text(encoding="utf-8").rstrip()
  sections.append(f"## Source: `docs/{rel_path}`")
  sections.append("")
  sections.append(source_text)
  sections.append("")
  if i < len(ordered_unique_rel_paths) - 1:
    sections.append("\\newpage")
    sections.append("")

compiled_markdown = "\n".join(sections).rstrip() + "\n"
compiled_md_path.write_text(compiled_markdown, encoding="utf-8")

subprocess.run(
  [
    "pandoc",
    str(compiled_md_path),
    "--from",
    "gfm",
    "--to",
    "docx",
    "--output",
    str(compiled_docx_path),
  ],
  check=True,
)

print(f"Wrote {compiled_md_path}")
print(f"Wrote {compiled_docx_path}")
