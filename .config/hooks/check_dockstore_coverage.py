#!/usr/bin/env python3
"""Every release workflow must be registered in .dockstore.yml.

`workflows/utilities/` is excluded: those files are subworkflows imported by
other workflows, not independently runnable entrypoints, so they are
deliberately unregistered. If that ever stops being true, edit EXCLUDE_PREFIXES.
"""
import re
import sys
from pathlib import Path

EXCLUDE_PREFIXES = ("workflows/utilities/",)

def main() -> int:
  dockstore = Path(".dockstore.yml")
  if not dockstore.is_file():
    print("error: .dockstore.yml not found", file=sys.stderr)
    return 1

  # Paths in .dockstore.yml are absolute-from-repo-root ("/workflows/..."),
  # so strip the leading slash to compare against on-disk relative paths.
  registered = {
    p.lstrip("/")
    for p in re.findall(r"primaryDescriptorPath:\s*(\S+)", dockstore.read_text())
  }

  missing = sorted(
    str(p)
    for p in Path("workflows").rglob("wf_*.wdl")
    if not str(p).startswith(EXCLUDE_PREFIXES) and str(p) not in registered
  )

  if missing:
    print("Workflows missing a .dockstore.yml entry:", file=sys.stderr)
    for path in missing:
      print(f"  {path}", file=sys.stderr)
    print(
      "\nAdd a `- subclass: WDL` block with primaryDescriptorPath: /<path> and\n"
      "testParameterFiles: /tests/inputs/empty.json, or add the path to\n"
      "EXCLUDE_PREFIXES in this script if it is an internal subworkflow.",
      file=sys.stderr,
    )
    return 1
  return 0

if __name__ == "__main__":
  raise SystemExit(main())
