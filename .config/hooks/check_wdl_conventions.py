#!/usr/bin/env python3
"""PHB WDL conventions that need more than a line-oriented regex.

  1. every `runtime` block sets BOTH `disks` (GCP/Terra) and `disk` (TES/Cromwell)
  2. every docker image string is pinned to an explicit tag

Rules that ARE expressible as a single-line regex live as pygrep hooks in
.config/.pre-commit-config.yaml instead — keep them there.
"""
import re
import sys

# Non-greedy body up to the first line that is only a closing brace at the
# indentation the block opened at — WDL runtime blocks are never nested.
RUNTIME_BLOCK = re.compile(r"runtime\s*\{(.*?)\n\s*\}", re.DOTALL)
# Only declarations pin an image; `docker: "~{docker}"` inside a runtime block
# just forwards the input, so skip anything containing a placeholder.
DOCKER_STRING = re.compile(r'docker\w*\s*[:=]\s*"([^"]+)"')
PLACEHOLDER = re.compile(r"[~$]\{")

def check(path: str) -> list[str]:
  errors = []
  text = open(path, encoding="utf-8").read()

  for match in RUNTIME_BLOCK.finditer(text):
    body = match.group(1)
    line = text[: match.start()].count("\n") + 1
    for key, why in (("disks", "GCP/Terra"), ("disk", "TES/Cromwell")):
      if not re.search(rf"^\s*{key}\s*:", body, re.MULTILINE):
        errors.append(f"{path}:{line}: runtime block is missing `{key}:` ({why})")

  for match in DOCKER_STRING.finditer(text):
    image = match.group(1)
    if PLACEHOLDER.search(image):
      continue
    line = text[: match.start()].count("\n") + 1
    # A tag lives after the last "/", so "us-docker.pkg.dev/foo/bar:1.0" is
    # pinned but "ubuntu" and "quay.io/org/tool" are not.
    if ":" not in image.rsplit("/", 1)[-1] and "@sha256" not in image:
      errors.append(f"{path}:{line}: docker image `{image}` has no explicit tag")

  return errors

def main(paths: list[str]) -> int:
  errors = [e for p in paths for e in check(p)]
  for e in errors:
    print(e, file=sys.stderr)
  return 1 if errors else 0

if __name__ == "__main__":
  raise SystemExit(main(sys.argv[1:]))
