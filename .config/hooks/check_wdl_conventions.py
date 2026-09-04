#!/usr/bin/env python3
"""PHB WDL conventions that need more than a line-oriented regex.

  1. every `runtime` block sets BOTH `disks` (GCP/Terra) and `disk` (TES/Cromwell)
  2. every docker image string is pinned to an explicit tag, in Theiagen's GAR
  3. WDL structure is indented 2 spaces per nesting level

(3) is formatting: it auto-fixes, and is always correct to apply. (1) and (2)
are release policy — half-written work is legitimately allowed to fail them —
so they only report, and run at the manual stage.

  --fix      formatting only (3); the commit-stage hook
  --report   policy only (1, 2); the CI hook
  neither    both

Exits non-zero when it fixes something, so the commit stops and you re-stage.
Single-line regex rules live as pygrep hooks in .pre-commit-config.yaml.
"""
import re
import sys

RUNTIME_BLOCK = re.compile(r"runtime\s*\{(.*?)\n\s*\}", re.DOTALL)
DOCKER_STRING = re.compile(r'docker\w*\s*[:=]\s*"([^"]+)"')
# Mirror third-party images here rather than referencing quay.io or Docker Hub,
# which can delete or retag out from under a workflow.
GAR_PREFIX = "us-docker.pkg.dev/general-theiagen/"
# `docker: "~{docker}"` forwards the input rather than pinning an image.
PLACEHOLDER = re.compile(r"[~$]\{")
# Capture the opener so we know what closes it; treating a `command {` block as
# heredoc style would exempt the rest of the file from check 3.
COMMAND_OPEN = re.compile(r"command\s*(<<<|\{)")
COMMAND_CLOSE = {"<<<": ">>>", "{": "}"}
STRING_LITERAL = re.compile(r'"(?:[^"\\]|\\.)*"')

def check(path: str, fix: bool, report: bool) -> list[str]:
  errors = []
  text = open(path, encoding="utf-8").read()

  # Fix first; indentation never adds or removes lines, so the line numbers
  # reported below stay valid either way.
  if fix:
    fixed_text, fixes = fix_indentation(path, text)
    if fixes:
      with open(path, "w", encoding="utf-8") as handle:
        handle.write(fixed_text)
      text = fixed_text
      errors += fixes

  if not report:
    return errors

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
    # A tag lives after the last "/": "gcr.io/foo/bar:1.0" is pinned, "ubuntu" is not.
    if ":" not in image.rsplit("/", 1)[-1] and "@sha256" not in image:
      errors.append(f"{path}:{line}: docker image `{image}` has no explicit tag")
    if not image.startswith(GAR_PREFIX):
      errors.append(
        f"{path}:{line}: docker image `{image}` is not in GAR; "
        f"mirror it under {GAR_PREFIX} first")

  return errors

def _structure(line: str) -> str:
  """The line with string literals and trailing comments blanked out, so only
  syntactic brackets are counted. A `String x = "}"` would otherwise skew every
  indent below it."""
  return STRING_LITERAL.sub('""', line).split("#", 1)[0]

def fix_indentation(path: str, text: str) -> tuple[str, list[str]]:
  """Indent 2 spaces per nesting level; returns the corrected text and one
  message per line rewritten.

  Depth is `{}` nesting plus one level per enclosing `input:` — that label
  opens a block visually but not syntactically, so braces alone under-predict
  every call argument by one level.

  Command-block bodies are exempt (their shell and Python have their own
  conventions). Lines continuing an unbalanced `(` or `[` are hand-aligned, so
  they are only held to an even indent.
  """
  errors = []
  lines = text.split("\n")
  brace = 0          # {} nesting depth
  expr = 0           # ()/[] depth; > 0 means the previous line ran on
  input_levels = []  # brace depths at which an `input:` added a level
  closer = None

  for lineno, line in enumerate(lines, 1):
    stripped = line.strip()
    if closer:
      if stripped.startswith(closer):
        closer = None
      continue
    opener = COMMAND_OPEN.search(line)
    if opener:
      closer = COMMAND_CLOSE[opener.group(1)]
      continue
    if not stripped:
      continue

    body = _structure(line)
    indent = len(line) - len(line.lstrip(" "))

    if expr:
      expected = indent - (indent % 2)
    else:
      # A line opening with a closer belongs one level out from the body.
      level = brace - (1 if stripped[0] in "})]" else 0)
      while input_levels and level < input_levels[-1]:
        input_levels.pop()
      expected = 2 * (level + len(input_levels))

    if indent != expected:
      # Leading indent only; trailing whitespace is another hook's job.
      lines[lineno - 1] = " " * expected + line.lstrip(" ")
      errors.append(f"{path}:{lineno}: fixed indent {indent} -> {expected}")

    # Depths advance for every line, including skipped ones — a skipped line
    # that opens a block must still push, or everything below it drifts.
    brace = max(brace + body.count("{") - body.count("}"), 0)
    expr = max(expr + body.count("(") - body.count(")")
               + body.count("[") - body.count("]"), 0)
    if stripped == "input:":
      input_levels.append(brace)

  return "\n".join(lines), errors

def main(argv: list[str]) -> int:
  fix = "--report" not in argv
  report = "--fix" not in argv
  paths = [a for a in argv if not a.startswith("--")]
  errors = [e for p in paths for e in check(p, fix, report)]
  for e in errors:
    print(e, file=sys.stderr)
  return 1 if errors else 0

if __name__ == "__main__":
  raise SystemExit(main(sys.argv[1:]))
