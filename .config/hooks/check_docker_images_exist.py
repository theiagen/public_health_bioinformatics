#!/usr/bin/env python3
"""Confirm every docker image in a WDL exists in Theiagen's GAR.

Fetches the image MANIFEST only via the OCI distribution API — nothing is
pulled and no docker daemon is needed. GAR serves these anonymously, so there
is no token handshake.

Exit policy, so this never blocks work the author didn't break:
  404                    -> fail; the tag was typo'd or never pushed
  401/403                -> warn, pass; the repo is not publicly readable, and
                            missing and private look alike
  DNS/timeout/connection -> warn, pass; offline commits must still work

Only GAR images are found, by prefix. Images anywhere else are wdl-conventions'
to report.
"""
import concurrent.futures
import re
import sys
import urllib.error
import urllib.request

REGISTRY = "us-docker.pkg.dev"
# The trailing class is what is legal in an image reference, so a "~{...}"
# placeholder simply does not match and needs no separate guard.
GAR_IMAGE = re.compile(rf'"{re.escape(REGISTRY)}/(general-theiagen/[\w./:@-]+)"')
# Without these the registry returns a v1 manifest, or 404s on a multi-arch
# image that only has an index.
ACCEPT = ", ".join([
  "application/vnd.oci.image.index.v1+json",
  "application/vnd.oci.image.manifest.v1+json",
  "application/vnd.docker.distribution.manifest.list.v2+json",
  "application/vnd.docker.distribution.manifest.v2+json",
])
TIMEOUT = 20
WORKERS = 16


def probe(image):
  """-> None if the image exists, else a message. `image` is the path after the
  registry host, e.g. "general-theiagen/staphb/bwa:0.7.18"."""
  if "@" in image:  # digest pin; the manifest endpoint takes it as a reference
    repo, reference = image.split("@", 1)
  elif ":" in image.rsplit("/", 1)[-1]:
    repo, reference = image.rsplit(":", 1)
  else:
    repo, reference = image, "latest"

  request = urllib.request.Request(
    f"https://{REGISTRY}/v2/{repo}/manifests/{reference}",
    method="HEAD",
    headers={"Accept": ACCEPT},
  )
  try:
    urllib.request.urlopen(request, timeout=TIMEOUT).close()
    return None
  except urllib.error.HTTPError as err:
    if err.code == 404:
      return "not found in registry"
    if err.code in (401, 403):
      return f"warning: HTTP {err.code} (not publicly readable)"
    return f"warning: HTTP {err.code}"
  except Exception as err:  # DNS, timeout, TLS, connection reset
    return f"warning: {type(err).__name__}: {err}"


def collect(paths):
  """-> {image: [(path, line), ...]} for every GAR image reference."""
  found = {}
  for path in paths:
    text = open(path, encoding="utf-8").read()
    for match in GAR_IMAGE.finditer(text):
      line = text[: match.start()].count("\n") + 1
      found.setdefault(match.group(1), []).append((path, line))
  return found


def main(paths):
  images = collect(paths)
  if not images:
    return 0

  with concurrent.futures.ThreadPoolExecutor(WORKERS) as pool:
    problems = dict(zip(images, pool.map(probe, images)))

  failed = False
  for image, problem in sorted(problems.items()):
    if problem is None:
      continue
    failed = failed or not problem.startswith("warning")
    for path, line in images[image]:
      print(f"{path}:{line}: {REGISTRY}/{image}: {problem}", file=sys.stderr)
  return 1 if failed else 0


if __name__ == "__main__":
  raise SystemExit(main(sys.argv[1:]))
