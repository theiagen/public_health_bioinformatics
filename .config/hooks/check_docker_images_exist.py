#!/usr/bin/env python3
"""Confirm every pinned docker image in a WDL actually exists in its registry.

Asks the registry for the image MANIFEST only (a few hundred bytes of JSON) via
the OCI distribution API — nothing is pulled, and no docker daemon is needed.

Exit policy, chosen so this never blocks work the author didn't break:
  404                     -> FAIL. The tag was typo'd or never pushed.
  401/403 we can't clear  -> warn, pass. We cannot tell missing from
                             unauthorized, so we do not guess.
  DNS/timeout/connection  -> warn, pass. Offline commits must still work.

All 161 us-docker.pkg.dev/general-theiagen images resolve anonymously; the
handful on Docker Hub need the bearer-token handshake below.

Coverage limit: quay.io answers 401 for a repo that does not exist, not 404 —
so on quay this hook can confirm an image IS there but can never prove one is
missing. GAR and Docker Hub both return a real 404.
"""
import concurrent.futures
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

# Manifest media types, newest first. Without these the registry returns a
# v1 manifest (or 404s on multi-arch images that only have an index).
ACCEPT = ", ".join([
  "application/vnd.oci.image.index.v1+json",
  "application/vnd.oci.image.manifest.v1+json",
  "application/vnd.docker.distribution.manifest.list.v2+json",
  "application/vnd.docker.distribution.manifest.v2+json",
])
TIMEOUT = 20
WORKERS = 16

# `docker\w*` deliberately also matches `String bwa_docker = "..."` outputs and
# similar. Do not widen it further — an input merely NAMED like a registry
# (docker_registry, docker_user) would then be probed as if it were an image.
DOCKER_STRING = re.compile(r'docker\w*\s*[:=]\s*"([^"]+)"')
PLACEHOLDER = re.compile(r"[~$]\{")


def split_image(image):
  """"us-docker.pkg.dev/foo/bar:1.0" -> ("us-docker.pkg.dev", "foo/bar", "1.0").

  The first path segment is a registry host only if it looks like one (has a
  dot or a port); otherwise the whole string is a Docker Hub repo, and a
  single-segment repo lives under the implicit "library/" namespace.
  """
  head = image.split("/")[0]
  if "/" in image and ("." in head or ":" in head or head == "localhost"):
    host, repo = head, image.split("/", 1)[1]
  else:
    host = "registry-1.docker.io"
    repo = image if "/" in image else f"library/{image}"

  if "@" in repo:  # digest pin; the manifest endpoint takes it in the tag slot
    repo, reference = repo.split("@", 1)
  elif ":" in repo.rsplit("/", 1)[-1]:
    repo, reference = repo.rsplit(":", 1)
  else:
    reference = "latest"
  return host, repo, reference


def _head(url, token=None):
  headers = {"Accept": ACCEPT}
  if token:
    headers["Authorization"] = f"Bearer {token}"
  req = urllib.request.Request(url, method="HEAD", headers=headers)
  return urllib.request.urlopen(req, timeout=TIMEOUT)


def _anonymous_token(challenge):
  """Answer a `WWW-Authenticate: Bearer realm=...,service=...,scope=...` 401.

  This is the standard registry handshake — anonymous only, no credentials are
  read from the environment or ~/.docker/config.json.
  """
  if not challenge.lower().startswith("bearer "):
    return None
  params = dict(re.findall(r'(\w+)="([^"]*)"', challenge))
  realm = params.pop("realm", None)
  if not realm:
    return None
  url = f"{realm}?{urllib.parse.urlencode(params)}"
  with urllib.request.urlopen(url, timeout=TIMEOUT) as resp:
    body = json.load(resp)
  return body.get("token") or body.get("access_token")


def probe(image):
  """-> (status, detail) where status is "ok", "missing", "unknown"."""
  host, repo, reference = split_image(image)
  url = f"https://{host}/v2/{repo}/manifests/{reference}"
  try:
    try:
      _head(url).close()
      return "ok", None
    except urllib.error.HTTPError as err:
      if err.code != 401:
        raise
      token = _anonymous_token(err.headers.get("WWW-Authenticate", ""))
      if not token:
        return "unknown", "registry requires credentials"
      _head(url, token).close()
      return "ok", None
  except urllib.error.HTTPError as err:
    if err.code == 404:
      return "missing", "not found in registry"
    if err.code in (401, 403):
      return "unknown", f"HTTP {err.code} (private or unauthorized)"
    return "unknown", f"HTTP {err.code}"
  except Exception as err:  # DNS, timeout, TLS, connection reset
    return "unknown", f"{type(err).__name__}: {err}"


def collect(paths):
  """-> {image: [(path, line), ...]} for every non-interpolated docker string."""
  found = {}
  for path in paths:
    text = open(path, encoding="utf-8").read()
    for match in DOCKER_STRING.finditer(text):
      image = match.group(1)
      if PLACEHOLDER.search(image):  # `docker: "~{docker}"` forwards the input
        continue
      line = text[: match.start()].count("\n") + 1
      found.setdefault(image, []).append((path, line))
  return found


def main(paths):
  images = collect(paths)
  if not images:
    return 0

  with concurrent.futures.ThreadPoolExecutor(WORKERS) as pool:
    results = dict(zip(images, pool.map(probe, images)))

  failed = False
  for image, (status, detail) in sorted(results.items()):
    if status == "ok":
      continue
    stream = sys.stderr
    for path, line in images[image]:
      if status == "missing":
        failed = True
        print(f"{path}:{line}: docker image `{image}` {detail}", file=stream)
      else:
        print(f"{path}:{line}: warning: could not verify `{image}` — {detail}",
              file=stream)
  return 1 if failed else 0


if __name__ == "__main__":
  raise SystemExit(main(sys.argv[1:]))
