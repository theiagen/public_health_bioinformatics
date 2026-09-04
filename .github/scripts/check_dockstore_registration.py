#!/usr/bin/env python3
"""Verify that every release workflow is registered in Dockstore.

The rule this enforces:

    Every workflows/**/*.wdl OUTSIDE workflows/utilities/ must have a matching
    primaryDescriptorPath entry in .dockstore.yml.

Files under workflows/utilities/ are subworkflows: they are imported by other
workflows and are never registered on their own. Note the rule is
one-directional -- living outside utilities/ means "must be registered", but
living inside it is the only thing that exempts a file. A registered workflow
may still be imported by another workflow (e.g. wf_snippy_variants.wdl), which
is why "is it imported?" cannot be used as the test.

Run locally:
    python3 .github/scripts/check_dockstore_registration.py
"""

import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: pip install pyyaml")

REPO = Path(__file__).resolve().parents[2]
DOCKSTORE = REPO / ".dockstore.yml"
WORKFLOWS = REPO / "workflows"
SUBWORKFLOW_DIR = WORKFLOWS / "utilities"


def main() -> int:
    config = yaml.safe_load(DOCKSTORE.read_text())
    entries = config.get("workflows") or []

    # Dockstore paths are absolute-from-repo-root ("/workflows/x/y.wdl").
    registered = {}
    for entry in entries:
        path = entry.get("primaryDescriptorPath")
        if path:
            registered[path.lstrip("/")] = entry.get("name", "<unnamed>")

    on_disk = {
        p.relative_to(REPO).as_posix()
        for p in WORKFLOWS.rglob("*.wdl")
    }
    subworkflows = {
        p.relative_to(REPO).as_posix()
        for p in SUBWORKFLOW_DIR.rglob("*.wdl")
    }
    release = on_disk - subworkflows

    unregistered = sorted(release - registered.keys())
    dangling = sorted(set(registered) - on_disk)
    registered_subworkflows = sorted(subworkflows & registered.keys())

    if unregistered:
        print("::error::Release workflows missing from .dockstore.yml")
        print(
            f"\n{len(unregistered)} workflow(s) live outside workflows/utilities/ "
            "but are not registered in .dockstore.yml.\n"
            "Add an entry for each, or move it under workflows/utilities/ if it "
            "is a subworkflow:\n"
        )
        for path in unregistered:
            name = Path(path).stem.removeprefix("wf_")
            print(f" - name: {name}_PHB")
            print(f"   subclass: WDL")
            print(f"   primaryDescriptorPath: /{path}")
            print(f"   testParameterFiles:")
            print(f"    - /tests/inputs/empty.json")
        print()

    if dangling:
        print("::error::.dockstore.yml points at files that do not exist")
        print(
            f"\n{len(dangling)} primaryDescriptorPath value(s) do not resolve. "
            "This usually means a workflow was moved or renamed without updating "
            ".dockstore.yml:\n"
        )
        for path in dangling:
            print(f"  {registered[path]}: /{path}")
        print()

    if registered_subworkflows:
        print("::error::Subworkflows must not be registered in Dockstore")
        print(
            f"\n{len(registered_subworkflows)} file(s) under workflows/utilities/ "
            "are registered. Either move the workflow out of utilities/, or drop "
            "its .dockstore.yml entry:\n"
        )
        for path in registered_subworkflows:
            print(f"  {registered[path]}: /{path}")
        print()

    if unregistered or dangling or registered_subworkflows:
        return 1

    print(
        f"OK: {len(release)} release workflow(s) registered, "
        f"{len(subworkflows)} subworkflow(s) correctly excluded."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
