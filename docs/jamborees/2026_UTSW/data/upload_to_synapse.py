"""Idempotently upload a local file to Synapse and return the resulting syn ID.

Used to mirror jamboree artifacts into Synapse on a per-file basis. Folder structure
under the parent (default: syn64423137 = tf_perturb_seq project) is created on demand.

Authentication: reads SYNAPSE_AUTH_TOKEN from the environment. Run via
`zsh -ic '<this script>'` so the token from ~/.zshrc is loaded.

Usage:
    python upload_to_synapse.py \
        --src /tmp/tf_metadata.tsv \
        --remote-path 2026_UTSW/reference/tf_metadata.tsv \
        --description "Comprehensive TF metadata table for the 2026 UTSW jamboree."
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import synapseclient
from synapseclient import File, Folder

DEFAULT_PARENT = "syn64423137"


def get_or_create_folder(syn: synapseclient.Synapse, name: str, parent_id: str) -> str:
    for child in syn.getChildren(parent_id):
        if child["type"].endswith(".Folder") and child["name"] == name:
            print(f"  using existing folder {child['id']} :: {name}")
            return child["id"]
    folder = syn.store(Folder(name, parent=parent_id))
    print(f"  created folder {folder.id} :: {name}")
    return folder.id


def ensure_folder_path(syn: synapseclient.Synapse, parent_id: str, parts: list[str]) -> str:
    cur = parent_id
    for p in parts:
        cur = get_or_create_folder(syn, p, cur)
    return cur


def upload_file(
    syn: synapseclient.Synapse,
    src: Path,
    parent_id: str,
    description: str | None = None,
) -> str:
    f = File(str(src), parent=parent_id, description=description)
    stored = syn.store(f)
    print(f"  uploaded file {stored.id} :: {src.name}")
    return stored.id


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, type=Path, help="Local path to upload")
    ap.add_argument("--remote-path", required=True, help="Slash-separated path under the parent (e.g., 2026_UTSW/reference/tf_metadata.tsv)")
    ap.add_argument("--parent", default=DEFAULT_PARENT, help=f"Synapse parent ID (default {DEFAULT_PARENT})")
    ap.add_argument("--description", default=None, help="File description (optional)")
    args = ap.parse_args()

    if not args.src.exists():
        sys.exit(f"src not found: {args.src}")

    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        sys.exit("SYNAPSE_AUTH_TOKEN not set in env (run via `zsh -ic` to inherit ~/.zshrc)")

    syn = synapseclient.Synapse(silent=True)
    syn.login(authToken=token)
    print(f"Logged in as: {syn.getUserProfile()['userName']}")

    parts = args.remote_path.strip("/").split("/")
    if len(parts) < 1:
        sys.exit("--remote-path must include at least a filename")
    folder_parts, filename = parts[:-1], parts[-1]
    if filename != args.src.name:
        print(f"NOTE: remote filename '{filename}' differs from local '{args.src.name}'; using local name on upload.")
        filename = args.src.name

    target_parent = ensure_folder_path(syn, args.parent, folder_parts)
    syn_id = upload_file(syn, args.src, target_parent, description=args.description)

    full_remote = "/".join(folder_parts + [args.src.name])
    print()
    print(f"DONE  syn_id={syn_id}  remote_path={full_remote}  url=https://www.synapse.org/Synapse:{syn_id}")


if __name__ == "__main__":
    main()
