#!/usr/bin/env python
"""Upload/refresh a private Kaggle DATASET from a single file (for shipping precomputed inputs to a kernel).

GEMME tables (12k+ variants x N proteins) are too large to embed in a notebook, so they ride as a Kaggle
dataset the kernel attaches. `dataset_create_new` the first time, `dataset_create_version` after.

  uv run python scripts/kaggle_upload_dataset.py <file> <dataset-slug> ["title"]
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__); return 2
    src = Path(sys.argv[1]); slug = sys.argv[2]; title = sys.argv[3] if len(sys.argv) > 3 else slug
    from kaggle.api.kaggle_api_extended import KaggleApi
    api = KaggleApi(); api.authenticate()
    user = api.get_config_value("username") or "emanueleebrahimi"
    ref = f"{user}/{slug}"
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        shutil.copyfile(src, tdp / src.name)
        (tdp / "dataset-metadata.json").write_text(json.dumps(
            {"title": title, "id": ref, "licenses": [{"name": "CC0-1.0"}]}, indent=2))
        # exists?
        exists = False
        try:
            api.dataset_status(ref); exists = True
        except Exception:
            exists = False
        if exists:
            print(f"updating dataset {ref} ...")
            api.dataset_create_version(str(tdp), version_notes="refresh", dir_mode="zip", quiet=False)
        else:
            print(f"creating dataset {ref} ...")
            api.dataset_create_new(str(tdp), dir_mode="zip", public=False, quiet=False)
    print(f"dataset ref: {ref}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
