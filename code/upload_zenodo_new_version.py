"""Create a new Zenodo version superseding 10.5281/zenodo.22746385 (requires ZENODO_TOKEN)."""

from __future__ import annotations

import json
import os
import sys
import zipfile
from pathlib import Path

import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
CONCEPT_OR_RECORD_ID = "22746385"  # current version record id
API = "https://zenodo.org/api"


def req(method: str, url: str, token: str, data: dict | None = None, raw: bytes | None = None, content_type: str | None = None):
    headers = {"Authorization": f"Bearer {token}"}
    body = None
    if data is not None:
        body = json.dumps(data).encode("utf-8")
        headers["Content-Type"] = "application/json"
    if raw is not None:
        body = raw
        if content_type:
            headers["Content-Type"] = content_type
    request = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request) as resp:
            payload = resp.read().decode("utf-8")
            return json.loads(payload) if payload else {}
    except urllib.error.HTTPError as e:
        err = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"{method} {url} failed: {e.code} {err}") from e


def make_zip(out_zip: Path) -> Path:
    if out_zip.exists():
        out_zip.unlink()
    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in ROOT.rglob("*"):
            if path.is_file() and ".git" not in path.parts:
                if path.name.endswith("_1.py") or path.name.endswith("_1.md"):
                    continue
                zf.write(path, arcname=str(path.relative_to(ROOT.parent / "corrected_verified_deposit_v2").as_posix() if False else path.relative_to(ROOT)))
                # store relative to archive root
    # rewrite cleanly
    out_zip.unlink(missing_ok=True)
    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in ROOT.rglob("*"):
            if not path.is_file():
                continue
            if any(part.startswith(".") for part in path.parts):
                continue
            zf.write(path, arcname=path.relative_to(ROOT).as_posix())
    return out_zip


def main() -> None:
    token = os.environ.get("ZENODO_TOKEN") or os.environ.get("ZENODO_ACCESS_TOKEN")
    if not token:
        print("Set ZENODO_TOKEN first.", file=sys.stderr)
        sys.exit(1)

    zip_path = ROOT / "ad-gamma-digital-twin-corrected_v2.0.0.zip"
    make_zip(zip_path)
    print(f"Built {zip_path} ({zip_path.stat().st_size} bytes)")

    # Start new version from existing record
    newver = req("POST", f"{API}/deposit/depositions/{CONCEPT_OR_RECORD_ID}/actions/newversion", token)
    latest_url = newver["links"]["latest_draft"]
    draft = req("GET", latest_url, token)
    dep_id = draft["id"]
    bucket = draft["links"]["bucket"]
    print(f"Draft deposition id={dep_id}")

    # Delete old files in draft
    for f in draft.get("files", []):
        req("DELETE", f["links"]["self"], token)
        print(f"Removed old file {f.get('filename')}")

    # Upload new zip
    with open(zip_path, "rb") as fh:
        raw = fh.read()
    upload_url = f"{bucket}/{zip_path.name}"
    req("PUT", upload_url, token, raw=raw, content_type="application/octet-stream")
    print(f"Uploaded {zip_path.name}")

    metadata = {
        "metadata": {
            "title": "Corrected verified archive for population-level gamma-entrainment digital twin analyses (v2.0.0)",
            "upload_type": "dataset",
            "description": (
                "<p><strong>CORRECTED archive superseding v1.0.0. Do not use v1.0.0.</strong></p>"
                "<p>Raw published Etter/Iaccarino values with provenance; analysis specification "
                "separated from results; Etter n=7 inclusion coded; claim boundary limited to "
                "frequency-titration consistency (not personalization validation).</p>"
            ),
            "creators": [
                {"name": "Sanker, Vivek", "affiliation": "Stanford University", "orcid": "0000-0003-0615-8397"},
                {"name": "Awasthi, Shriya", "affiliation": "Stanford University"},
                {"name": "Desai, Atman", "affiliation": "Stanford University"},
                {"name": "Sha, Sharon", "affiliation": "Stanford University"},
            ],
            "version": "2.0.0",
            "license": "mit",
            "keywords": ["Alzheimer", "gamma entrainment", "digital twin", "correction", "reproducibility"],
        }
    }
    req("PUT", f"{API}/deposit/depositions/{dep_id}", token, data=metadata)
    published = req("POST", f"{API}/deposit/depositions/{dep_id}/actions/publish", token)
    print(json.dumps({"doi": published.get("doi"), "links": published.get("links", {})}, indent=2))


if __name__ == "__main__":
    main()
