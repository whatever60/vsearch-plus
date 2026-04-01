#!/usr/bin/env python3
"""Download the latest RDP Classifier code and pretrained assets into the runtime root."""

import argparse
import json
import os
import re
import subprocess
import tarfile
import zipfile
from pathlib import Path
from urllib.request import Request, urlopen


BASE_URL = "https://sourceforge.net/projects/rdp-classifier/files/"
CLASSIFIER_PAGE_URL = BASE_URL + "rdp-classifier/"
TRAINING_PAGE_URL = BASE_URL + "RDP_Classifier_TrainingData/"
PACKAGE_ROOT = Path(__file__).resolve().parent.parent


def fetch_html(url: str) -> str:
    """Fetch an HTML page as UTF-8 text."""
    request = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(request) as response:
        return response.read().decode("utf-8", errors="replace")


def version_tuple(version_text: str) -> tuple[int, ...]:
    """Convert semantic version text to a tuple for sorting."""
    return tuple(int(part) for part in version_text.split("."))


def pick_latest_classifier_filename(html: str) -> str:
    """Pick the latest classifier zip filename from a SourceForge listing page."""
    versions = sorted(
        set(re.findall(r"rdp_classifier_(\d+(?:\.\d+)*)\.zip", html)),
        key=version_tuple,
    )
    latest_version = versions[-1]
    return f"rdp_classifier_{latest_version}.zip"


def pick_latest_training_filename(html: str, suffix: str) -> str:
    """Pick the latest numbered training filename for a given suffix."""
    pattern = rf"RDPClassifier_16S_trainsetNo(\d+)_{suffix}\.zip"
    numbers = sorted(set(re.findall(pattern, html)), key=int)
    latest_number = numbers[-1]
    return f"RDPClassifier_16S_trainsetNo{latest_number}_{suffix}.zip"


def download_file(url: str, destination: Path, force: bool) -> None:
    """Download a URL to destination unless file exists and force is false."""
    if destination.exists() and (not force):
        print(f"Using existing file: {destination}")
        return

    print(f"Downloading {url}")
    subprocess.run(
        [
            "wget",
            "-q",
            "-O",
            str(destination),
            url,
        ],
        check=True,
    )


def extract_zip_archive(archive_path: Path, destination: Path, force: bool) -> None:
    """Extract a zip archive into destination."""
    if destination.exists() and (not force):
        print(f"Using existing extracted folder: {destination}")
        return

    if destination.exists() and force:
        for path in sorted(destination.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            else:
                path.rmdir()

    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path, "r") as archive:
        archive.extractall(destination)


def extract_tgz_archive(archive_path: Path, destination: Path, force: bool) -> None:
    """Extract a tar.gz archive into destination."""
    if destination.exists() and (not force):
        print(f"Using existing extracted folder: {destination}")
        return

    if destination.exists() and force:
        for path in sorted(destination.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            else:
                path.rmdir()

    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "r:gz") as archive:
        archive.extractall(destination)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download and re-extract.",
    )
    return parser.parse_args()


def resolve_runtime_root() -> Path:
    """Resolve the runtime root for downloaded assets and generated files."""
    if "VSEARCH_PLUS_RUNTIME_ROOT" in os.environ:
        return Path(os.environ["VSEARCH_PLUS_RUNTIME_ROOT"]).resolve()

    if (PACKAGE_ROOT / "cpp").exists() and (PACKAGE_ROOT / "java").exists():
        return PACKAGE_ROOT.resolve()

    if "XDG_DATA_HOME" in os.environ:
        return (Path(os.environ["XDG_DATA_HOME"]) / "vsearch-plus").resolve()

    return (Path.home() / ".local" / "share" / "vsearch-plus").resolve()


def manifest_path_text(runtime_root: Path, path: Path) -> str:
    """Render a manifest path relative to the runtime root when possible."""
    resolved_path = path.resolve()
    resolved_runtime_root = runtime_root.resolve()
    if resolved_path.is_relative_to(resolved_runtime_root):
        return resolved_path.relative_to(resolved_runtime_root).as_posix()
    return str(resolved_path)


def main() -> None:
    """Entrypoint for downloading and extracting latest RDP assets."""
    args = parse_args()
    runtime_root = resolve_runtime_root()
    data_root = runtime_root / "data" / "rdp_classifier"
    data_downloads_dir = data_root / "downloads"
    data_extract_dir = data_root / "extracted"
    java_root = runtime_root / "extern" / "java" / "rdp_classifier"
    java_downloads_dir = java_root / "downloads"
    java_extract_dir = java_root / "extracted"
    data_downloads_dir.mkdir(parents=True, exist_ok=True)
    data_extract_dir.mkdir(parents=True, exist_ok=True)
    java_downloads_dir.mkdir(parents=True, exist_ok=True)
    java_extract_dir.mkdir(parents=True, exist_ok=True)

    classifier_html = fetch_html(CLASSIFIER_PAGE_URL)
    training_html = fetch_html(TRAINING_PAGE_URL)

    classifier_zip_name = pick_latest_classifier_filename(classifier_html)
    classifier_version = classifier_zip_name.replace("rdp_classifier_", "").replace(".zip", "")
    raw_training_zip_name = pick_latest_training_filename(training_html, "rawtrainingdata")
    qiime_training_zip_name = pick_latest_training_filename(training_html, "QiimeFormat")

    classifier_url = BASE_URL + f"rdp-classifier/{classifier_zip_name}/download"
    pretrained_data_url = BASE_URL + "rdp-classifier/data.tgz/download"
    raw_training_url = BASE_URL + f"RDP_Classifier_TrainingData/{raw_training_zip_name}/download"
    qiime_training_url = BASE_URL + f"RDP_Classifier_TrainingData/{qiime_training_zip_name}/download"

    classifier_zip_path = java_downloads_dir / classifier_zip_name
    pretrained_data_path = data_downloads_dir / "data.tgz"
    raw_training_zip_path = data_downloads_dir / raw_training_zip_name
    qiime_training_zip_path = data_downloads_dir / qiime_training_zip_name

    download_file(classifier_url, classifier_zip_path, args.force)
    download_file(pretrained_data_url, pretrained_data_path, args.force)
    download_file(raw_training_url, raw_training_zip_path, args.force)
    download_file(qiime_training_url, qiime_training_zip_path, args.force)

    classifier_extract_dir = java_extract_dir / f"rdp_classifier_{classifier_version}"
    pretrained_extract_dir = data_extract_dir / "data"
    raw_training_extract_dir = data_extract_dir / raw_training_zip_name.replace(".zip", "")
    qiime_training_extract_dir = data_extract_dir / qiime_training_zip_name.replace(".zip", "")

    extract_zip_archive(classifier_zip_path, classifier_extract_dir, args.force)
    extract_tgz_archive(pretrained_data_path, pretrained_extract_dir, args.force)
    extract_zip_archive(raw_training_zip_path, raw_training_extract_dir, args.force)
    extract_zip_archive(qiime_training_zip_path, qiime_training_extract_dir, args.force)

    manifest = {
        "classifier_zip": classifier_zip_name,
        "classifier_version": classifier_version,
        "raw_training_zip": raw_training_zip_name,
        "qiime_training_zip": qiime_training_zip_name,
        "paths": {
            "root": manifest_path_text(runtime_root, data_root),
            "downloads": manifest_path_text(runtime_root, data_downloads_dir),
            "extracted": manifest_path_text(runtime_root, data_extract_dir),
            "data_root": manifest_path_text(runtime_root, data_root),
            "java_root": manifest_path_text(runtime_root, java_root),
            "java_downloads": manifest_path_text(runtime_root, java_downloads_dir),
            "java_extracted": manifest_path_text(runtime_root, java_extract_dir),
            "classifier": manifest_path_text(runtime_root, classifier_extract_dir),
            "pretrained_data": manifest_path_text(runtime_root, pretrained_extract_dir),
            "raw_training": manifest_path_text(runtime_root, raw_training_extract_dir),
            "qiime_training": manifest_path_text(runtime_root, qiime_training_extract_dir),
        },
        "urls": {
            "classifier": classifier_url,
            "pretrained_data": pretrained_data_url,
            "raw_training": raw_training_url,
            "qiime_training": qiime_training_url,
        },
    }

    manifest_path = data_root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print("\nDownload and extraction complete.")
    print(f"Data root: {data_root}")
    print(f"Java root: {java_root}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
