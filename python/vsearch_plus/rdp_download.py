#!/usr/bin/env python3
"""Download the latest RDP Classifier release and pretrained assets from SourceForge."""

import argparse
import json
import re
import subprocess
import tarfile
import zipfile
from pathlib import Path
from urllib.request import Request, urlopen


BASE_URL = "https://sourceforge.net/projects/rdp-classifier/files/"
CLASSIFIER_PAGE_URL = BASE_URL + "rdp-classifier/"
TRAINING_PAGE_URL = BASE_URL + "RDP_Classifier_TrainingData/"


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
        "--output-root",
        default="data/third_party/rdp_classifier",
        help="Root output directory for downloaded/extracted files.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download and re-extract.",
    )
    return parser.parse_args()


def main() -> None:
    """Entrypoint for downloading and extracting latest RDP assets."""
    args = parse_args()
    root_dir = Path(args.output_root)
    downloads_dir = root_dir / "downloads"
    extracted_dir = root_dir / "extracted"
    downloads_dir.mkdir(parents=True, exist_ok=True)
    extracted_dir.mkdir(parents=True, exist_ok=True)

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

    classifier_zip_path = downloads_dir / classifier_zip_name
    pretrained_data_path = downloads_dir / "data.tgz"
    raw_training_zip_path = downloads_dir / raw_training_zip_name
    qiime_training_zip_path = downloads_dir / qiime_training_zip_name

    download_file(classifier_url, classifier_zip_path, args.force)
    download_file(pretrained_data_url, pretrained_data_path, args.force)
    download_file(raw_training_url, raw_training_zip_path, args.force)
    download_file(qiime_training_url, qiime_training_zip_path, args.force)

    classifier_extract_dir = extracted_dir / f"rdp_classifier_{classifier_version}"
    pretrained_extract_dir = extracted_dir / "data"
    raw_training_extract_dir = extracted_dir / raw_training_zip_name.replace(".zip", "")
    qiime_training_extract_dir = extracted_dir / qiime_training_zip_name.replace(".zip", "")

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
            "root": str(root_dir),
            "downloads": str(downloads_dir),
            "extracted": str(extracted_dir),
            "classifier": str(classifier_extract_dir),
            "pretrained_data": str(pretrained_extract_dir),
            "raw_training": str(raw_training_extract_dir),
            "qiime_training": str(qiime_training_extract_dir),
        },
        "urls": {
            "classifier": classifier_url,
            "pretrained_data": pretrained_data_url,
            "raw_training": raw_training_url,
            "qiime_training": qiime_training_url,
        },
    }

    manifest_path = root_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print("\nDownload and extraction complete.")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
