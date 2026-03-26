#!/usr/bin/env python3
"""Run native paired-end TAV taxonomy assignment via Java NB core extension."""

import argparse
import json
import os
import shlex
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for native paired-end TAV classification."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="R1 input FASTA/FASTQ or interleaved input file.")
    parser.add_argument("--input2", help="R2 input FASTA/FASTQ file (required unless --interleaved).")
    parser.add_argument(
        "--interleaved",
        action="store_true",
        help="Interpret --input as one interleaved paired stream.",
    )
    parser.add_argument("--output", required=True, help="Output taxonomy file (single TAV-level output).")
    parser.add_argument(
        "--rdp-root",
        default="data/third_party/rdp_classifier",
        help="RDP root directory containing manifest.json.",
    )
    parser.add_argument("--rdp-jar", help="Override path to stock classifier.jar.")
    parser.add_argument("--train-prop", help="Override path to pretrained rRNAClassifier.properties.")
    parser.add_argument(
        "--gene",
        default="16srrna",
        help="Gene model when --train-prop is not set (default: 16srrna).",
    )
    parser.add_argument(
        "--format",
        default="allrank",
        choices=["allrank", "fixrank", "filterbyconf", "db"],
        help="Stock-like output format.",
    )
    parser.add_argument(
        "--conf",
        type=float,
        default=0.8,
        help="Confidence cutoff for filterbyconf/db output semantics.",
    )
    parser.add_argument(
        "--min-words",
        type=int,
        default=5,
        help="Minimum bootstrap words (stock RDP semantics, minimum 5).",
    )
    parser.add_argument(
        "--shortseq-outfile",
        help="Optional output file for short/unclassifiable paired query IDs.",
    )
    parser.add_argument(
        "--java-opt",
        action="append",
        default=[],
        help="Extra JVM option. Repeat to pass multiple options.",
    )
    parser.add_argument("--java-bin", default="java", help="Java executable name or path.")
    parser.add_argument("--javac-bin", default="javac", help="Javac executable name or path.")
    return parser.parse_args()


def load_manifest(rdp_root: Path) -> dict:
    """Load RDP downloader manifest from disk."""
    manifest_path = rdp_root / "manifest.json"
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def resolve_rdp_paths(args: argparse.Namespace) -> tuple[Path, str | None]:
    """Resolve classifier jar path and effective training property path."""
    if args.rdp_jar:
        jar_path = Path(args.rdp_jar)
    else:
        manifest = load_manifest(Path(args.rdp_root))
        classifier_root = Path(manifest["paths"]["classifier"])
        jar_path = sorted(classifier_root.rglob("classifier.jar"))[0]

    if args.train_prop:
        train_prop = args.train_prop
    else:
        manifest = load_manifest(Path(args.rdp_root))
        pretrained_root = Path(manifest["paths"]["pretrained_data"])
        prop_pattern = f"classifier/{args.gene}/rRNAClassifier.properties"
        prop_candidates = sorted(pretrained_root.rglob(prop_pattern))
        train_prop = str(prop_candidates[0]) if prop_candidates else None

    return jar_path, train_prop


def collect_runtime_jars(classifier_jar: Path) -> list[str]:
    """Collect runtime classpath jars from stock RDP distribution."""
    lib_dir = classifier_jar.parent.parent / "lib"
    lib_jars = sorted(str(path) for path in lib_dir.glob("*.jar"))
    return [str(classifier_jar)] + lib_jars


def compile_java_extension(
    repo_root: Path,
    classpath_jars: list[str],
    javac_bin: str,
) -> Path:
    """Compile local Java paired-classifier sources into build classes."""
    java_src_dir = repo_root / "java" / "src"
    java_build_dir = repo_root / "java" / "build" / "classes"
    java_build_dir.mkdir(parents=True, exist_ok=True)

    source_files = sorted(str(path) for path in java_src_dir.rglob("*.java"))
    compile_cmd = [
        javac_bin,
        "-cp",
        os.pathsep.join(classpath_jars),
        "-d",
        str(java_build_dir),
        *source_files,
    ]
    subprocess.run(compile_cmd, check=True)
    return java_build_dir


def build_java_command(
    args: argparse.Namespace,
    classpath_entries: list[str],
    train_prop: str | None,
) -> list[str]:
    """Build the Java execution command for paired native classifier."""
    cmd = [args.java_bin]
    for item in args.java_opt:
        for token in shlex.split(item):
            cmd.append(token)

    cmd.extend(
        [
            "-cp",
            os.pathsep.join(classpath_entries),
            "org.vsearchplus.rdp.PairedClassifierMain",
            "--input",
            args.input,
            "--output",
            args.output,
            "--format",
            args.format,
            "--conf",
            str(args.conf),
            "--min-words",
            str(args.min_words),
        ]
    )

    if args.interleaved:
        cmd.append("--interleaved")
    else:
        cmd.extend(["--input2", args.input2])

    if train_prop:
        cmd.extend(["--train-prop", train_prop])
    else:
        cmd.extend(["--gene", args.gene])

    if args.shortseq_outfile:
        cmd.extend(["--shortseq-outfile", args.shortseq_outfile])

    return cmd


def main() -> None:
    """Compile and execute native paired-end TAV classification."""
    args = parse_args()

    if args.interleaved and args.input2:
        raise RuntimeError("--input2 is not allowed when --interleaved is set")
    if (not args.interleaved) and (not args.input2):
        raise RuntimeError("--input2 is required unless --interleaved is set")

    classifier_jar, train_prop = resolve_rdp_paths(args)
    runtime_jars = collect_runtime_jars(classifier_jar)

    repo_root = Path(__file__).resolve().parents[2]
    build_classes = compile_java_extension(repo_root, runtime_jars, args.javac_bin)

    classpath_entries = [str(build_classes)] + runtime_jars
    java_cmd = build_java_command(args, classpath_entries, train_prop)
    subprocess.run(java_cmd, check=True)

    print(f"RDP jar: {classifier_jar}")
    if train_prop:
        print(f"Pretrained model property: {train_prop}")
    else:
        print(f"Gene model: {args.gene}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
