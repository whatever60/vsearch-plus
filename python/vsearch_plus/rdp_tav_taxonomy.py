#!/usr/bin/env python3
"""Run native paired-end TAV taxonomy assignment using the RDP classifier."""

import argparse
import csv
import json
import re
import shlex
import subprocess
import tempfile
from pathlib import Path


RANK_ORDER = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIX = {
    "domain": "d",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
    "species": "s",
}
SIZE_PATTERN = re.compile(r"(?:^|;)size=(\d+)(?:;|$)")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments for TAV taxonomy assignment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-fasta", help="Left-anchor TAV FASTA input.")
    parser.add_argument("--right-fasta", help="Right-anchor TAV FASTA input.")
    parser.add_argument(
        "--tav-catalog",
        help="TAV catalog TSV input with tav_id, abundance, left_anchor, right_anchor columns.",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for outputs, for example data/real_out2/tav_taxonomy.",
    )
    parser.add_argument(
        "--rdp-root",
        default="data/third_party/rdp_classifier",
        help="RDP root folder containing manifest.json from python/get_rdp_classifier.py.",
    )
    parser.add_argument("--rdp-jar", help="Override path to classifier.jar.")
    parser.add_argument(
        "--train-prop",
        help="Override path to rRNAClassifier.properties to force pretrained model selection.",
    )
    parser.add_argument(
        "--gene",
        default="16srrna",
        help="RDP gene name used only when --train-prop is not provided.",
    )
    parser.add_argument(
        "--min-conf",
        type=float,
        default=0.8,
        help="Confidence threshold used for paired aggregation.",
    )
    parser.add_argument(
        "--pair-filter",
        choices=["any", "both"],
        default="both",
        help="Aggregation rule: any keeps one end if confident, both requires confident agreement.",
    )
    parser.add_argument(
        "--min-words",
        type=int,
        default=5,
        help="RDP bootstrap minimum words (maps to classifier -w).",
    )
    parser.add_argument(
        "--java-opt",
        action="append",
        default=[],
        help="Extra JVM option. Repeat this flag to add more options.",
    )
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="Keep generated left/right FASTA files when using --tav-catalog.",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read FASTA records as (header_without_gt, sequence) tuples."""
    records: list[tuple[str, str]] = []
    header = ""
    seq_chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        if line[0] == ">":
            if header:
                records.append((header, "".join(seq_chunks)))
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
    if header:
        records.append((header, "".join(seq_chunks)))
    return records


def extract_tav_id(seq_id: str) -> str:
    """Extract canonical TAV ID from a FASTA sequence identifier token."""
    return seq_id.split()[0].split(";")[0]


def parse_size(seq_id: str) -> int:
    """Parse abundance from ';size=N' in sequence ID and default to 1."""
    match = SIZE_PATTERN.search(seq_id)
    if match:
        return int(match.group(1))
    return 1


def build_records_from_pair_fastas(left_fasta: Path, right_fasta: Path) -> list[dict[str, object]]:
    """Build paired TAV records from aligned left and right FASTA files."""
    left_records = read_fasta(left_fasta)
    right_records = read_fasta(right_fasta)
    if len(left_records) != len(right_records):
        raise RuntimeError("Left and right FASTA record counts differ")

    records: list[dict[str, object]] = []
    for index, left_record in enumerate(left_records):
        right_record = right_records[index]
        left_header = left_record[0]
        right_header = right_record[0]
        left_seq = left_record[1]
        right_seq = right_record[1]

        left_seq_id = left_header.split()[0]
        right_seq_id = right_header.split()[0]
        left_tav_id = extract_tav_id(left_seq_id)
        right_tav_id = extract_tav_id(right_seq_id)
        if left_tav_id != right_tav_id:
            raise RuntimeError(
                f"TAV ID mismatch at row {index + 1}: {left_tav_id} vs {right_tav_id}"
            )

        abundance = parse_size(left_seq_id)
        records.append(
            {
                "tav_id": left_tav_id,
                "abundance": abundance,
                "left_seq_id": left_seq_id,
                "right_seq_id": right_seq_id,
                "left_seq": left_seq,
                "right_seq": right_seq,
            }
        )
    return records


def build_records_from_catalog(catalog_path: Path) -> list[dict[str, object]]:
    """Build paired TAV records from a TAV catalog TSV."""
    records: list[dict[str, object]] = []
    with catalog_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tav_id = row["tav_id"]
            abundance = int(row["abundance"])
            left_seq = row["left_anchor"]
            right_seq = row["right_anchor"]
            seq_id = f"{tav_id};size={abundance}"
            records.append(
                {
                    "tav_id": tav_id,
                    "abundance": abundance,
                    "left_seq_id": seq_id,
                    "right_seq_id": seq_id,
                    "left_seq": left_seq,
                    "right_seq": right_seq,
                }
            )
    return records


def write_fasta(records: list[dict[str, object]], side: str, output_path: Path) -> None:
    """Write one side of paired records to FASTA for RDP classification."""
    lines: list[str] = []
    for record in records:
        seq_id_key = f"{side}_seq_id"
        seq_key = f"{side}_seq"
        lines.append(f">{record[seq_id_key]}")
        lines.append(str(record[seq_key]))
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_manifest(rdp_root: Path) -> dict:
    """Load the downloader manifest from an RDP root folder."""
    manifest_path = rdp_root / "manifest.json"
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def resolve_rdp_inputs(args: argparse.Namespace) -> tuple[Path, str | None]:
    """Resolve classifier jar path and optional training properties path."""
    if args.rdp_jar:
        jar_path = Path(args.rdp_jar)
    else:
        manifest = load_manifest(Path(args.rdp_root))
        classifier_root = Path(manifest["paths"]["classifier"])
        jar_candidates = sorted(classifier_root.rglob("classifier.jar"))
        jar_path = jar_candidates[0]

    if args.train_prop:
        train_prop = args.train_prop
    else:
        manifest = load_manifest(Path(args.rdp_root))
        pretrained_root = Path(manifest["paths"]["pretrained_data"])
        pattern = f"classifier/{args.gene}/rRNAClassifier.properties"
        prop_candidates = sorted(pretrained_root.rglob(pattern))
        if prop_candidates:
            train_prop = str(prop_candidates[0])
        else:
            train_prop = None

    return jar_path, train_prop


def build_java_cmd(args: argparse.Namespace, jar_path: Path) -> list[str]:
    """Build the JVM command prefix including extra JVM options."""
    cmd = ["java"]
    for item in args.java_opt:
        for token in shlex.split(item):
            cmd.append(token)
    cmd.extend(["-jar", str(jar_path), "classify"])
    return cmd


def run_rdp(
    args: argparse.Namespace,
    jar_path: Path,
    train_prop: str | None,
    input_fasta: Path,
    output_tsv: Path,
) -> None:
    """Run one RDP classification pass for one FASTA file."""
    cmd = build_java_cmd(args, jar_path)
    cmd.extend(["-f", "allrank", "-w", str(args.min_words), "-o", str(output_tsv)])
    if train_prop:
        cmd.extend(["-t", train_prop])
    else:
        cmd.extend(["-g", args.gene])
    cmd.append(str(input_fasta))
    subprocess.run(cmd, check=True)


def parse_allrank_tsv(path: Path) -> dict[str, dict[str, dict[str, float | str]]]:
    """Parse RDP allrank output into seq_id -> rank -> {taxon, conf}."""
    parsed: dict[str, dict[str, dict[str, float | str]]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        parts = line.split("\t")
        seq_id = parts[0]
        rank_map: dict[str, dict[str, float | str]] = {}
        index = 2
        while index + 2 < len(parts):
            taxon = parts[index]
            rank = parts[index + 1].lower()
            conf = float(parts[index + 2])
            rank_map[rank] = {"taxon": taxon, "conf": conf}
            index += 3
        parsed[seq_id] = rank_map
    return parsed


def complete_rank_map(rank_map: dict[str, dict[str, float | str]]) -> dict[str, dict[str, float | str]]:
    """Fill missing ranks with blank taxonomy and zero confidence."""
    complete: dict[str, dict[str, float | str]] = {}
    for rank in RANK_ORDER:
        if rank in rank_map:
            complete[rank] = {
                "taxon": str(rank_map[rank]["taxon"]),
                "conf": float(rank_map[rank]["conf"]),
            }
        else:
            complete[rank] = {"taxon": "", "conf": 0.0}
    return complete


def aggregate_rank(
    left_rank: dict[str, float | str],
    right_rank: dict[str, float | str],
    min_conf: float,
    pair_filter: str,
) -> dict[str, float | str]:
    """Aggregate one taxonomic rank from left and right end assignments."""
    left_taxon = str(left_rank["taxon"])
    right_taxon = str(right_rank["taxon"])
    left_conf = float(left_rank["conf"])
    right_conf = float(right_rank["conf"])
    left_pass = left_taxon != "" and left_conf >= min_conf
    right_pass = right_taxon != "" and right_conf >= min_conf

    if pair_filter == "both":
        if left_pass and right_pass and left_taxon == right_taxon:
            return {"taxon": left_taxon, "conf": min(left_conf, right_conf), "source": "both"}
        return {"taxon": "", "conf": 0.0, "source": "none"}

    if left_pass and right_pass:
        if left_taxon == right_taxon:
            return {"taxon": left_taxon, "conf": max(left_conf, right_conf), "source": "both"}
        if left_conf >= right_conf:
            return {"taxon": left_taxon, "conf": left_conf, "source": "left_conflict"}
        return {"taxon": right_taxon, "conf": right_conf, "source": "right_conflict"}

    if left_pass:
        return {"taxon": left_taxon, "conf": left_conf, "source": "left"}
    if right_pass:
        return {"taxon": right_taxon, "conf": right_conf, "source": "right"}
    return {"taxon": "", "conf": 0.0, "source": "none"}


def aggregate_lineage(
    left_map: dict[str, dict[str, float | str]],
    right_map: dict[str, dict[str, float | str]],
    min_conf: float,
    pair_filter: str,
) -> dict[str, dict[str, float | str]]:
    """Aggregate complete lineage across ranks with prefix-consistent truncation."""
    paired: dict[str, dict[str, float | str]] = {}
    blocked = False
    for rank in RANK_ORDER:
        if blocked:
            paired[rank] = {"taxon": "", "conf": 0.0, "source": "none"}
            continue

        rank_result = aggregate_rank(left_map[rank], right_map[rank], min_conf, pair_filter)
        paired[rank] = rank_result
        if str(rank_result["taxon"]) == "":
            blocked = True
    return paired


def lineage_string(rank_map: dict[str, dict[str, float | str]]) -> str:
    """Convert a rank map into semicolon-delimited lineage text."""
    pieces: list[str] = []
    for rank in RANK_ORDER:
        taxon = str(rank_map[rank]["taxon"])
        if taxon == "":
            break
        pieces.append(f"{RANK_PREFIX[rank]}__{taxon}")
    if pieces:
        return ";".join(pieces)
    return "Unclassified"


def write_paired_table(
    records: list[dict[str, object]],
    left_results: dict[str, dict[str, dict[str, float | str]]],
    right_results: dict[str, dict[str, dict[str, float | str]]],
    output_path: Path,
    min_conf: float,
    pair_filter: str,
) -> list[dict[str, object]]:
    """Write paired taxonomy table and return cached paired lineage records."""
    paired_rows: list[dict[str, object]] = []
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        header = [
            "tav_id",
            "abundance",
            "left_seq_id",
            "right_seq_id",
            "left_lineage",
            "right_lineage",
            "paired_lineage",
        ]
        for rank in RANK_ORDER:
            header.extend(
                [
                    f"left_{rank}",
                    f"left_{rank}_conf",
                    f"right_{rank}",
                    f"right_{rank}_conf",
                    f"paired_{rank}",
                    f"paired_{rank}_conf",
                    f"paired_{rank}_source",
                ]
            )
        writer.writerow(header)

        for record in records:
            left_seq_id = str(record["left_seq_id"])
            right_seq_id = str(record["right_seq_id"])
            if left_seq_id not in left_results:
                raise RuntimeError(f"Missing left RDP result for {left_seq_id}")
            if right_seq_id not in right_results:
                raise RuntimeError(f"Missing right RDP result for {right_seq_id}")

            left_map = complete_rank_map(left_results[left_seq_id])
            right_map = complete_rank_map(right_results[right_seq_id])
            paired_map = aggregate_lineage(left_map, right_map, min_conf, pair_filter)
            left_lineage = lineage_string(left_map)
            right_lineage = lineage_string(right_map)
            paired_lineage = lineage_string(paired_map)

            row = [
                str(record["tav_id"]),
                str(record["abundance"]),
                left_seq_id,
                right_seq_id,
                left_lineage,
                right_lineage,
                paired_lineage,
            ]
            for rank in RANK_ORDER:
                row.extend(
                    [
                        str(left_map[rank]["taxon"]),
                        f"{float(left_map[rank]['conf']):.4f}",
                        str(right_map[rank]["taxon"]),
                        f"{float(right_map[rank]['conf']):.4f}",
                        str(paired_map[rank]["taxon"]),
                        f"{float(paired_map[rank]['conf']):.4f}",
                        str(paired_map[rank]["source"]),
                    ]
                )
            writer.writerow(row)

            paired_rows.append(
                {
                    "tav_id": str(record["tav_id"]),
                    "abundance": int(record["abundance"]),
                    "paired_map": paired_map,
                }
            )

    return paired_rows


def write_rank_counts(paired_rows: list[dict[str, object]], output_path: Path) -> None:
    """Write aggregated counts by rank and assigned taxon."""
    counts: dict[tuple[str, str], dict[str, int]] = {}
    for row in paired_rows:
        abundance = int(row["abundance"])
        paired_map = row["paired_map"]
        for rank in RANK_ORDER:
            taxon = str(paired_map[rank]["taxon"])
            if taxon == "":
                continue
            key = (rank, taxon)
            if key not in counts:
                counts[key] = {"tav_count": 0, "total_abundance": 0}
            counts[key]["tav_count"] += 1
            counts[key]["total_abundance"] += abundance

    rank_index = {rank: idx for idx, rank in enumerate(RANK_ORDER)}
    sorted_keys = sorted(
        counts,
        key=lambda key: (rank_index[key[0]], -counts[key]["total_abundance"], key[1]),
    )

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "taxon", "tav_count", "total_abundance"])
        for key in sorted_keys:
            writer.writerow([
                key[0],
                key[1],
                str(counts[key]["tav_count"]),
                str(counts[key]["total_abundance"]),
            ])


def main() -> None:
    """Run the full paired-end TAV taxonomy assignment workflow."""
    args = parse_args()
    if args.tav_catalog:
        records = build_records_from_catalog(Path(args.tav_catalog))
        if args.keep_intermediate:
            intermediate_dir = Path(f"{args.output_prefix}.intermediate")
            intermediate_dir.mkdir(parents=True, exist_ok=True)
            left_input = intermediate_dir / "left_from_catalog.fa"
            right_input = intermediate_dir / "right_from_catalog.fa"
            write_fasta(records, "left", left_input)
            write_fasta(records, "right", right_input)
            temp_context = None
        else:
            temp_context = tempfile.TemporaryDirectory(prefix="vsearch_plus_rdp_")
            temp_dir = Path(temp_context.name)
            left_input = temp_dir / "left_from_catalog.fa"
            right_input = temp_dir / "right_from_catalog.fa"
            write_fasta(records, "left", left_input)
            write_fasta(records, "right", right_input)
    else:
        if not args.left_fasta or not args.right_fasta:
            raise RuntimeError("Either provide --tav-catalog or both --left-fasta and --right-fasta")
        records = build_records_from_pair_fastas(Path(args.left_fasta), Path(args.right_fasta))
        left_input = Path(args.left_fasta)
        right_input = Path(args.right_fasta)
        temp_context = None

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    left_output = Path(f"{output_prefix}.left.allrank.tsv")
    right_output = Path(f"{output_prefix}.right.allrank.tsv")
    paired_output = Path(f"{output_prefix}.paired.tsv")
    counts_output = Path(f"{output_prefix}.rank_counts.tsv")

    jar_path, train_prop = resolve_rdp_inputs(args)

    run_rdp(args, jar_path, train_prop, left_input, left_output)
    run_rdp(args, jar_path, train_prop, right_input, right_output)

    left_results = parse_allrank_tsv(left_output)
    right_results = parse_allrank_tsv(right_output)
    paired_rows = write_paired_table(
        records,
        left_results,
        right_results,
        paired_output,
        args.min_conf,
        args.pair_filter,
    )
    write_rank_counts(paired_rows, counts_output)

    if temp_context is not None:
        temp_context.cleanup()

    print(f"RDP jar: {jar_path}")
    if train_prop:
        print(f"Pretrained model property: {train_prop}")
    else:
        print(f"Gene model: {args.gene}")
    print(f"Input records: {len(records)}")
    print(f"Left output: {left_output}")
    print(f"Right output: {right_output}")
    print(f"Paired output: {paired_output}")
    print(f"Rank counts: {counts_output}")


if __name__ == "__main__":
    main()
