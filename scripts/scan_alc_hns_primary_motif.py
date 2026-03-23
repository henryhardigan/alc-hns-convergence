#!/usr/bin/env python3
"""Scan reviewed Swiss-Prot for the Alc/H-NS primary chemistry-aware motif."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


HEADER_RE = re.compile(r"^(?P<db>sp|tr)\|(?P<acc>[^|]+)\|(?P<entry>\S+)")
ORGANISM_RE = re.compile(r"\bOS=(.+?) OX=")
TAXON_RE = re.compile(r"\bOX=(\d+)\b")
GN_RE = re.compile(r"\bGN=([^ ]+)")

PRIMARY_MOTIF_TEXT = (
    "E-X-R-[NT]-R-[RK]-L-[QE]-Q-X-R-[DE]-M-[LI]{1,2}-[AP]-D-[GA]-[IV]-[DE]-[EP]-X-[EK]-[VFL]{2,3}-X{0,1}-N-[QS]-L-A"
)
PRIMARY_PATTERN = (
    r"(?P<p1>E)"
    r"(?P<p2>.)"
    r"(?P<p3>R)"
    r"(?P<p4>[NT])"
    r"(?P<p5>R)"
    r"(?P<p6>[RK])"
    r"(?P<p7>L)"
    r"(?P<p8>[QE])"
    r"(?P<p9>Q)"
    r"(?P<p10>.)"
    r"(?P<p11>R)"
    r"(?P<p12>[DE])"
    r"(?P<p13>M)"
    r"(?P<p14>[LI]{1,2})"
    r"(?P<p15>[AP])"
    r"(?P<p16>D)"
    r"(?P<p17>[GA])"
    r"(?P<p18>[IV])"
    r"(?P<p19>[DE])"
    r"(?P<p20>[EP])"
    r"(?P<p21>.)"
    r"(?P<p22>[EK])"
    r"(?P<p23>[VFL]{2,3})"
    r"(?P<p24>.{0,1})"
    r"(?P<p25>N)"
    r"(?P<p26>[QS])"
    r"(?P<p27>L)"
    r"(?P<p28>A)"
)
PRIMARY_REGEX = re.compile(PRIMARY_PATTERN)


@dataclass(frozen=True)
class Record:
    accession: str
    entry: str
    organism: str
    taxon_id: str
    gene: str
    header: str
    sequence: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fasta-glob",
        default="sprot_chunks/sprot_*.fasta",
        help="Glob for local Swiss-Prot shards.",
    )
    parser.add_argument(
        "--outdir",
        default="results/alc_hns_convergence",
        help="Output directory for tables and summary.",
    )
    return parser.parse_args()


def iter_fasta_records(fasta_glob: str):
    for path in sorted(Path().glob(fasta_glob)):
        header = None
        seq_parts: list[str] = []
        with path.open() as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_parts)
                    header = line[1:]
                    seq_parts = []
                else:
                    seq_parts.append(line.upper())
            if header is not None:
                yield header, "".join(seq_parts)


def parse_record(header: str, sequence: str) -> Record | None:
    match = HEADER_RE.match(header)
    if match is None or match.group("db") != "sp":
        return None
    organism = ORGANISM_RE.search(header)
    taxon = TAXON_RE.search(header)
    gene = GN_RE.search(header)
    return Record(
        accession=match.group("acc"),
        entry=match.group("entry"),
        organism=organism.group(1) if organism else "",
        taxon_id=taxon.group(1) if taxon else "",
        gene=gene.group(1) if gene else "",
        header=header,
        sequence=sequence,
    )


def load_reviewed_records(fasta_glob: str) -> list[Record]:
    records: list[Record] = []
    for header, sequence in iter_fasta_records(fasta_glob):
        record = parse_record(header, sequence)
        if record is not None:
            records.append(record)
    if not records:
        raise SystemExit("No reviewed Swiss-Prot records found.")
    return records


def is_virus(record: Record) -> bool:
    header = record.header.lower()
    return "virus" in header or "phage" in header or "bacteriophage" in header


def is_hns_family(record: Record) -> bool:
    return (
        "DNA-binding protein H-NS" in record.header
        or "DNA-binding protein StpA" in record.header
        or "VicH" in record.header
    )


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    records = load_reviewed_records(args.fasta_glob)
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    viral_hits = 0
    bacterial_hits = 0
    family_hits = 0
    family_size = 0

    for record in records:
        in_family = is_hns_family(record)
        if in_family:
            family_size += 1
        for match in PRIMARY_REGEX.finditer(record.sequence):
            rows.append(
                {
                    "accession": record.accession,
                    "entry": record.entry,
                    "organism": record.organism,
                    "taxon_id": record.taxon_id,
                    "gene": record.gene,
                    "start_1based": match.start() + 1,
                    "end_1based": match.end(),
                    "matched_sequence": match.group(),
                    "p14_li_run": match.group("p14"),
                    "p20_ep": match.group("p20"),
                    "p23_vfl_run": match.group("p23"),
                    "p24_indel": match.group("p24"),
                    "is_hns_family": "yes" if in_family else "no",
                }
            )
            if is_virus(record):
                viral_hits += 1
            else:
                bacterial_hits += 1
            if in_family:
                family_hits += 1

    rows.sort(key=lambda row: (str(row["entry"]), int(row["start_1based"])))
    write_tsv(
        outdir / "primary_motif_hits_reviewed.tsv",
        rows,
        [
            "accession",
            "entry",
            "organism",
            "taxon_id",
            "gene",
            "start_1based",
            "end_1based",
            "matched_sequence",
            "p14_li_run",
            "p20_ep",
            "p23_vfl_run",
            "p24_indel",
            "is_hns_family",
        ],
    )

    summary = outdir / "primary_motif_summary.md"
    summary.write_text(
        "\n".join(
            [
                "# Alc/H-NS Primary Motif Scan Summary",
                "",
                "## Primary chemistry-aware motif",
                "",
                f"`{PRIMARY_MOTIF_TEXT}`",
                "",
                "Interpretive notes:",
                "",
                "- `X` denotes an arbitrary mismatched residue rather than a gap.",
                "- `[LI]{1,2}` and `[VFL]{2,3}` encode short hydrophobic runs instead of arbitrary spacers.",
                "- `X{0,1}` denotes the sole true indel tolerance in the motif.",
                "",
                "## Hit counts",
                "",
                f"- Reviewed Swiss-Prot hits: {len(rows)}",
                f"- Viral hits: {viral_hits}",
                f"- Bacterial hits: {bacterial_hits}",
                f"- Reviewed H-NS/StpA family proteins considered: {family_size}",
                f"- Exact primary motif hits within reviewed H-NS/StpA family: {family_hits}",
                "",
                "## Hits",
                "",
            ]
            + [
                f"- {row['entry']} ({row['accession']}) {row['start_1based']}-{row['end_1based']}: "
                f"`{row['matched_sequence']}` "
                f"[LI-run={row['p14_li_run']}; p20={row['p20_ep']}; VFL-run={row['p23_vfl_run']}; indel={row['p24_indel'] or '-'}]"
                for row in rows
            ]
        )
        + "\n"
    )

    print(f"outdir={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
