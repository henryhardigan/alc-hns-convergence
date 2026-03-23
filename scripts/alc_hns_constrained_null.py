#!/usr/bin/env python3
"""Composition-preserving null for the Alc/H-NS quadripartite constrained score."""

from __future__ import annotations

import argparse
import csv
import random
import re
from dataclasses import dataclass
from pathlib import Path

from Bio.Align import substitution_matrices


HEADER_RE = re.compile(r"^(?P<db>sp|tr)\|(?P<acc>[^|]+)\|(?P<entry>\S+)")
ORGANISM_RE = re.compile(r"\bOS=(.+?) OX=")
TAXON_RE = re.compile(r"\bOX=(\d+)\b")
GN_RE = re.compile(r"\bGN=([^ ]+)")

ALC_ACCESSION = "P04546"
ALC_SPAN_START_1BASED = 53
ALC_SPAN_LEN = 30
HNS_SPAN_START_1BASED = 52
HNS_SPAN_LEN = 29

DEFAULT_GAP_OPEN = -5.0
DEFAULT_GAP_EXTEND = -0.5

ALC_PARTS = (
    ("block1", 14),
    ("bridge1", 1),
    ("block2", 4),
    ("bridge2", 4),
    ("block3", 2),
    ("bridge3", 1),
    ("block4", 4),
)
HNS_PARTS = (
    ("block1", 14),
    ("bridge1", 2),
    ("block2", 4),
    ("bridge2", 3),
    ("block3", 2),
    ("bridge3", 0),
    ("block4", 4),
)


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
    parser.add_argument(
        "--shuffle-replicates",
        type=int,
        default=20000,
        help="Number of composition-preserving shuffles of the Alc quadripartite span.",
    )
    parser.add_argument("--seed", type=int, default=7, help="Random seed.")
    parser.add_argument("--gap-open", type=float, default=DEFAULT_GAP_OPEN)
    parser.add_argument("--gap-extend", type=float, default=DEFAULT_GAP_EXTEND)
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


def is_hns_family(record: Record) -> bool:
    return (
        "DNA-binding protein H-NS" in record.header
        or "DNA-binding protein StpA" in record.header
        or "VicH" in record.header
    )


def extract_window(record: Record, start_1based: int, length: int) -> str:
    start0 = start_1based - 1
    end0 = start0 + length
    if len(record.sequence) < end0:
        return ""
    return record.sequence[start0:end0]


def split_by_parts(seq: str, parts: tuple[tuple[str, int], ...]) -> dict[str, str]:
    out: dict[str, str] = {}
    cursor = 0
    for label, length in parts:
        out[label] = seq[cursor : cursor + length]
        cursor += length
    if cursor != len(seq):
        raise ValueError(f"Part lengths do not consume sequence: {cursor} != {len(seq)}")
    return out


def pair_score(matrix, seq_a: str, seq_b: str) -> float:
    return float(sum(matrix[a, b] for a, b in zip(seq_a, seq_b)))


def gap_cost(length: int, gap_open: float, gap_extend: float) -> float:
    if length <= 0:
        return 0.0
    return gap_open + gap_extend * (length - 1)


def global_affine_bridge(seq_a: str, seq_b: str, matrix, gap_open: float, gap_extend: float) -> tuple[float, str, str]:
    n = len(seq_a)
    m = len(seq_b)
    neg_inf = float("-inf")
    M = [[neg_inf] * (m + 1) for _ in range(n + 1)]
    X = [[neg_inf] * (m + 1) for _ in range(n + 1)]  # gap in seq_b
    Y = [[neg_inf] * (m + 1) for _ in range(n + 1)]  # gap in seq_a
    tb_M = [[None] * (m + 1) for _ in range(n + 1)]
    tb_X = [[None] * (m + 1) for _ in range(n + 1)]
    tb_Y = [[None] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0.0
    for i in range(1, n + 1):
        X[i][0] = gap_open if i == 1 else X[i - 1][0] + gap_extend
        tb_X[i][0] = ("M" if i == 1 else "X", i - 1, 0)
    for j in range(1, m + 1):
        Y[0][j] = gap_open if j == 1 else Y[0][j - 1] + gap_extend
        tb_Y[0][j] = ("M" if j == 1 else "Y", 0, j - 1)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = float(matrix[seq_a[i - 1], seq_b[j - 1]])
            prevs = [
                (M[i - 1][j - 1], "M"),
                (X[i - 1][j - 1], "X"),
                (Y[i - 1][j - 1], "Y"),
            ]
            prev_score, prev_state = max(prevs, key=lambda item: item[0])
            M[i][j] = prev_score + score
            tb_M[i][j] = (prev_state, i - 1, j - 1)

            prevs_x = [
                (M[i - 1][j] + gap_open, "M"),
                (X[i - 1][j] + gap_extend, "X"),
                (Y[i - 1][j] + gap_open, "Y"),
            ]
            prev_score_x, prev_state_x = max(prevs_x, key=lambda item: item[0])
            X[i][j] = prev_score_x
            tb_X[i][j] = (prev_state_x, i - 1, j)

            prevs_y = [
                (M[i][j - 1] + gap_open, "M"),
                (Y[i][j - 1] + gap_extend, "Y"),
                (X[i][j - 1] + gap_open, "X"),
            ]
            prev_score_y, prev_state_y = max(prevs_y, key=lambda item: item[0])
            Y[i][j] = prev_score_y
            tb_Y[i][j] = (prev_state_y, i, j - 1)

    end_state, best_score = max(
        (("M", M[n][m]), ("X", X[n][m]), ("Y", Y[n][m])),
        key=lambda item: item[1],
    )
    aligned_a: list[str] = []
    aligned_b: list[str] = []
    i, j, state = n, m, end_state
    while i > 0 or j > 0:
        if state == "M":
            prev_state, pi, pj = tb_M[i][j]
            aligned_a.append(seq_a[i - 1])
            aligned_b.append(seq_b[j - 1])
        elif state == "X":
            prev_state, pi, pj = tb_X[i][j]
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
        else:
            prev_state, pi, pj = tb_Y[i][j]
            aligned_a.append("-")
            aligned_b.append(seq_b[j - 1])
        i, j, state = pi, pj, prev_state

    return best_score, "".join(reversed(aligned_a)), "".join(reversed(aligned_b))


def best_global_bridge(seq_a: str, seq_b: str, matrix, gap_open: float, gap_extend: float) -> tuple[float, str, str]:
    if not seq_a and not seq_b:
        return 0.0, "", ""
    if not seq_a:
        return gap_cost(len(seq_b), gap_open, gap_extend), "-" * len(seq_b), seq_b
    if not seq_b:
        return gap_cost(len(seq_a), gap_open, gap_extend), seq_a, "-" * len(seq_a)
    return global_affine_bridge(seq_a, seq_b, matrix, gap_open, gap_extend)


def constrained_score(
    query_seq: str,
    target_seq: str,
    matrix,
    gap_open: float,
    gap_extend: float,
) -> tuple[float, str, str]:
    q_parts = split_by_parts(query_seq, ALC_PARTS)
    t_parts = split_by_parts(target_seq, HNS_PARTS)
    total = 0.0
    aligned_q: list[str] = []
    aligned_t: list[str] = []
    for label, _len in ALC_PARTS:
        q = q_parts[label]
        t = t_parts[label]
        if label.startswith("block"):
            score = pair_score(matrix, q, t)
            aq = q
            at = t
        else:
            score, aq, at = best_global_bridge(q, t, matrix, gap_open, gap_extend)
        total += score
        aligned_q.append(aq)
        aligned_t.append(at)
    return total, "".join(aligned_q), "".join(aligned_t)


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    records = load_reviewed_records(args.fasta_glob)
    matrix = substitution_matrices.load("BLOSUM62")
    rng = random.Random(args.seed)

    by_accession = {record.accession: record for record in records}
    alc = by_accession[ALC_ACCESSION]
    alc_span = extract_window(alc, ALC_SPAN_START_1BASED, ALC_SPAN_LEN)
    if len(alc_span) != ALC_SPAN_LEN:
        raise SystemExit("Could not extract Alc quadripartite span.")

    family_rows: list[dict[str, object]] = []
    for record in records:
        if not is_hns_family(record):
            continue
        hns_span = extract_window(record, HNS_SPAN_START_1BASED, HNS_SPAN_LEN)
        if len(hns_span) != HNS_SPAN_LEN:
            continue
        score, aligned_q, aligned_t = constrained_score(
            alc_span,
            hns_span,
            matrix,
            args.gap_open,
            args.gap_extend,
        )
        family_rows.append(
            {
                "accession": record.accession,
                "entry": record.entry,
                "organism": record.organism,
                "gene": record.gene,
                "hns_span_52_80": hns_span,
                "constrained_score_vs_alc_53_82": f"{score:.1f}",
                "aligned_alc": aligned_q,
                "aligned_hns": aligned_t,
            }
        )

    family_rows.sort(
        key=lambda row: (-float(row["constrained_score_vs_alc_53_82"]), str(row["entry"]))
    )
    write_tsv(
        outdir / "constrained_null_family_scores.tsv",
        family_rows,
        [
            "accession",
            "entry",
            "organism",
            "gene",
            "hns_span_52_80",
            "constrained_score_vs_alc_53_82",
            "aligned_alc",
            "aligned_hns",
        ],
    )

    observed_best = float(family_rows[0]["constrained_score_vs_alc_53_82"])
    best_family = family_rows[0]["entry"]

    null_rows: list[dict[str, object]] = []
    exceed_count = 0
    for replicate in range(1, args.shuffle_replicates + 1):
        shuffled = list(alc_span)
        rng.shuffle(shuffled)
        shuffled_seq = "".join(shuffled)
        best_score = None
        best_entry = None
        for row in family_rows:
            score, _aq, _at = constrained_score(
                shuffled_seq,
                row["hns_span_52_80"],
                matrix,
                args.gap_open,
                args.gap_extend,
            )
            if best_score is None or score > best_score:
                best_score = score
                best_entry = row["entry"]
        assert best_score is not None and best_entry is not None
        if best_score >= observed_best:
            exceed_count += 1
        null_rows.append(
            {
                "replicate": replicate,
                "shuffled_alc_span": shuffled_seq,
                "best_constrained_score_vs_hns_family": f"{best_score:.1f}",
                "best_entry": best_entry,
            }
        )

    write_tsv(
        outdir / "constrained_null_scores.tsv",
        null_rows,
        ["replicate", "shuffled_alc_span", "best_constrained_score_vs_hns_family", "best_entry"],
    )

    empirical_p = (exceed_count + 1) / (args.shuffle_replicates + 1)
    summary = outdir / "constrained_null_summary.md"
    summary.write_text(
        "\n".join(
            [
                "# Alc/H-NS Constrained Null Summary",
                "",
                "## Model",
                "",
                "Composition-preserving shuffles of the Alc quadripartite span were scored against the reviewed H-NS/StpA family using the quadripartite-constrained alignment model.",
                "",
                "Quadripartite motif:",
                "",
                "`E-X-R-[NT]-R-[RK]-L-[QE]-Q-X-R-[DE]-M-L ... D-[AG]-[IV]-[DE] ... [FL]-L ... N-[QS]-L-A`",
                "",
                "## Parameters",
                "",
                f"- gap open: {args.gap_open}",
                f"- gap extend: {args.gap_extend}",
                f"- shuffles: {args.shuffle_replicates}",
                f"- seed: {args.seed}",
                "",
                "## Observed family score",
                "",
                f"- observed best constrained score: {observed_best:.1f}",
                f"- best family entry: {best_family}",
                "",
                "## Null result",
                "",
                f"- empirical P(best constrained score >= observed): {empirical_p:.6f}",
                f"- exceedances: {exceed_count}/{args.shuffle_replicates}",
                "",
                "## Top family scores",
                "",
            ]
            + [
                f"- {row['entry']} ({row['accession']}): score={row['constrained_score_vs_alc_53_82']} `{row['hns_span_52_80']}`"
                for row in family_rows[:8]
            ]
        )
        + "\n"
    )

    print(f"outdir={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
