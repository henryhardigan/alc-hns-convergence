#!/usr/bin/env python3
"""Additional nulls supporting the Alc/H-NS host-clade convergence claim."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import random
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
CONSTRAINED_NULL_SCRIPT = SCRIPT_DIR / "alc_hns_constrained_null.py"


def load_constrained_module():
    spec = importlib.util.spec_from_file_location("alc_hns_constrained_null_mod", CONSTRAINED_NULL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {CONSTRAINED_NULL_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


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
        help="Output directory.",
    )
    parser.add_argument("--shuffle-replicates", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


def host_branch(entry: str) -> bool:
    return entry in {"HNS_ECO57", "HNS_ECOL6", "HNS_ECOLI", "HNS_SHIFL"}


def exact_family_branch_pvalue(family_size: int, host_branch_size: int, hit_count: int) -> float:
    # Probability that all hit_count hits land in the host branch under uniform placement.
    if hit_count > host_branch_size:
        return 0.0
    return math.comb(host_branch_size, hit_count) * math.comb(
        family_size - host_branch_size, 0
    ) / math.comb(family_size, hit_count)


def main() -> int:
    args = parse_args()
    mod = load_constrained_module()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed)
    records = mod.load_reviewed_records(args.fasta_glob)
    matrix = mod.substitution_matrices.load("BLOSUM62")

    by_accession = {record.accession: record for record in records}
    alc_record = by_accession[mod.ALC_ACCESSION]
    alc_span = mod.extract_window(alc_record, mod.ALC_SPAN_START_1BASED, mod.ALC_SPAN_LEN)
    if len(alc_span) != mod.ALC_SPAN_LEN:
        raise SystemExit("Could not extract Alc span.")

    family_rows = []
    for record in records:
        if not mod.is_hns_family(record):
            continue
        hns_span = mod.extract_window(record, mod.HNS_SPAN_START_1BASED, mod.HNS_SPAN_LEN)
        if len(hns_span) != mod.HNS_SPAN_LEN:
            continue
        score, aligned_alc, aligned_hns = mod.constrained_score(
            alc_span,
            hns_span,
            matrix,
            mod.DEFAULT_GAP_OPEN,
            mod.DEFAULT_GAP_EXTEND,
        )
        family_rows.append(
            {
                "entry": record.entry,
                "accession": record.accession,
                "organism": record.organism,
                "gene": record.gene,
                "hns_span_52_80": hns_span,
                "observed_score": score,
                "host_branch": "yes" if host_branch(record.entry) else "no",
            }
        )

    observed_best = max(float(row["observed_score"]) for row in family_rows)
    observed_top = [row for row in family_rows if float(row["observed_score"]) == observed_best]
    observed_top_host = [row for row in observed_top if row["host_branch"] == "yes"]
    observed_nonhost_best = max(float(row["observed_score"]) for row in family_rows if row["host_branch"] == "no")
    observed_margin = observed_best - observed_nonhost_best

    shuffle_rows = []
    any_host_top_count = 0
    all_host_top_count = 0
    margin_ge_observed_count = 0
    any_host_top_and_margin_ge_observed = 0
    top_score_ge_observed = 0

    family_spans = [(row["entry"], row["host_branch"], row["hns_span_52_80"]) for row in family_rows]

    for rep in range(1, args.shuffle_replicates + 1):
        shuffled = list(alc_span)
        rng.shuffle(shuffled)
        shuffled_seq = "".join(shuffled)
        scores = []
        for entry, is_host, hns_span in family_spans:
            score, _, _ = mod.constrained_score(
                shuffled_seq,
                hns_span,
                matrix,
                mod.DEFAULT_GAP_OPEN,
                mod.DEFAULT_GAP_EXTEND,
            )
            scores.append((entry, is_host, score))
        best = max(score for _, _, score in scores)
        top_entries = [(entry, is_host) for entry, is_host, score in scores if score == best]
        any_host = any(is_host == "yes" for _, is_host in top_entries)
        all_host = all(is_host == "yes" for _, is_host in top_entries)
        nonhost_best = max(score for _, is_host, score in scores if is_host == "no")
        margin = best - nonhost_best

        if any_host:
            any_host_top_count += 1
        if all_host:
            all_host_top_count += 1
        if margin >= observed_margin - 1e-9:
            margin_ge_observed_count += 1
        if any_host and margin >= observed_margin - 1e-9:
            any_host_top_and_margin_ge_observed += 1
        if best >= observed_best - 1e-9:
            top_score_ge_observed += 1

        shuffle_rows.append(
            {
                "replicate": rep,
                "shuffled_alc_span": shuffled_seq,
                "best_score": f"{best:.1f}",
                "top_entries": ",".join(entry for entry, _ in top_entries),
                "any_host_branch_top": "yes" if any_host else "no",
                "all_host_branch_top": "yes" if all_host else "no",
                "margin_vs_best_nonhost": f"{margin:.1f}",
            }
        )

    with (outdir / "claim_nulls_shuffle_branch.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "replicate",
                "shuffled_alc_span",
                "best_score",
                "top_entries",
                "any_host_branch_top",
                "all_host_branch_top",
                "margin_vs_best_nonhost",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(shuffle_rows)

    exact_hit_count = 4
    family_size = len(family_rows)
    host_branch_size = sum(row["host_branch"] == "yes" for row in family_rows)
    exact_p = exact_family_branch_pvalue(family_size, host_branch_size, exact_hit_count)

    summary_lines = [
        "# Alc/H-NS Claim-Completion Nulls",
        "",
        "## Exact family-branch null",
        "",
        f"- Reviewed H-NS/StpA family proteins: `{family_size}`",
        f"- Host branch proteins (E. coli/Shigella H-NS only): `{host_branch_size}`",
        f"- Exact primary-motif hits within family: `{exact_hit_count}`",
        f"- Exact hits observed in host branch: `{exact_hit_count}`",
        f"- Uniform-placement branch p-value: `{exact_p:.6f}`",
        "",
        "## Composition-preserving branch null",
        "",
        f"- Shuffles of the Alc 30-mer: `{args.shuffle_replicates}`",
        f"- Observed best constrained family score: `{observed_best:.1f}`",
        f"- Observed top entries: `{','.join(row['entry'] for row in observed_top)}`",
        f"- Observed best non-host score: `{observed_nonhost_best:.1f}`",
        f"- Observed host-vs-nonhost margin: `{observed_margin:.1f}`",
        "",
        "Empirical frequencies over shuffled Alc sequences:",
        f"- Any host-branch top hit: `{any_host_top_count}/{args.shuffle_replicates}`",
        f"- All top hits confined to host branch: `{all_host_top_count}/{args.shuffle_replicates}`",
        f"- Margin >= observed (`{observed_margin:.1f}`): `{margin_ge_observed_count}/{args.shuffle_replicates}`",
        f"- Any host-branch top hit AND margin >= observed: `{any_host_top_and_margin_ge_observed}/{args.shuffle_replicates}`",
        f"- Best score >= observed (`{observed_best:.1f}`): `{top_score_ge_observed}/{args.shuffle_replicates}`",
        "",
        "Interpretation:",
        "- The exact family-branch null asks whether all exact family hits would fall into the E. coli/Shigella H-NS branch by chance if family placement were uniform.",
        "- The shuffled-Alc branch null asks whether Alc-like composition alone tends to pick the host branch as the best constrained family match, and whether it does so with the observed score margin over non-host H-NS/StpA proteins.",
    ]
    (outdir / "claim_nulls_summary.md").write_text("\n".join(summary_lines) + "\n")

    with (outdir / "claim_nulls_metrics.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t")
        writer.writeheader()
        metrics = [
            ("family_size", family_size),
            ("host_branch_size", host_branch_size),
            ("exact_family_hit_count", exact_hit_count),
            ("exact_family_branch_pvalue", f"{exact_p:.6f}"),
            ("shuffle_replicates", args.shuffle_replicates),
            ("observed_best_score", f"{observed_best:.1f}"),
            ("observed_nonhost_best_score", f"{observed_nonhost_best:.1f}"),
            ("observed_margin", f"{observed_margin:.1f}"),
            ("any_host_branch_top_count", any_host_top_count),
            ("all_host_branch_top_count", all_host_top_count),
            ("margin_ge_observed_count", margin_ge_observed_count),
            ("any_host_branch_top_and_margin_ge_observed", any_host_top_and_margin_ge_observed),
            ("top_score_ge_observed_count", top_score_ge_observed),
        ]
        for metric, value in metrics:
            writer.writerow({"metric": metric, "value": value})

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
