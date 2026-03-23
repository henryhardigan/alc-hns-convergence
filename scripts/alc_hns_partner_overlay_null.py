from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, Superimposer


HOMO_CIF = Path("/Users/henryhardigan/Downloads/H-NS:H-NS/fold_2026_03_17_09_41_model_0.cif")
HETERO_CIF = Path("/Users/henryhardigan/Downloads/alc:H-NS/fold_2026_03_17_09_43_model_0.cif")
OUTDIR = Path(
    "/Users/henryhardigan/Desktop/src/linear-complement-alignment/results/alc_hns_convergence"
)
TSV_PATH = OUTDIR / "alc_hns_partner_overlay_null.tsv"
SUMMARY_MD = OUTDIR / "alc_hns_partner_overlay_null.md"

ALC_SEQ = "EYRNRRLEQARDMLPDAVEEMKVFLENQLA"
HNS_SEQ = "EERTRKLQQYREMLIADGIDPNELLNSLA"

# Preferred manual register used throughout the manuscript.
OBS_ALC_ALIGN = "EYRNRRLEQARDML-PDAVEEMKVFLENQLA"
OBS_HNS_ALIGN = "EERTRKLQQYREMLIADGIDPNEL-L-NSLA"


def pairs_from_alignment(alc_aln: str, hns_aln: str) -> list[tuple[int, int, str, str]]:
    pairs: list[tuple[int, int, str, str]] = []
    ai = hi = 0
    for a, h in zip(alc_aln, hns_aln):
        if a != "-":
            ai += 1
        if h != "-":
            hi += 1
        if a != "-" and h != "-":
            pairs.append((ai, hi, a, h))
    return pairs


def pairs_from_drop_positions(
    alc_seq: str,
    hns_seq: str,
    alc_drops: tuple[int, int],
    hns_drop: int,
) -> list[tuple[int, int, str, str]]:
    alc_keep = [i for i in range(1, len(alc_seq) + 1) if i not in alc_drops]
    hns_keep = [i for i in range(1, len(hns_seq) + 1) if i != hns_drop]
    assert len(alc_keep) == len(hns_keep) == 28
    return [(ai, hi, alc_seq[ai - 1], hns_seq[hi - 1]) for ai, hi in zip(alc_keep, hns_keep)]


def ca_residues(chain):
    return [res for res in chain if res.id[0] == " "]


def score_pairs(matrix, pairs: list[tuple[int, int, str, str]]) -> float:
    return float(sum(matrix[a, h] for _, _, a, h in pairs))


def rmsd_for_pairs(hns_partner, alc_partner, pairs: list[tuple[int, int, str, str]]) -> float:
    sup = Superimposer()
    fixed = [hns_partner[h_i - 1]["CA"] for _, h_i, _, _ in pairs]
    moving = [alc_partner[a_i - 1]["CA"] for a_i, _, _, _ in pairs]
    sup.set_atoms(fixed, moving)
    return float(sup.rms)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    blosum62 = substitution_matrices.load("BLOSUM62")
    parser = MMCIFParser(QUIET=True)
    homo = parser.get_structure("homo", str(HOMO_CIF))[0]
    hetero = parser.get_structure("hetero", str(HETERO_CIF))[0]
    hns_partner = ca_residues(homo["B"])
    alc_partner = ca_residues(hetero["B"])

    observed_pairs = pairs_from_alignment(OBS_ALC_ALIGN, OBS_HNS_ALIGN)
    observed_score = score_pairs(blosum62, observed_pairs)
    observed_rmsd = rmsd_for_pairs(hns_partner, alc_partner, observed_pairs)

    rows = []
    for alc_drops in itertools.combinations(range(1, len(ALC_SEQ) + 1), 2):
        for hns_drop in range(1, len(HNS_SEQ) + 1):
            pairs = pairs_from_drop_positions(ALC_SEQ, HNS_SEQ, alc_drops, hns_drop)
            score = score_pairs(blosum62, pairs)
            rmsd = rmsd_for_pairs(hns_partner, alc_partner, pairs)
            rows.append(
                {
                    "alc_drop_positions": f"{alc_drops[0]},{alc_drops[1]}",
                    "hns_drop_position": str(hns_drop),
                    "blosum62_score": score,
                    "partner_only_rmsd_A": rmsd,
                }
            )

    # sort for output readability
    rows.sort(key=lambda r: (r["partner_only_rmsd_A"], -r["blosum62_score"]))

    with TSV_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "alc_drop_positions",
                "hns_drop_position",
                "blosum62_score",
                "partner_only_rmsd_A",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **row,
                    "blosum62_score": f"{row['blosum62_score']:.1f}",
                    "partner_only_rmsd_A": f"{row['partner_only_rmsd_A']:.3f}",
                }
            )

    all_rmsds = [r["partner_only_rmsd_A"] for r in rows]
    all_scores = [r["blosum62_score"] for r in rows]
    better_or_equal_rmsd = sum(r <= observed_rmsd + 1e-9 for r in all_rmsds)
    better_or_equal_score = sum(s >= observed_score - 1e-9 for s in all_scores)

    # chemistry-favored subset: BLOSUM score at least as good as the observed register
    score_filtered = [r["partner_only_rmsd_A"] for r in rows if r["blosum62_score"] >= observed_score - 1e-9]
    # near-optimal chemistry subset: within 2 BLOSUM points of observed
    near_filtered = [r["partner_only_rmsd_A"] for r in rows if r["blosum62_score"] >= observed_score - 2.0]

    def empirical_leq(values: list[float], threshold: float) -> tuple[int, int]:
        count = sum(v <= threshold + 1e-9 for v in values)
        return count, len(values)

    obs_rank_all = empirical_leq(all_rmsds, observed_rmsd)
    obs_rank_score = empirical_leq(score_filtered, observed_rmsd)
    obs_rank_near = empirical_leq(near_filtered, observed_rmsd)

    best_rmsd_row = min(rows, key=lambda r: r["partner_only_rmsd_A"])
    best_score_row = max(rows, key=lambda r: r["blosum62_score"])

    SUMMARY_MD.write_text(
        "\n".join(
            [
                "# Alc/H-NS Partner-Overlay Null",
                "",
                "Null design:",
                "- Enumerated every order-preserving 28-position alignment between the 30-residue Alc partner segment and the 29-residue H-NS partner segment.",
                "- This is equivalent to dropping 2 Alc positions and 1 H-NS position, then pairing the remaining residues in order.",
                "- Total alignments considered: `12615`.",
                "- Each alignment was scored by BLOSUM62 and evaluated by best-fit CA RMSD of the partner chains only.",
                "",
                "Observed preferred register:",
                f"- BLOSUM62 score: `{observed_score:.1f}`",
                f"- Partner-only RMSD: `{observed_rmsd:.3f} A`",
                "",
                "Extremes over the full alignment space:",
                f"- Best (lowest) partner-only RMSD among all 12,615 alignments: `{best_rmsd_row['partner_only_rmsd_A']:.3f} A` (drops Alc {best_rmsd_row['alc_drop_positions']}; H-NS {best_rmsd_row['hns_drop_position']})",
                f"- Best BLOSUM62 score among all 12,615 alignments: `{best_score_row['blosum62_score']:.1f}` (partner-only RMSD `{best_score_row['partner_only_rmsd_A']:.3f} A`)",
                "",
                "Observed-register rank by partner-only RMSD:",
                f"- Against all alignments: `{obs_rank_all[0]}/{obs_rank_all[1]}` alignments have RMSD <= observed",
                f"- Among alignments with BLOSUM62 score >= observed: `{obs_rank_score[0]}/{obs_rank_score[1]}`",
                f"- Among alignments within 2 BLOSUM62 points of observed: `{obs_rank_near[0]}/{obs_rank_near[1]}`",
                "",
                "Interpretation:",
                "- This null asks whether the preferred convergence register produces an unusually good structural overlay relative to the full space of order-preserving 28-position mappings.",
                "- The chemistry-filtered subsets are the most relevant comparison for a convergent-mimicry claim, because they hold sequence plausibility roughly constant while testing structural specificity.",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
