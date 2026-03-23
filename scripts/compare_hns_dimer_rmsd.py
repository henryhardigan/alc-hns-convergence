from __future__ import annotations

import csv
import math
from pathlib import Path

from Bio.PDB import MMCIFParser, Superimposer


ROOT = Path(__file__).resolve().parents[1]
HOMO_CIF = ROOT / "models" / "hns_hns" / "model_0.cif"
HETERO_CIF = ROOT / "models" / "alc_hns" / "model_0.cif"

# Preferred manual convergence register from the manuscript work.
ALC_ALIGN = "EYRNRRLEQARDML-PDAVEEMKVFLENQLA"
HNS_ALIGN = "EERTRKLQQYREMLIADGIDPNEL-L-NSLA"

OUTDIR = ROOT / "analysis"
SUMMARY_TSV = OUTDIR / "hns_homodimer_vs_alc_heterodimer_rmsd.tsv"
SUMMARY_MD = OUTDIR / "hns_homodimer_vs_alc_heterodimer_rmsd.md"


def aligned_position_pairs() -> list[tuple[int, int, str, str]]:
    pairs: list[tuple[int, int, str, str]] = []
    ai = hi = 0
    for a, h in zip(ALC_ALIGN, HNS_ALIGN):
        if a != "-":
            ai += 1
        if h != "-":
            hi += 1
        if a != "-" and h != "-":
            pairs.append((ai, hi, a, h))
    return pairs


def ca_residues(chain):
    return [res for res in chain if res.id[0] == " "]


def ca_atoms(residues, one_based_indices: list[int]):
    return [residues[i - 1]["CA"] for i in one_based_indices]


def transformed_partner_rmsd(
    hns_receptor,
    hns_partner,
    hetero_receptor,
    hetero_partner,
    pairs,
):
    sup = Superimposer()
    sup.set_atoms(
        [res["CA"] for res in hns_receptor],
        [res["CA"] for res in hetero_receptor],
    )
    rot, tran = sup.rotran
    ssq = []
    for alc_i, hns_i, _, _ in pairs:
        ref = hns_partner[hns_i - 1]["CA"].coord
        mov = hetero_partner[alc_i - 1]["CA"].coord @ rot + tran
        ssq.append(((mov - ref) ** 2).sum())
    return sup.rms, math.sqrt(sum(ssq) / len(ssq))


def partner_shape_rmsd(hns_partner, hetero_partner, pairs):
    sup = Superimposer()
    sup.set_atoms(
        [hns_partner[hns_i - 1]["CA"] for _, hns_i, _, _ in pairs],
        [hetero_partner[alc_i - 1]["CA"] for alc_i, _, _, _ in pairs],
    )
    return sup.rms


def combined_state_rmsd(hns_receptor, hns_partner, hetero_receptor, hetero_partner, pairs):
    sup = Superimposer()
    sup.set_atoms(
        [res["CA"] for res in hns_receptor],
        [res["CA"] for res in hetero_receptor],
    )
    rot, tran = sup.rotran
    ssq = []
    for i in range(len(hns_receptor)):
        ref = hns_receptor[i]["CA"].coord
        mov = hetero_receptor[i]["CA"].coord @ rot + tran
        ssq.append(((mov - ref) ** 2).sum())
    for alc_i, hns_i, _, _ in pairs:
        ref = hns_partner[hns_i - 1]["CA"].coord
        mov = hetero_partner[alc_i - 1]["CA"].coord @ rot + tran
        ssq.append(((mov - ref) ** 2).sum())
    return math.sqrt(sum(ssq) / len(ssq))


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    parser = MMCIFParser(QUIET=True)
    homo = parser.get_structure("homo", str(HOMO_CIF))[0]
    hetero = parser.get_structure("hetero", str(HETERO_CIF))[0]

    h_a = ca_residues(homo["A"])
    h_b = ca_residues(homo["B"])
    t_a = ca_residues(hetero["A"])
    t_b = ca_residues(hetero["B"])
    pairs = aligned_position_pairs()

    receptor_a_rmsd, partner_b_rmsd = transformed_partner_rmsd(h_a, h_b, t_a, t_b, pairs)
    receptor_b_rmsd, partner_a_rmsd = transformed_partner_rmsd(h_b, h_a, t_a, t_b, pairs)
    partner_shape = partner_shape_rmsd(h_b, t_b, pairs)
    combined = combined_state_rmsd(h_a, h_b, t_a, t_b, pairs)

    rows = [
        {
            "comparison": "hetero_A_to_homo_A__partner_to_homo_B",
            "receptor_ca_rmsd_A": f"{receptor_a_rmsd:.3f}",
            "partner_ca_rmsd_A": f"{partner_b_rmsd:.3f}",
            "aligned_partner_positions": str(len(pairs)),
        },
        {
            "comparison": "hetero_A_to_homo_B__partner_to_homo_A",
            "receptor_ca_rmsd_A": f"{receptor_b_rmsd:.3f}",
            "partner_ca_rmsd_A": f"{partner_a_rmsd:.3f}",
            "aligned_partner_positions": str(len(pairs)),
        },
        {
            "comparison": "partner_shape_best_fit_only",
            "receptor_ca_rmsd_A": "",
            "partner_ca_rmsd_A": f"{partner_shape:.3f}",
            "aligned_partner_positions": str(len(pairs)),
        },
        {
            "comparison": "combined_state_rmsd_after_receptor_superposition",
            "receptor_ca_rmsd_A": "",
            "partner_ca_rmsd_A": f"{combined:.3f}",
            "aligned_partner_positions": str(29 + len(pairs)),
        },
    ]

    with SUMMARY_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "comparison",
                "receptor_ca_rmsd_A",
                "partner_ca_rmsd_A",
                "aligned_partner_positions",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    SUMMARY_MD.write_text(
        "\n".join(
            [
                "# H-NS/Alc Dimer-State RMSD",
                "",
                f"- Homodimer model: `{HOMO_CIF}`",
                f"- Heterodimer model: `{HETERO_CIF}`",
                f"- Preferred Alc/H-NS register positions compared on partner chain: `{len(pairs)}`",
                "",
                "Key results:",
                f"- H-NS receptor superposition RMSD: `{receptor_a_rmsd:.3f} A`",
                f"- Partner-chain RMSD after receptor superposition: `{partner_b_rmsd:.3f} A`",
                f"- Best-fit partner monomer RMSD over aligned positions: `{partner_shape:.3f} A`",
                f"- Combined dimer-state RMSD (29 receptor + {len(pairs)} partner positions): `{combined:.3f} A`",
                "",
                "Interpretation:",
                "- The H-NS chain itself overlays reasonably after superposition, but the partner chain does not reproduce the H-NS:H-NS partner state tightly.",
                "- This supports a related interface geometry rather than a near-identical dimer state.",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
