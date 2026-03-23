from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser, Superimposer


ROOT = Path(__file__).resolve().parents[1]
HOMO_CIF = ROOT / "models" / "hns_hns" / "model_0.cif"
HETERO_CIF = ROOT / "models" / "alc_hns" / "model_0.cif"
OUTDIR = ROOT / "figures"
OUTPNG = OUTDIR / "figure3_structure_overlay.png"

ALC_ALIGN = "EYRNRRLEQARDML-PDAVEEMKVFLENQLA"
HNS_ALIGN = "EERTRKLQQYREMLIADGIDPNEL-L-NSLA"


def aligned_position_pairs() -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    ai = hi = 0
    for a, h in zip(ALC_ALIGN, HNS_ALIGN):
        if a != "-":
            ai += 1
        if h != "-":
            hi += 1
        if a != "-" and h != "-":
            pairs.append((ai, hi))
    return pairs


def residues(chain):
    return [res for res in chain if res.id[0] == " "]


def ca_coords(reslist):
    return np.array([res["CA"].coord for res in reslist], dtype=float)


def transform(coords: np.ndarray, rot: np.ndarray, tran: np.ndarray) -> np.ndarray:
    return coords @ rot + tran


def pca_project(*coord_sets: np.ndarray) -> list[np.ndarray]:
    stacked = np.vstack(coord_sets)
    centered = stacked - stacked.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    basis = vt[:2].T
    projected = []
    start = 0
    for arr in coord_sets:
        n = len(arr)
        projected.append((stacked[start : start + n] - stacked.mean(axis=0)) @ basis)
        start += n
    return projected


def set_equal_2d(ax, arrs: list[np.ndarray], pad: float = 1.0) -> None:
    pts = np.vstack(arrs)
    xmin, ymin = pts.min(axis=0)
    xmax, ymax = pts.max(axis=0)
    cx = (xmin + xmax) / 2
    cy = (ymin + ymax) / 2
    half = max(xmax - xmin, ymax - ymin) / 2 + pad
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)


def plot_trace(ax, xy: np.ndarray, color: str, label: str, lw: float = 3.0, alpha: float = 1.0):
    ax.plot(xy[:, 0], xy[:, 1], "-", color=color, lw=lw, alpha=alpha, solid_capstyle="round", label=label)
    ax.scatter(xy[:, 0], xy[:, 1], s=18, color=color, alpha=alpha, zorder=3)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    parser = MMCIFParser(QUIET=True)
    homo = parser.get_structure("homo", str(HOMO_CIF))[0]
    hetero = parser.get_structure("hetero", str(HETERO_CIF))[0]

    h_a = residues(homo["A"])
    h_b = residues(homo["B"])
    t_a = residues(hetero["A"])
    t_b = residues(hetero["B"])

    sup_receptor = Superimposer()
    sup_receptor.set_atoms([res["CA"] for res in h_a], [res["CA"] for res in t_a])
    rot_r, tran_r = sup_receptor.rotran
    h_a_xy, h_b_xy, t_b_on_receptor_xy = pca_project(
        ca_coords(h_a),
        ca_coords(h_b),
        transform(ca_coords(t_b), rot_r, tran_r),
    )

    pairs = aligned_position_pairs()
    h_partner = np.array([h_b[hi - 1]["CA"].coord for _, hi in pairs], dtype=float)
    t_partner = np.array([t_b[ai - 1]["CA"].coord for ai, _ in pairs], dtype=float)
    sup_partner = Superimposer()
    sup_partner.set_atoms([h_b[hi - 1]["CA"] for _, hi in pairs], [t_b[ai - 1]["CA"] for ai, _ in pairs])
    rot_p, tran_p = sup_partner.rotran
    h_partner_xy, t_partner_xy = pca_project(
        h_partner,
        transform(t_partner, rot_p, tran_p),
    )

    plt.rcParams.update(
        {
            "font.size": 11,
            "font.family": "DejaVu Sans",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.03, h_pad=0.08, hspace=0.08, wspace=0.04)

    ax = axes[0]
    plot_trace(ax, h_a_xy, "#333333", "H-NS receptor", lw=3.4)
    plot_trace(ax, h_b_xy, "#3B82F6", "H-NS partner", lw=3.0)
    plot_trace(ax, t_b_on_receptor_xy, "#D946EF", "Alc partner", lw=3.0, alpha=0.95)
    set_equal_2d(ax, [h_a_xy, h_b_xy, t_b_on_receptor_xy])
    ax.set_title("Receptor-Aligned Overlay", weight="bold", pad=18)
    ax.text(
        0.02,
        0.02,
        f"Receptor Cα RMSD = {sup_receptor.rms:.2f} Å\nPartner Cα RMSD = 4.99 Å",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
    )
    ax.axis("off")

    ax = axes[1]
    plot_trace(ax, h_partner_xy, "#3B82F6", "H-NS partner", lw=3.2)
    plot_trace(ax, t_partner_xy, "#D946EF", "Alc partner", lw=3.2, alpha=0.95)
    set_equal_2d(ax, [h_partner_xy, t_partner_xy])
    ax.set_title("Partner-Only Best Fit", weight="bold", pad=18)
    ax.text(
        0.02,
        0.02,
        f"28 aligned Cα positions\nRMSD = {sup_partner.rms:.2f} Å",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
    )
    ax.axis("off")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.955),
        frameon=False,
        fontsize=10,
    )

    fig.suptitle(
        "AlphaFold Motif Models: H-NS:H-NS vs Alc:H-NS",
        y=1.06,
        fontsize=13,
        weight="bold",
    )
    fig.savefig(OUTPNG, dpi=300, bbox_inches="tight", facecolor="white")
    print(OUTPNG)


if __name__ == "__main__":
    main()
