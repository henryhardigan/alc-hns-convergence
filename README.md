# Alc-HNS Convergence

Repository for the manuscript, figures, motif alignment, AlphaFold motif-model files, and supporting analysis outputs for:

**Unusually specific convergence of bacteriophage T4 Alc on the E. coli/Shigella H-NS dimerization-site-2 state**

## Overview

This project documents an unusually specific local convergence between a 30-residue segment of bacteriophage T4 Alc and the E. coli/Shigella H-NS dimerization-site-2 state. The main biological interpretation advanced in the manuscript is that this convergence is consistent with competitive interference with H-NS filament assembly in enterobacterial hosts.

The repository includes:

- a direct Alc/H-NS sequence register and family alignment
- proteome-scale motif rarity results against reviewed Swiss-Prot
- family-restricted null-model summaries
- AlphaFold-Multimer motif models for Alc:H-NS and H-NS:H-NS
- a partner-only structural overlay supporting related geometry

## Repository Layout

- `manuscript/`
  - [`alc_hns_convergence_manuscript.docx`](manuscript/alc_hns_convergence_manuscript.docx)
  - [`alc_hns_convergence_manuscript.pdf`](manuscript/alc_hns_convergence_manuscript.pdf)
- `figures/`
  - `figure1_alc_hns_alignment.pdf`
  - `figure2_hns_family_alignment.pdf`
  - `figure3_structure_overlay.png`
- `data/`
  - `alc_hns_family_alignment.fasta`
- `models/`
  - AlphaFold Server job requests, confidence summaries, and top CIF models for `alc_hns/` and `hns_hns/`
- `analysis/`
  - markdown and TSV summaries of motif rarity, family nulls, partner-overlay nulls, and dimer-state RMSD
- `scripts/`
  - core analysis scripts copied from the original working project

## Reproducibility Notes

- `scripts/make_alc_hns_structure_figure.py` is repo-relative and directly runnable from this repository.
- `scripts/compare_hns_dimer_rmsd.py` is also repo-relative and writes its output into `analysis/`.
- Several other scripts assume access to a large reviewed Swiss-Prot FASTA collection that is not bundled here.
- The `analysis/` summaries capture the main numerical results used in the manuscript even when the full background dataset is not present.

## Suggested Citation

Hardigan HL. *Unusually specific convergence of bacteriophage T4 Alc on the E. coli/Shigella H-NS dimerization-site-2 state*. Repository version and Zenodo archive.
