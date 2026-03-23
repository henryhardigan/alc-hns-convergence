# Alc/H-NS Partner-Overlay Null

Null design:
- Enumerated every order-preserving 28-position alignment between the 30-residue Alc partner segment and the 29-residue H-NS partner segment.
- This is equivalent to dropping 2 Alc positions and 1 H-NS position, then pairing the remaining residues in order.
- Total alignments considered: `12615`.
- Each alignment was scored by BLOSUM62 and evaluated by best-fit CA RMSD of the partner chains only.

Observed preferred register:
- BLOSUM62 score: `67.0`
- Partner-only RMSD: `3.440 A`

Extremes over the full alignment space:
- Best (lowest) partner-only RMSD among all 12,615 alignments: `2.980 A` (drops Alc 1,21; H-NS 1)
- Best BLOSUM62 score among all 12,615 alignments: `67.0` (partner-only RMSD `3.440 A`)

Observed-register rank by partner-only RMSD:
- Against all alignments: `2020/12615` alignments have RMSD <= observed
- Among alignments with BLOSUM62 score >= observed: `1/1`
- Among alignments within 2 BLOSUM62 points of observed: `4/5`

Interpretation:
- This null asks whether the preferred convergence register produces an unusually good structural overlay relative to the full space of order-preserving 28-position mappings.
- The chemistry-filtered subsets are the most relevant comparison for a convergent-mimicry claim, because they hold sequence plausibility roughly constant while testing structural specificity.
