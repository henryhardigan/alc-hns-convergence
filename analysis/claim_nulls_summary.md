# Alc/H-NS Claim-Completion Nulls

## Exact family-branch null

- Reviewed H-NS/StpA family proteins: `20`
- Host branch proteins (E. coli/Shigella H-NS only): `4`
- Exact primary-motif hits within family: `4`
- Exact hits observed in host branch: `4`
- Uniform-placement branch p-value: `0.000206`

## Composition-preserving branch null

- Shuffles of the Alc 30-mer: `20000`
- Observed best constrained family score: `51.5`
- Observed top entries: `HNS_ECO57,HNS_ECOL6,HNS_ECOLI,HNS_SHIFL`
- Observed best non-host score: `49.5`
- Observed host-vs-nonhost margin: `2.0`

Empirical frequencies over shuffled Alc sequences:
- Any host-branch top hit: `1165/20000`
- All top hits confined to host branch: `640/20000`
- Margin >= observed (`2.0`): `348/20000`
- Any host-branch top hit AND margin >= observed: `348/20000`
- Best score >= observed (`51.5`): `0/20000`

Interpretation:
- The exact family-branch null asks whether all exact family hits would fall into the E. coli/Shigella H-NS branch by chance if family placement were uniform.
- The shuffled-Alc branch null asks whether Alc-like composition alone tends to pick the host branch as the best constrained family match, and whether it does so with the observed score margin over non-host H-NS/StpA proteins.
