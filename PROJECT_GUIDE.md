# Project Guide

## Aim
Compare chemical structures computationally using Morgan fingerprints, MACCS keys and Tanimoto similarity.

## Workflow
1. Read molecule names and SMILES from CSV.
2. Parse SMILES with RDKit and flag invalid structures.
3. Generate Morgan fingerprints (radius 2, 2048 bits) and MACCS keys.
4. Calculate pairwise Tanimoto coefficients.
5. Rank unique molecule pairs.
6. Plot similarity heatmaps and molecule structures.
7. Export CSV files and a PDF summary.

## Tanimoto coefficient
For binary fingerprints A and B:

T(A,B) = c / (a + b - c)

where a and b are the number of on-bits in A and B, and c is the number of shared on-bits.

## Discussion questions
- What information does a molecular fingerprint encode?
- How does Morgan radius affect the local atom environment?
- How do Morgan fingerprints differ from MACCS keys?
- Why is Tanimoto similarity useful for binary fingerprints?
- Why does structural similarity not guarantee similar biological activity?

## Caution
Results depend on dataset quality and fingerprint parameters. Do not interpret similarity as evidence of activity, safety, binding, or therapeutic equivalence.
