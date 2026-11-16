# Molecular Fingerprinting & Similarity Analysis Using RDKit

Independent educational cheminformatics portfolio project.

## What it demonstrates
- Morgan fingerprints (radius 2; ECFP4-style) with 2048 bits
- MACCS structural keys
- Pairwise Tanimoto similarity
- Similarity matrices and heatmaps
- Molecular structure grid
- CSV summaries and PDF report

## Files
- `data/molecule_dataset.csv`: 20 example molecules with SMILES
- `Molecular_Fingerprinting_Similarity_Analysis.ipynb`: Google Colab/Jupyter notebook
- `molecular_similarity_analysis.py`: runnable Python script
- `requirements.txt`: dependencies
- `PROJECT_GUIDE.md`: method and interpretation guide
- `output/`: generated results (created after running)

## Google Colab
1. Upload the notebook and `molecule_dataset.csv`.
2. Run `!pip install -q rdkit pandas numpy matplotlib seaborn reportlab`.
3. Upload the CSV when prompted in the notebook if it is not in `data/`.
4. Run all cells. The final cell downloads a ZIP of generated outputs.

## Local execution
```bash
pip install -r requirements.txt
python molecular_similarity_analysis.py
```

## Interpretation
Tanimoto scores range from 0 to 1 for these binary fingerprints. A score of 1 means identical fingerprint bit vectors. Similarity depends on fingerprint type and settings.

## Limitations
A high fingerprint similarity score does not prove similar biological activity, potency, safety, binding affinity, or clinical usefulness. This is an educational demonstration, not experimental validation or a validated drug-discovery workflow.

**Project type:** Independent educational portfolio project; not client work.
