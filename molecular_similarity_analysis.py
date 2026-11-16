from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem, MACCSkeys, Draw
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)
DATA = ROOT / "data" / "molecule_dataset.csv"
RADIUS, N_BITS = 2, 2048

def main():
    print("RDKit:", rdBase.rdkitVersion)
    df = pd.read_csv(DATA)
    records, invalid = [], []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["SMILES"]))
        if mol is None:
            invalid.append({"Name": row["Name"], "SMILES": row["SMILES"]})
        else:
            records.append((str(row["Name"]), str(row["SMILES"]), mol))
    if len(records) < 2:
        raise ValueError("At least two valid SMILES are required.")

    names = [r[0] for r in records]
    mols = [r[2] for r in records]
    morgan = [AllChem.GetMorganFingerprintAsBitVect(m, RADIUS, nBits=N_BITS) for m in mols]
    maccs = [MACCSkeys.GenMACCSKeys(m) for m in mols]

    def matrix(fps):
        return [DataStructs.BulkTanimotoSimilarity(fp, fps) for fp in fps]

    mm, km = matrix(morgan), matrix(maccs)
    morgan_df = pd.DataFrame(mm, index=names, columns=names)
    maccs_df = pd.DataFrame(km, index=names, columns=names)
    morgan_df.to_csv(OUT / "morgan_similarity_matrix.csv")
    maccs_df.to_csv(OUT / "maccs_similarity_matrix.csv")

    pairs = []
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            pairs.append([names[i], names[j], mm[i][j], km[i][j]])
    pair_df = pd.DataFrame(pairs, columns=["Molecule_1","Molecule_2","Morgan_Tanimoto","MACCS_Tanimoto"])
    pair_df = pair_df.sort_values("Morgan_Tanimoto", ascending=False)
    pair_df.to_csv(OUT / "pairwise_similarity_scores.csv", index=False)

    summary = pd.DataFrame({
        "Name": names,
        "Canonical_SMILES": [Chem.MolToSmiles(m) for m in mols],
        "Morgan_Bits": [N_BITS]*len(mols),
        "Morgan_On_Bits": [fp.GetNumOnBits() for fp in morgan],
        "MACCS_On_Bits": [fp.GetNumOnBits() for fp in maccs]
    })
    summary.to_csv(OUT / "fingerprint_summary.csv", index=False)
    if invalid:
        pd.DataFrame(invalid).to_csv(OUT / "invalid_smiles_rows.csv", index=False)

    for label, frame, cmap in [
        ("morgan", morgan_df, "viridis"), ("maccs", maccs_df, "mako")
    ]:
        plt.figure(figsize=(12,10))
        sns.heatmap(frame, cmap=cmap, vmin=0, vmax=1, square=True)
        plt.title(f"{label.upper()} Fingerprint Tanimoto Similarity")
        plt.xticks(rotation=75, ha="right")
        plt.yticks(rotation=0)
        plt.tight_layout()
        plt.savefig(OUT / f"{label}_similarity_heatmap.png", dpi=300)
        plt.close()

    img = Draw.MolsToGridImage(mols[:20], molsPerRow=4, subImgSize=(250,200),
                               legends=names[:20], useSVG=False)
    img.save(OUT / "molecule_structure_grid.png")

    styles = getSampleStyleSheet()
    doc = SimpleDocTemplate(str(OUT/"Molecular_Similarity_Analysis_Report.pdf"), pagesize=A4)
    story = [
        Paragraph("Molecular Fingerprinting & Similarity Analysis", styles["Title"]),
        Spacer(1,10),
        Paragraph(f"Valid molecules analysed: {len(records)}. Invalid SMILES skipped: {len(invalid)}.", styles["BodyText"]),
        Spacer(1,10),
        Paragraph("Top 10 molecule pairs by Morgan Tanimoto similarity", styles["Heading2"])
    ]
    table_data = [["Molecule 1","Molecule 2","Morgan Tanimoto"]]
    for _, row in pair_df.head(10).iterrows():
        table_data.append([row.Molecule_1, row.Molecule_2, f"{row.Morgan_Tanimoto:.3f}"])
    table = Table(table_data, repeatRows=1, colWidths=[145,145,100])
    table.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#263B5A")),
        ("TEXTCOLOR",(0,0),(-1,0),colors.white),
        ("GRID",(0,0),(-1,-1),0.4,colors.grey),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.whitesmoke,colors.lightgrey])
    ]))
    story += [table, Spacer(1,12), Image(str(OUT/"morgan_similarity_heatmap.png"), width=460, height=380),
              Spacer(1,10),
              Paragraph("Tanimoto similarity is a fingerprint-based structural comparison. "
                        "It does not establish biological activity, potency, safety or therapeutic equivalence.",
                        styles["BodyText"])]
    doc.build(story)
    print(pair_df.head(10).to_string(index=False))
    print("All outputs saved to:", OUT)

if __name__ == "__main__":
    main()
