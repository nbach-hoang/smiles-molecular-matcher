import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Draw

def calculate_similarity(target_smiles, candidate_smiles):
    try:
        target_mol = Chem.MolFromSmiles(target_smiles)
        candidate_mol = Chem.MolFromSmiles(candidate_smiles)

        if not target_mol or not candidate_mol:
            return None

        target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2, nBits=2048)
        candidate_fp = AllChem.GetMorganFingerprintAsBitVect(candidate_mol, 2, nBits=2048)

        similarity = DataStructs.TanimotoSimilarity(target_fp, candidate_fp)
        return round(similarity, 4)
    except Exception:
        return None

st.set_page_config(page_title="Virtual Screening Tool", layout="centered")
st.title("Drug Discovery Tool: Molecular Similarity Matcher")
st.write("Screen chemical libraries against a target drug using Morgan Fingerprints & Tanimoto Coefficients.")

with st.sidebar:
    st.header("Quick Testing Data")
    st.write("**Target Molecule (Ibuprofen):**")
    st.code("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
    st.write("**Candidate Molecules to Paste:**")
    st.code("CC(C(=O)O)C1=CC=CC(=C1)C(=O)C2=CC=CC=C2\nCC(=O)OC1=CC=CC=C1C(=O)O\nCCO")
    st.caption("Top is Ketoprofen (similar), middle is Aspirin (less similar), bottom is Ethanol (not similar).")

# User Inputs
target_input = st.text_input("Enter Target Drug SMILES String:", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
candidates_input = st.text_area(
    "Enter Candidate SMILES Strings (One per line):",
    "CC(C(=O)O)C1=CC=CC(=C1)C(=O)C2=CC=CC=C2\nCC(=O)OC1=CC=CC=C1C(=O)O\nCCO"
)

if st.button("Run Virtual Screening Alignment"):
    # Split the textarea input into an individual list of strings
    candidate_list = [line.strip() for line in candidates_input.split("\n") if line.strip()]

    if not target_input or not candidate_list:
        st.warning("Please ensure both target and candidate fields contain valid chemical structures.")
    else:
        results = []

        for candidate in candidate_list:
            score = calculate_similarity(target_input, candidate)
            if score is not None:
                results.append({"Candidate SMILES": candidate, "Similarity Score": score})
            else:
                results.append({"Candidate SMILES": candidate, "Similarity Score": "Error (Invalid SMILES)"})

        results = sorted(results,
                         key=lambda x: x["Similarity Score"] if isinstance(x["Similarity Score"], float) else -1,
                         reverse=True)

        st.subheader("Screening Rankings:")
        for idx, item in enumerate(results):
            score = item["Similarity Score"]
            smiles = item["Candidate SMILES"]

            if isinstance(score, float):
                if score > 0.7:
                    st.success(f" **Rank #{idx + 1}: Score: {score}** — High Structural Match!\n`{smiles}`")
                elif score > 0.3:
                    st.info(f" **Rank #{idx + 1}: Score: {score}** — Moderate Match.\n`{smiles}`")
                else:
                    st.warning(f" **Rank #{idx + 1}: Score: {score}** — Weak Match.\n`{smiles}`")
            else:
                st.error(f" **Rank #{idx + 1}: {score}** for entry: `{smiles}`")
st.write('---')
st.title('Visual Structure')
if target_input:
    visual_mol = Chem.MolFromSmiles(target_input)

    if visual_mol:
        molecular_image = Draw.MolToImage(visual_mol, size=(300, 300))

        st.image(molecular_image, caption=f"Visual structure layout for: {target_input}", width=300)
    else:
        st.caption("Please enter a valid target SMILES string above to generate a 2D render.")
st.write('---')
st.title('Description')
st.write('A web application using Python programming language along with RDKit for automation of ligand-based virtual screening for drug discovery at initial stages,this application takes small molecules in form of SMILES and encodes them into 2048 bit long Morgan fingerprints (ECFP4) and calculates vector intersections by utilizing Tanimoto Similarity Coefficient.')
st.write('It ranks candidate compound collections based on the lead compound in real time and hence provCreated a web application using Python programming language along with RDKit for automation of ligand-based virtual screening for drug discovery at initial stages. This application takes small molecules in form of SMILES and encodes them into 2048 bit long Morgan fingerprints (ECFP4) and calculates vector intersections by utilizing Tanimoto Similarity Coefficient.')
