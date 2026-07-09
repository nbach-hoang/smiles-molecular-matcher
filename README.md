**Try The App Web Version Here: https://smiles-molecular-matcher.streamlit.app/
**
A web application using Python programming language along with RDKit for automation of ligand-based virtual screening for drug discovery at initial stages, this application takes small molecules in form of SMILES and encodes them into 2048 bit long Morgan fingerprints (ECFP4) and calculates vector intersections by utilizing Tanimoto Similarity Coefficient. 

It ranks candidate compound collections based on the lead compound in real time and hence provides an approach similar to high throughput chemical filtering in the pharmaceuticals domain. It is implemented using fully open-source technology stack using Streamlit.

**The program follows this pipeline:
**
User enters target SMILES
          │
          ▼
Convert to RDKit molecule
          │
          ▼
Generate Morgan fingerprint (2048-bit ECFP4)
          │
          ▼
Repeat for every candidate molecule
          │
          ▼
Calculate Tanimoto similarity
          │
          ▼
Sort by similarity score
          │
          ▼
Display ranked results
          │
          ▼
Render a 2D image of the target molecule

This is a simplified ligand-based virtual screening (LBVS) application. It does not predict biological activity or docking interactions; instead, it ranks candidate compounds based on 2D structural similarity to a known lead compound using Morgan (ECFP4) fingerprints and the Tanimoto coefficient, which is a common first-pass approach in early-stage drug discovery.
