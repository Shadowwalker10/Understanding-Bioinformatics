import pickle
import streamlit as st
from PIL import Image
from typing import List
from moleculeprocessor import MoleculeProcessor
import subprocess
import base64
import os
import pandas as pd
import numpy as np

def filedownload(df):
    csv = df.to_csv(index = False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href = "data:file/csv;base64,{b64}" download = "predictions.csv">Download csv File</a>'
    return href

def predictor(input_data:pd.DataFrame):
    predictions = pickle.load(file = open("./predictor.pkl", "rb")).predict(input_data)
    st.header("**Prediction Output**")
    final_prediction = pd.concat([pd.Series(load_data["Chembl ID"],name = "Chembl ID"),
                                 pd.Series(predictions, name  = "Predicted pEC50")], axis = 1)
    st.write(final_prediction)
    st.markdown(filedownload(final_prediction), unsafe_allow_html = True)

## Image
image = Image.open("logo.png")
st.image(image, use_column_width = True)

## Page title
st.markdown("""
# Bioactivity Prediction App

This activity allows the user to predict the bioactivity towards inhibiting the
Acetylcholinesterase Enzyme.
""")



## Including a sidebar

with st.sidebar.header("Upload Your CSV Data"):
    uploaded_file = st.file_uploader(label = "Upload your input file",type = ['txt'])
    st.sidebar.markdown("""
    [Example Input File]("./example_file.txt")
    """)

## Prediction
if st.sidebar.button("Predict"):
    load_data = pd.read_table(uploaded_file, sep = " ",header = None, names = ["Chembl ID", "Smiles"])
    st.header("**Original Input Data**")
    st.write(load_data)

    with st.spinner("Calculating Fingerprint..."):
        ## Fingerprint Calculation
        processor = MoleculeProcessor()
        fingerprint = processor.process_smile(load_data["Smiles"].values)
        ## Checking if infinity present in any row
        rows_with_inf = fingerprint.applymap(np.isinf).any(axis = 1)
        final_fingerprint = fingerprint[~rows_with_inf]
        final_fingerprint.to_csv("./pubchem.csv", index = False)

    ## Display the calculated descriptor
    st.header('**Calculated Molecular Descriptor**')
    desc = pd.read_csv("./pubchem.csv")
    st.write(desc)
    st.write(desc.shape)
    

    ## Apply PCA for dimensionality reduction
    pca_fingerprint = pd.DataFrame(pickle.load(open("./pca.pkl", "rb")).transform(final_fingerprint))
    pca_fingerprint.to_csv("./pca_descriptors.csv", index = False)
    
    st.header("**Descriptors after dimensionality reduction**")
    desc = pd.read_csv("./pca_descriptors.csv")
    st.write(desc)
    st.write(desc.shape)

    predictor(desc)

else:
    st.info("Upload Input Data!!!")