import streamlit as st
import pandas as pd
import requests
import os

url = "https://raw.githubusercontent.com/HPPProteome/HPP_TargetList/main/lists/Supplementary_Table_3.xlsx"
gene_file = "Supplemental_Table_3.xlsx"


if not os.path.exists(gene_file):
    print("Downloading Supplemental Table 3 from GitHub")
    response = requests.get(url)


    if response.status_code == 200:
        with open(gene_file, "wb") as file:
            file.write(response.content)
        print(f"File downloaded successfully and saved as '{gene_file}'")
    else:
        print(f"Failed to download file. HTTP status code: {response.status_code}")

else:
    print("Data file found")

gene_data = pd.read_excel(gene_file, engine="openpyxl")

uniprot_search = "https://www.uniprot.org/uniprotkb//entry"

for index, row in gene_data.iterrows():
    id = row["UniProtKB id"]
    link = f"https://www.uniprot.org/uniprotkb/{id}/entry"
    gene_data.at[index, "UniProtKB id"] = "[%s](%s)" % (id, link)

markdown_table = gene_data.to_markdown(index=False)
st.title("76 PE5 Genes")
st.markdown("Below is the full table. Click on any of the gene identifiers to go to its corrosponding database entry.")

st.markdown(markdown_table, unsafe_allow_html=True)
