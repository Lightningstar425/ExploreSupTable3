import sys
print(sys.path)
import streamlit as st
import pandas as pd
import requests
import os
import re


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

raw_gene_data = pd.read_excel(gene_file, engine="openpyxl")
gene_data = pd.read_excel(gene_file, engine="openpyxl")

uniprot_search = "https://www.uniprot.org/uniprotkb/{}/entry"
ensemble_search = "https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={}"
genecard_searc = "https://www.genecards.org/cgi-bin/carddisp.pl?gene={}"
peptideAtals_search = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=592&protein_name={}&action=QUERY"
HPA_search = "https://www.proteinatlas.org/{}"
alphaFold_search = "https://alphafold.ebi.ac.uk/entry/{}"
NCBI_search = "https://www.ncbi.nlm.nih.gov/protein/{}/"
pub_med_search = "https://pubmed.ncbi.nlm.nih.gov/?term={}"

#Set links
for index, row in gene_data.iterrows():
    id = row["UniProtKB id"]
    uniprot_link = f"https://www.uniprot.org/uniprotkb/{id}/entry"
    ensemble_link = f"https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={row["ENSEMBL id"]}"
    genecard_link = f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={row["UniProtKB symbol"]}"
    peptideAtlas_link = f"https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=592&protein_name={id}&action=QUERY"
    HPA_link = f"https://www.proteinatlas.org/{row["ENSEMBL id"]}"
    alphaFold_link = f"https://alphafold.ebi.ac.uk/entry/{id}"
    NCBI_link = f"https://www.ncbi.nlm.nih.gov/protein/{id}/"

    gene_data.at[index, "UniProtKB id"] = "[%s](%s)" % (id, uniprot_link)
    gene_data.at[index, "ENSEMBL id"] = "[%s](%s)" % (row["ENSEMBL id"], ensemble_link)
    gene_data.at[index, "UniProtKB symbol"] = "[%s](%s)" % (row['UniProtKB symbol'], genecard_link)
    gene_data.at[index, "PeptideAtlas status"] = "[%s](%s)" % (row['PeptideAtlas status'], peptideAtlas_link)
    gene_data.at[index, "HPA tissue-specific RNA evidence"] = "[%s](%s)" % (row['HPA tissue-specific RNA evidence'], HPA_link)
    gene_data.at[index, "AlphaFold predicted structure average pLDDT"] = "[%s](%s)" % (row['AlphaFold predicted structure average pLDDT'], alphaFold_link)
    gene_data.at[index, "NCBI Entrez Gene summary"] = "[%s](%s)" % (row['NCBI Entrez Gene summary'], NCBI_link)

#markdown_table = gene_data.to_markdown(index=False)
st.set_page_config(layout="wide")

st.title("76 PE5 Gene Table")

st.header("Filter table")

#Create 7 colunms for sorting: 
cols = st.columns(6)

with cols[0]:
    st.markdown("Assessment Category")
    show_CPC = st.checkbox("CPC", value=True)
    show_MoFLPC = st.checkbox("MoFLPC", value=True)
    show_TxEOEP = st.checkbox("TxEOEP", value=True)
    show_LNC = st.checkbox("LNC", value=True)
    show_PPC = st.checkbox("PPC", value=True)
    show_Both = st.checkbox("Both", value=True)

selected_assessment = []
if show_CPC:
    selected_assessment.append("CPC")
if show_MoFLPC:
    selected_assessment.append("MoFLPC")
if show_TxEOEP:
    selected_assessment.append("TxEOEP")
if show_LNC:
    selected_assessment.append("LNC")
if show_PPC:
    selected_assessment.append("PPC")
if show_Both:
    selected_assessment.append("Both")

with cols[1]:
    st.markdown("New Suggested PE")
    show_PE1 = st.checkbox("PE1", value=True)
    show_PE2 = st.checkbox("PE2", value=True)
    show_PE3 = st.checkbox("PE3", value=True)
    show_PE5 = st.checkbox("PE5", value=True)

selected_pe = []
if show_PE1:
    selected_pe.append('1')
if show_PE2:
    selected_pe.append('2')
if show_PE3:
    selected_pe.append('3')
if show_PE5:
    selected_pe.append('5')


with cols[2]:
    st.markdown("Present In:")
    show_gencode = st.checkbox("GENCODE", value=True)
    show_peptideAtlas = st.checkbox("Peptide Atlas", value=False)
    show_iTasser = st.checkbox("I-TASSER C>0", value=False)
    show_pfind = st.checkbox("pFind PE5 list", value=False)

gencode = []
if show_gencode:
    gencode.append("yes")
else:
    gencode.append("no")

peptideAtlas = []
if show_peptideAtlas:
    peptideAtlas.append("yes")
else:
    peptideAtlas.append("no")

iTasser = []
if show_iTasser:
    iTasser.append("yes")
else:
    iTasser.append("no")

p_find = []
if show_pfind:
    p_find.append("yes")
else:
    p_find.append("no")

with cols[3]:
    st.markdown("PeptideAtlas status:")
    show_canonical = st.checkbox("Canonical", value=True)
    show_weak = st.checkbox("Weak", value=True)
    show_subsumed = st.checkbox("Subsumed", value=True)
    show_margDis = st.checkbox("Marginally Distinguished", value=True)
    show_PAnotDetected = st.checkbox("Not Detected", value=True)
    show_PAna = st.checkbox("N/A", value=True)

selected_PA = []
if show_canonical:
    selected_PA.append("canonical")
if show_weak:
    selected_PA.append("weak")
if show_subsumed:
    selected_PA.append("subsumed")
if show_margDis:
    selected_PA.append("marginally dist")
if show_PAnotDetected:
    selected_PA.append("not detected")
if show_PAna:
    selected_PA.append("-")


with cols[4]:
    st.markdown("HPA tissue-specific RNA evidence:")
    show_tissE = st.checkbox("Tissue enriched", value=True)
    show_groupE = st.checkbox("Group enriched", value = True)
    show_lowTissE = st.checkbox("Low Tissue Specificity", value=True)
    show_HPAnotDetected = st.checkbox("Not Detected ", value=True)
    show_HPAna = st.checkbox("N/A ", value=True)

selected_HPA = []
if show_tissE:
    selected_HPA.append("Tissue enriched")
if show_groupE:
    selected_HPA.append("Group enriched")
if show_lowTissE:
    selected_HPA.append("low tissue specificity")
if show_HPAnotDetected:
    selected_HPA.append("not detected")
if show_HPAna:
    selected_HPA.append("-")


with cols[5]:
    st.markdown("AlphaFold predicted structure average pLDDT:")
    show_Vhigh = st.checkbox("Very High", value=True)
    show_high = st.checkbox("High", value=True)
    show_low = st.checkbox("Low", value=True)
    show_Vlow = st.checkbox("Very Low", value=True)
    show_AFna = st.checkbox("N/A  ", value=True)

selected_AF = []
if show_Vhigh:
    selected_AF.append("Very High")
if show_high:
    selected_AF.append('High')
if show_low:
    selected_AF.append("Low")
if show_Vlow:
    selected_AF.append("Very Low")
if show_AFna:
    selected_AF.append("-")

search = st.text_input("Search by ENSG number or ID:")
gene_data['Suggested New PE'] = gene_data['Suggested New PE'].astype(str)
gene_data['PeptideAtlas status'] = gene_data['PeptideAtlas status'].astype(str)

pa_mask = gene_data['PeptideAtlas status'].apply(lambda x: any(phrase in str(x) for phrase in selected_PA))
hpa_mask = gene_data['HPA tissue-specific RNA evidence'].apply(lambda x: any(phrase in str(x) for phrase in selected_HPA))
af_mask = gene_data['AlphaFold predicted structure average pLDDT'].apply(lambda x: any(phrase in str(x)for phrase in selected_AF))

wanted_data = gene_data[gene_data['Assess- ment Category'].isin(selected_assessment) & gene_data['Suggested New PE'].isin(selected_pe) & gene_data["GEN CODE?"].isin(gencode) &
                        gene_data['Peptide Atlas?'].isin(peptideAtlas) & gene_data['I-TASSER C>0?'].isin(iTasser) &
                        gene_data['pFind PE5 list?'].isin(p_find) &
                        pa_mask & hpa_mask & (gene_data["UniProtKB id"].str.contains(search, case=False, na=False) |
                        gene_data["ENSEMBL id"].str.contains(search, case=False, na=False))]


#wanted_data = gene_data[gene_data["UniProtKB id"].str.contains(search, case=False, na=False) |
 #                       gene_data["ENSEMBL id"].str.contains(search, case=False, na=False) &
  #                      gene_data['Suggested New PE'].isin(selected_pe)]

col_width = 20  # number of characters per column
data_headers = gene_data.columns.tolist() 

header = "| " + " | ".join([str(h).ljust(col_width) for h in data_headers]) + " |"
divider = "| " + " | ".join(["-" * col_width] * len(data_headers)) + " |"

markdown_table = "\n".join([header, divider])

for index, row in wanted_data.iterrows():
        


        row_values = [str(val).ljust(col_width) for val in row.tolist()]
        

        data_row = "| " + " | ".join(row_values) + " |"

        markdown_table += "\n" + data_row

st.markdown(markdown_table)

                  

