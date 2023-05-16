import streamlit as st
import pandas as pd
import numpy as np
import rdkit
# import py3Dmol
from sklearn.metrics import roc_auc_score, roc_curve, auc
import matplotlib.pyplot as plt
import altair as alt
from io import BytesIO
import base64
import re
# st.set_page_config(layout="wide")
from io import StringIO
import seaborn as sns
import requests
import streamlit.components.v1 as components
import urllib
import io
import biotite.structure.io.pdbx as pdbx
import webbrowser
UPLINK = "https://www.uniprot.org/uniprotkb/"


st.title("VF-Proteome scrapwork")

@st.cache_resource
def process_csv(df, df2):
    df = df.copy()
    cst =['ligand', 'collection_key', 'scenario', 'score_average', 'score_min', 'attr_smi', 'attr_heavy_atom_count']
    receptors = list(set(df.columns) - set(cst))

    upid_to_gene = dict(zip(df2['UniProtID'], df2['GeneName']))
    gene_to_upid = dict(zip(df2['GeneName'], df2['UniProtID']))
    pdbid_to_gene = dict(zip(df2['PDBID'], df2['GeneName']))

    gene_names = set()
    for r in receptors:
        l = r.split('_replica_')
        num = 0
        if len(l) != 1:
            num = l[-1]
        name = l[0]
        pattern = 'AF-(.+?)-F1'
        match = re.search(pattern, name)
        if match:
            try:
                df = df.rename(columns={r: upid_to_gene[match.group(1)] + '_' + str(num)})
                gene_names.add(upid_to_gene[match.group(1)])
            except:
                print(f"Warning: No correspondance to any gene name was found for receptor named: {name}")
                gene_names.add(name)
        else:
            try:
                df = df.rename(columns={r: pdbid_to_gene[name] + '_' + str(num)})
                gene_names.add(pdbid_to_gene[name])
            except:
                print(f"Warning: No correspondance to any gene name was found for receptor named: {name}")
                gene_names.add(name)

    for gn in gene_names:
        df[gn + '_avg'] = df.filter(regex='^'+gn).mean(axis=1)
    
    df.to_csv('final.csv', index=False)
    df_final = pd.read_csv('final.csv')

    return df_final, upid_to_gene, gene_to_upid, pdbid_to_gene


st.title("Data loading")
upfile1 = st.file_uploader("Select your docking result file")
upfile2 = st.file_uploader("Select your UniProt mapping file")
    

st.title("Results")

if upfile1 and upfile2 is not None:
    df = pd.read_csv(upfile1)
    df2 = pd.read_csv(upfile2)

    try:
        df_final, upid_to_gene, gene_to_upid, pdbid_to_gene = process_csv(df, df2)
        df_final = df_final.set_index('ligand')
        st.write("Successfully processed results")

    except:
        st.write("Failed to process the docking results file")

    st.download_button(
            label="Download CSV",
            data=df_final.to_csv().encode('utf-8'),
            file_name='final.csv',
            mime='text/csv',
        )
    
    df2 = df_final.loc[:, df_final.columns.str.endswith('_avg')].transpose()
    df2 = df2.rename(index = lambda x: x.replace('_avg', ''))
    # df2['UniProt ID'] = [f'<a target="_blank" href="{UPLINK}{gene_to_upid[i]}/entry">{gene_to_upid[i]}</a>' for i in df2.index.values]
    df2['UniProt ID'] = [gene_to_upid[i] for i in df2.index.values]
    df2.index.name = 'Gene Name'

    col_list = list(df2.columns)
    col_list.pop(col_list.index('UniProt ID'))

    st.header("Average score distribution")
    options = st.multiselect(
    'Selected ligand ID', col_list, col_list[0])
    
    if len(options) != 0:
        sns.set_style('darkgrid')
        fig = sns.histplot(df2[options])
        fig.set(xlabel='Binding score (kcal/mol)')
        st.pyplot(fig.figure)
    
    st.header("Ligand-centric analysis")
    option = st.selectbox('Selected ligand ID', col_list)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Top targets")
        cutoff = st.slider('Select range', 0, len(df2.index), 5)
        top = df2[['UniProt ID', option]].sort_values(by= [option])
        top = top.rename(columns= {f'{option}' : 'Score'}) #Binding affinity (kcal/mol)
        top.reset_index(inplace=True)
        top.index.name = 'Rank'
        st.dataframe(top.head(cutoff), use_container_width=True)

        txt = ''
        for g in top['Gene Name'].head(cutoff):
            url = UPLINK + gene_to_upid[g]
            txt += f'[{g}]({url}) '
        expander = st.expander("UniProt links")        
        expander.write(txt)
    
    with col2:

        st.subheader("Query target ranking")
        with st.form("entry_form"):
            
            g = st.text_input("Type the name of genes/UniProt ID you screened (separated by commas):", f"{top['Gene Name'][0]}").upper()
            g_list = g.split(',')
            g_list = [x.strip(' ') for x in g_list]

            genre = st.radio("Input type", ('Gene name', 'UniProt ID'))

            submitted = st.form_submit_button("Query")
            if submitted:
                try:
                    if genre == 'Gene Name':
                        res = top.loc[top['Gene Name'].isin(g_list)]
                        st.dataframe(res, use_container_width=True)
                        
                        txt = ''
                        for g in g_list:
                            url = UPLINK + gene_to_upid[g]
                            txt += f'[{g}]({url}) '
                        expander1 = st.expander("UniProt links")        
                        expander1.write(txt)

                    elif genre == 'UniProt ID':
                        res = top.loc[top['UniProt ID'].isin(g_list)]
                        st.dataframe(res, use_container_width=True)

                        txt = ''
                        for g in g_list:
                            url = UPLINK + g
                            txt += f'[{g}]({url}) '
                        expander1 = st.expander("UniProt links")        
                        expander1.write(txt)


                except:
                    st.write("Gene name not found")

    st.subheader("True targets")
    upfile = st.file_uploader("Select your ground truth target file")
    if upfile is not None:
        st.write(f"Selected the file: {upfile.name}")
        
        tmp = pd.read_csv(upfile)
        gene_to_label = dict(zip(tmp['Gene Name'], tmp[option]))

        gt = top['Gene Name'].map(gene_to_label)
        norm = (top['Score']-top['Score'].max())/(top['Score'].min()-top['Score'].max())

        st.subheader("ROC-AUC score")

        fpr, tpr, _ = roc_curve(gt, norm)
        auc = auc(fpr, tpr)

        roc_df = pd.DataFrame({'False Positive Rate': fpr, 'True Positive Rate': tpr})

        fig = plt.figure(figsize=(8, 6))
        sns.lineplot(x='False Positive Rate', y='True Positive Rate', data=roc_df, label=f'AUC = {auc:.2f}')
        plt.plot([0, 1], [0, 1], 'r--')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend(loc='lower right')
        st.pyplot(fig)

st.title("TODO")
"Sequence similarity -- Ligand info (?) "

# -- True target module: ROC AUC (with interactive cutoff selection ?) is this done well ? Ask Alex


def get_alphafold_data(entry_name):
        """downloads alphafold .cif file for the input protein and extracts sequence and model confidence array for each amino acid position
        Args:
            entry_name (str): uniprot entry name of the protein
        Returns:
            tuple of sequence and confidence array: the protein sequence from the alphafold website with the model confidence array. 
                            If the entry is not available on the server, "x" insted of the sequence and an array of [0] as the confidence is returned
        """
        entry = entry_name
        
        connection = urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/files/AF-{entry}-F1-model_v4.cif")
        databytes = connection.read()
        connection.close()
        cif_txt = databytes.decode("utf8")

        f = io.StringIO(cif_txt)
        cif = pdbx.PDBxFile.read(f)
        
        confidence = pd.DataFrame(cif.get_category("ma_qa_metric_local")).metric_value.astype(float).values
        sequence = cif.get_category("entity_poly")["pdbx_seq_one_letter_code"]

        return sequence, confidence, cif

