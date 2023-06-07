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

@st.cache_resource
def process_csv(df, df2):
    df = df.copy()
    cst =['ligand', 'collection_key', 'scenario']
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

@st.cache_resource
def merge_csv(csv_list):

    to_drop = ['score_average', 'score_min', 'attr_smi', 'attr_heavy_atom_count']

    if len(csv_list) == 1:
        return pd.read_csv(csv_list[0])
    else:
        df_list = []
        for i in range(len(csv_list)):
            df_list.append(pd.read_csv(csv_list[i]).drop(columns=to_drop))
        df_csv_concat = pd.concat(df_list)
        df_csv_concat = df_csv_concat.groupby(['ligand', 'collection_key', 'scenario'], as_index=False).sum(min_count=1)
        return df_csv_concat

st.title("VF-Stream")

st.title("Data loading")
upfile1 = st.file_uploader("Select your docking result file(s)", accept_multiple_files=True)
upfile2 = st.file_uploader("Select your UniProt mapping file")
    
st.title("Results")

if upfile1 and upfile2 is not None:
    df = merge_csv(upfile1)
    df2 = pd.read_csv(upfile2)
    df_final, upid_to_gene, gene_to_upid, pdbid_to_gene = process_csv(df, df2)

    try:
        df_final, upid_to_gene, gene_to_upid, pdbid_to_gene = process_csv(df, df2)
        df_final = df_final.set_index('ligand')
        st.write("Successfully processed results")

    except:
        st.write("Failed to process the docking results file")

    expander = st.expander("Processed results")        
    expander.dataframe(df_final)
    st.download_button(
            label="Download CSV",
            data=df_final.to_csv().encode('utf-8'),
            file_name='final.csv',
            mime='text/csv',
        )
    
    scenarios = list(df_final['scenario'].unique())
    
    #########################AVERAGE SCORE DIST#########################

    st.header("Average score distribution")
    
    scenario_options = st.multiselect(
    'Selected docking scenario', scenarios, scenarios[0])

    if len(scenario_options) > 1:

        scenario = df_final.copy().reset_index()
        scenario['ligand'] = scenario['ligand'] + '_' + scenario['scenario']
        scenario = scenario.set_index('ligand')
    else:
        scenario = df_final.copy()

    scenario = scenario.loc[scenario['scenario'].isin(scenario_options)]

    df2 = scenario.loc[:, scenario.columns.str.endswith('_avg')].transpose()
    df2 = df2.rename(index = lambda x: x.replace('_avg', ''))
    # df2['UniProt ID'] = [f'<a target="_blank" href="{UPLINK}{gene_to_upid[i]}/entry">{gene_to_upid[i]}</a>' for i in df2.index.values]
    df2['UniProt ID'] = [gene_to_upid[i] for i in df2.index.values]
    df2.index.name = 'Gene Name'

    col_list = list(df2.columns)
    col_list.pop(col_list.index('UniProt ID'))
    
    ligand_options = st.multiselect(
    'Selected ligand ID', col_list, col_list[0])
    
    if len(ligand_options) != 0:
        sns.set_style('darkgrid')
        fig = sns.histplot(df2[ligand_options])
        fig.set(xlabel='Binding score (kcal/mol)')
        st.pyplot(fig.figure)

    #########################LIGAND-CENTRIC ANALYSIS#########################

    st.header("Ligand-centric analysis")
    scenario_option = st.selectbox('Selected docking scenario', scenarios)

    scenario = df_final.loc[df_final['scenario'] == scenario_option]

    df2 = scenario.loc[:, scenario.columns.str.endswith('_avg')].transpose()
    df2 = df2.rename(index = lambda x: x.replace('_avg', ''))
    # df2['UniProt ID'] = [f'<a target="_blank" href="{UPLINK}{gene_to_upid[i]}/entry">{gene_to_upid[i]}</a>' for i in df2.index.values]
    df2['UniProt ID'] = [gene_to_upid[i] for i in df2.index.values]
    df2.index.name = 'Gene Name'

    col_list = list(df2.columns)
    col_list.pop(col_list.index('UniProt ID'))

    ligand_option = st.selectbox('Selected ligand ID', col_list)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Top targets")
        cutoff = st.slider('Select range', 0, len(df2.index), 5)
        top = df2[['UniProt ID', ligand_option]].sort_values(by= [ligand_option])
        top = top.rename(columns= {f'{ligand_option}' : 'Score'}) #Binding affinity (kcal/mol)
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

            genre = st.radio("Input type", ('Gene Name', 'UniProt ID'))

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
                    st.write("Gene name/UniProt ID not found")

    st.subheader("True targets ranking")
    upfile = st.file_uploader("Select your ground truth target file")
    if upfile is not None:
        st.write(f"Selected the file: {upfile.name}")

        genre = st.selectbox("True target input type", ('Gene Name', 'UniProt ID'))

        try:
            tmp = pd.read_csv(upfile)

            tt_list = list(tmp[genre].loc[tmp[ligand_option] == 1])
            res = top.loc[top[genre].isin(tt_list)]
            st.dataframe(res, use_container_width=True)

            genre_to_label = dict(zip(tmp[genre], tmp[ligand_option]))
            gt = top[genre].map(genre_to_label)
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
        except:
            st.markdown(f"Failed to process the ground truth target file with input type :red[{genre}]")

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

