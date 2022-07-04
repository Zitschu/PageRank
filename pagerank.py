"""
This script implements the PageRank and Topic-Sensitive PageRank for a biomedical dataset.
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import streamlit as st
from PIL import Image

def generate_personalization_dict(df):
    """ Generates a dict containing only certain proteins of a certain area.
    :param url: URL to a csv containing target proteins as a String.
    :return: Returns a dict where the key is the protein ID and the value is 1
    """
    data = df

    topic_dict = {key:1 for key in data['geneid']}

    return topic_dict


def convert_pr_to_df(pr_dict):
    """ Converts a PageRank dict to a DataFrame.
    :param pr_dict: dict containing the PageRank.
    :return: Returns a dataframe
    """
    return pd.DataFrame([{"geneid": geneid, "page_rank": page_rank} for geneid, page_rank in pr_dict.items()])


def get_missing_gene_ids(topic_dict, df):
    """ Compares the Genes in a topic_dict with those in a dataframe and creates a list containing Gene Ids that can't be found in the dataframe.
    :param topic_dict: dict containing the Gene IDs of a certain topic.
    :param df: Dataframe that is used for calculating the PageRank.
    :return: Returns a list of Ints
    """
    missing_genes = []

    for key in topic_dict.keys():
        if (key in df['Protein_1']) or (key in df['Protein_2']):
            l = 1
        else:
            missing_genes.append(key)
    
    return missing_genes


def pr_diff(df1, df2):
    """ Compares dataframes containing PageRanks for the same dataset.
    :param df1: Dataframe containing the basic PageRank.
    :param df2: Dataframe containing the personalized PageRank.
    :return: Returns a merged dataframe with a new column for the difference in PageRank
    """
    merged_df = df2.merge(df1, on="geneid", suffixes=("_pers", "_basic"))
    merged_df['page_rank_diff'] = (merged_df["page_rank_pers"] - merged_df["page_rank_basic"])
    merged_df['page_rank_diff_perc (%)'] = (merged_df['page_rank_diff']/merged_df["page_rank_basic"])*100
    
    return merged_df.head(10)


def main():
    df = pd.read_csv('datasets\PP-Pathways_ppi.csv', header=None, names=["Protein_1", "Protein_2"])
    topic_df = pd.read_csv('datasets\personalization\Fatty_acid_metabolism_Proteins.csv')
    topic_df = topic_df.dropna()
    topic_df['geneid'] = topic_df['geneid'].astype(int)

    topic_dict = generate_personalization_dict(topic_df)
    
    missing_genes = get_missing_gene_ids(topic_dict, df)

    for key in missing_genes:
        topic_dict.pop(key, None) # drop missing Ids from Topic Dict since they aren't relevant for the Topic Sensitive PageRank.

    ppi_graph = nx.from_pandas_edgelist(df, source="Protein_1", target="Protein_2")
    
    # pos = nx.spring_layout(ppi_graph, k=9*1/np.sqrt(len(ppi_graph.nodes())), iterations=20)
    # plt.figure(3, figsize=(20, 20))
    # nx.draw(ppi_graph, pos=pos)
    # nx.draw_networkx_labels(ppi_graph, pos=pos)
    # plt.savefig("graph.png", dpi=500)

    pr = nx.pagerank(ppi_graph)
    pr_df = convert_pr_to_df(pr)
    pr_df_sorted = pr_df.sort_values('page_rank', ascending=False)
    pr_df_sorted = pr_df_sorted.reset_index(drop=True)

    topic_sensitive_pr = nx.pagerank(ppi_graph, personalization=topic_dict)
    topic_sensitive_pr_df = convert_pr_to_df(topic_sensitive_pr)
    topic_sensitive_pr_df_sorted = topic_sensitive_pr_df.sort_values('page_rank', ascending=False)
    topic_sensitive_pr_df_sorted = topic_sensitive_pr_df_sorted.reset_index(drop=True)

    network_image = Image.open('graph.png')
    # network_image_resized = network_image.resize((500,500))

    st.header("PageRank on a Human Protein-Protein Interaction Network")
    st.subheader("Dataset(s)")
    st.markdown("There are two datasets used in this project. The ['Human Protein-Protein Interaction Network' from Stanford](http://snap.stanford.edu/biodata/datasets/10000/10000-PP-Pathways.html) provides the data for a Graph of different Human Protein interacting with each other. Additionally, we'll use a dataset containing Gene IDs for Proteins in [Fatty Acid Metabolism](https://pubchem.ncbi.nlm.nih.gov/pathway/Reactome:R-HSA-8978868#section=Proteins&fullscreen=true). This second dataset is necessary to calculate a Topic-Sensitive PageRank.")
    st.text("Human Protein-Protein Interaction Network Dataset")
    column1, column2 = st.columns(2)
    column1.dataframe(df)
    column2.image(network_image, caption="Human Protein-Protein Interaction Network", output_format="JPEG")
    st.markdown("Undoubtably, the graph is very cluttered because there are so many nodes (proteins) interacting with each other. To make it even more confusion, the edgelist allows for an undirected graph only. However in this graph undirected edges are replaced by two directional edges in both ways. I left the image in here, just to give a general idea how complex this network is. Don't expect to gain any information from it though. The table should be easier to understand.")
    st.text("Fatty Acid Metabolism Dataset")
    st.dataframe(topic_df)

    st.subheader("PageRank")
    st.markdown("The following shows the basic PageRank for the Human Protein-Protein Interaction Network. Note, that the Dataset provides data that only allows for an undirected graph. Therefore, this implementation converts undirected edges into bi-directional edges.")
    st.dataframe(pr_df)

    st.subheader("Topic-Sensitive PageRank")
    st.markdown("Topic-Sensitive PageRank calculates the PageRank for a subset of the Human Protein-Protein Interaction Network. I decided to look at Proteins that are involved in the Fatty Acid Metabolism and see what their importance in this subdomain looks like.")
    st.dataframe(topic_sensitive_pr_df)

    st.subheader("Comparison of Basic PageRank and Topic-Sensitive PageRank")
    st.markdown("We will take a look at the 10 highest PageRanks respectively. On the left hand-side, you can see the basic PageRank and on the right hand-side the Topic-Sensitive PageRank.")
    col1, col2 = st.columns(2)
    col1.dataframe(pr_df_sorted.head(10))
    col2.dataframe(topic_sensitive_pr_df_sorted.head(10))
    st.markdown("As you can see, there are huge differences in the basic PageRank and the Topic-Sensitive PageRank. Protein 351 has almost doubled its PageRank and has the highest PageRank in Fatty Acid Metabolism. Some other Proteins don't even appear in the Topic-Sensitive PageRank, e.g. 7514 and 801.")
    st.markdown("Below you find the two tables merged. The Proteins are sorted according to their Topic Sensitive PageRank. In Column 'page_rank_diff' you can find the difference between the PageRanks. Note that some PageRanks have a difference of up to 0.0035. The percentage change from the Basic PageRank to the Topic-Sensitive PageRank can be found in column 'page_rank_diff_perc (%)'. Some changes are as high as 5.287% and as low as -10%")
    st.dataframe(pr_diff(pr_df_sorted, topic_sensitive_pr_df_sorted))

if __name__ == '__main__':
    main()