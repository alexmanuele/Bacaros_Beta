import pandas as pd

def qiime2_tsv_to_taxa_list(file):
    df = pd.read_table(file)
    taxa = df[['Taxon']].drop_duplicates()['Taxon']
    return taxa
