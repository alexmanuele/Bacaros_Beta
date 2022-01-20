import numpy as np
import pandas as pd
import itertools as it

def qiime2_tsv_to_taxa_list(file):
    df = pd.read_table(file)
    taxa = df[['Taxon']].drop_duplicates()['Taxon']
    return taxa
    
# For each taxon, create a list of the rank values.
# Unassigned ranks will be left as "0"
# The cardinality of the list is given as L. E.g. Species, L=7.
def taxon_to_list(taxon, L):
    assert L > 0
    assert type(L) == int
    n = np.zeros(L, dtype='object')
    for i, rank in enumerate(taxon.split(";")):
        if i > L:
            break
        n[i] = rank
    return n

# Compute the distance as described above
def compare_taxa_lists(taxon1, taxon2):
    assert taxon1.shape == taxon2.shape
    distance = 0
    i = taxon1.shape[0] - 1
    while i >= 0:
        # If unassigned no match.
        if taxon1[i] == 0:
            distance += 1
        # If they arent zero and they match, stop
        elif taxon1[i] == taxon2[i]:
            break
        #Otherwise, there's no match.
        else:
            distance += 1
        i -= 1
    return distance

def compare_plots(plotA, plotB, L):
    #Make an empty distance matrix with taxaA rows and taxaB columns
    matr = np.zeros((len(plotA), len(plotB)))
    dist = pd.DataFrame(data=matr, index=plotA, columns=plotB)
    #calculate pairwise distances
    pairs = it.product(plotA, plotB)
    for pair in pairs:
        #pairwise distance
        w_ij = compare_taxa_lists(taxon_to_list(pair[0], L), taxon_to_list(pair[1], L))
        dist.loc[pair[0], pair[1]] = w_ij
    row_minima = np.zeros(dist.shape[0])
    col_minima = np.zeros(dist.shape[1])

    for i in range(dist.shape[0]):
        row_minima[i] = dist.iloc[i].min()

    for i, col in enumerate(dist.columns):
        col_minima[i] = dist[col].min()

    TD = (row_minima.sum() + col_minima.sum()) / (dist.shape[0] + dist.shape[1])

    delta_S = 1 - (TD / (L-1))
    return delta_S

def calculate_beta(samples, L):
    distances = []
    pairs = it.combinations(samples, 2)
    for pair in pairs:
        delta_S = compare_plots(pair[0], pair[1], L)
        distances.append(delta_S)
    return np.asarray(distances).mean()
