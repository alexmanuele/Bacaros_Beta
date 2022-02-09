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
        if i >= L:
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

# Delta S
def delta_S(plotA, plotB, L):
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

    delta_S = 1 - (TD / L)
    return delta_S

#delta T
def sample_to_frame(sample, L):
    records = []
    for species in sample:
        taxa = taxon_to_list(species, L)
        for i, t in enumerate(taxa):
            for check in ['unassigned', 'undefined']:
                if t != 0:
                    if check in t.lower():
                        taxa[i] = 0
        records.append({k:v for k, v in zip(range(len(taxa)), taxa)})

    return pd.DataFrame.from_records(records).replace(0, np.nan)
# This is the delta T metric
def delta_T(sample1, sample2, L):
    frame1 = sample_to_frame(sample1, L)
    frame2 = sample_to_frame(sample2, L)

    common = pd.concat([frame1, frame2]).drop_duplicates()
    shared_nodes = 0
    all_nodes = 0
    for col in common.columns:
        nodes = common[col].dropna().unique()
        f1_nodes = frame1[col].dropna().unique()
        f2_nodes= frame2[col].dropna().unique()
        for node in nodes:
            all_nodes += 1
            if node in f1_nodes and node in f2_nodes:
                shared_nodes += 1
    return shared_nodes / all_nodes

def calculate_beta(samples, L, metric):
    #The beta diversity is just an average of the delta_S metrics
    #We will also save a matrix of delta_S values.

    #Empty matrix to put values in
    matr = np.ones((len(samples), len(samples)))
    deltaMatrix = pd.DataFrame(data=matr,
                         index=[s['name'] for s in samples],
                         columns=[s['name'] for s in samples])
    distances = []
    records = []
    dfunc = {'s': delta_S,
             't': delta_T}[metric]
    #For each sample pair:
    pairs = it.combinations(samples, 2)
    for pair in pairs:
        #calculate delta_S for the pair.
        delta = dfunc(pair[0]['taxa'], pair[1]['taxa'], L)
        distances.append(delta)
        #Keep dicts to help populate the matrix.
        records.append({'pair0': pair[0]['name'], 'pair1': pair[1]['name'], 'delta': delta})

    for record in records:
        deltaMatrix.loc[record['pair0'], record['pair1']]= record['delta']
        deltaMatrix.loc[record['pair1'], record['pair0']]= record['delta']
    #Return the matrix and the average.
    return deltaMatrix.T, np.asarray(distances).mean()
