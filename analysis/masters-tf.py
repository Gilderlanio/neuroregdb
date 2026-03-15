import pandas as pd
import networkx as nx
from itertools import combinations

df = pd.read_csv('regulons_TFs_masters.txt', sep='\t')

graphs = {}

for tissue, subdf in df.groupby("tissue"):
    G = nx.DiGraph()
    for _, row in subdf.iterrows():
        G.add_edge(
            row["tf"],
            row["target"],
            weight=row["mode"]
        )
    graphs[tissue] = G

regulons = {}

for tissue, G in graphs.items():
    regulons[tissue] = {}
    for tf in {u for u, v in G.edges()}:
        regulons[tissue][tf] = set(G.successors(tf))

regulon_size = []

for tissue, tf_dict in regulons.items():
    for tf, targets in tf_dict.items():
        regulon_size.append({
            "tf": tf,
            "tissue": tissue,
            "regulon_size": len(targets)
        })

regulon_size_df = pd.DataFrame(regulon_size)

results = []

tissues = list(regulons.keys())
tfs = set(df["tf"])

for tf in tfs:
    for t1, t2 in combinations(tissues, 2):

        targets_1 = regulons.get(t1, {}).get(tf, set())
        targets_2 = regulons.get(t2, {}).get(tf, set())

        size_1 = len(targets_1)
        size_2 = len(targets_2)

        if size_1 == 0 and size_2 == 0:
            continue

        intersection = targets_1 & targets_2
        union = targets_1 | targets_2

        intersection_size = len(intersection)
        union_size = len(union)

        # Metrics
        jaccard = intersection_size / union_size if union_size > 0 else 0
        overlap_coeff = (
            intersection_size / min(size_1, size_2)
            if min(size_1, size_2) > 0 else 0
        )

        containment_t1_in_t2 = intersection_size / size_1 if size_1 > 0 else 0
        containment_t2_in_t1 = intersection_size / size_2 if size_2 > 0 else 0

        results.append({
            "tf": tf,
            "tissue_1": t1,
            "tissue_2": t2,
            "regulon_size_tissue_1": size_1,
            "regulon_size_tissue_2": size_2,
            "overlap_size": intersection_size,
            "overlap_coefficient": overlap_coeff,
            "containment_t1_in_t2": containment_t1_in_t2,
            "containment_t2_in_t1": containment_t2_in_t1
        })

results_df = pd.DataFrame(results)

results_df = results_df[results_df["overlap_size"] > 0]

regulon_size_df.to_csv(
    'regulons_TFs_masters_regulon_size.csv',
    index=False
)

results_df.to_csv(
    'regulons_TFs_masters_overlap_metrics.csv',
    index=False
)

