import numpy as np
import pandas as pd
from bct_for import fort as bf


def normalize_metrics_for(W, k, iteration=1):
    assert type(W) == pd.DataFrame
    M = W.copy().values.astype(float)
    Mf = np.asfortranarray(M)
    globa, noda = bf.normalize_metrics(Mf, k, iteration)
    graph_measures = ["degree centrality", "eigen centrality",
                      "mean strength", "clustering coefficient",
                      "network density", "rich club coefficient",
                      "rich club indiv", "betweenness centrality",
                      "global efficiency",
                      "characteristic path length", "radius",
                      "diameter"]
    node_measures = ["degree centrality", "eigen centrality",
                     "mean strength", "clustering coefficient",
                     "betweenness centrality"]
    norm_met = pd.DataFrame(globa, columns=["unnormalized", "normalized"],
                            index=graph_measures)
    node_stats = pd.DataFrame(noda.T, columns=node_measures,
                              index=W.index)
    return norm_met, node_stats
