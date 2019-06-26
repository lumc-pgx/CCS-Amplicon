#! /usr/bin/env python
from __future__ import print_function, division
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform, euclidean
import markov_clustering as mc
import click
import json

def similarity_from_dists(dists, percentile):
    """
    Convert a distance matrix to a similarity matrix using Gaussian RBF.
    
    :param dists: distance matrix
    :type dists: numpy array-like
    :param percentile: Parameter used to tune the RBF.
    :type percentile: real number
    :return: the similarity matrix
    :rtype: numpy array 
    """
    # determine sigma based on the 'percentile' percentile of the distances
    sigma = np.percentile(dists, percentile)
    # RBF
    return np.exp(-dists ** 2 / (2. * sigma ** 2))


def distances_from_embeddings(embeddings):
    """
    Convert t-SNE embeddings into a euclidean distance matrix
    
    :param embeddings: embedded data
    :type embeddings: numpy array-like
    :return: distance matrix
    :rtype: numpy array
    """
    return squareform(pdist(embeddings))


def find_clusters(embeddings, similarity_threshold, inflation):
    """
    Identify clusters within embedded data.
    
    :param embeddings: embedded data
    :type embeddings: numpy array-like
    :param similarity_threshold: Can be used to tune the clustering
                                 performance.
    :type similarity_threshold: real number
    :param inflation: Markov clustering inflation. Used to control the
                      granularity of the clustering. Low values give 
                      fewer, larger clusters. Higher values give more,
                      smaller clusters.
    :type inflation: real number
    :return: The identified clusters
    :rtype: dict
    """
    dists = distances_from_embeddings(embeddings)
    similarity = similarity_from_dists(dists, similarity_threshold)
    results = mc.run_mcl(similarity, inflation=inflation)
    clusters = sorted(mc.get_clusters(results), key=lambda x: len(x), reverse=True)
    return clusters


def sort_cluster_elements(cluster, embeddings, info):
    # prioritize by distance from cluster center
    n_dim = len(embeddings[0])
    
    cluster_center = [
        np.mean([embeddings[c][i] for c in cluster for _ in range(info.iloc[c]["np"])])
        for i in range(n_dim)
    ]
    
    offset = [euclidean(cluster_center, embeddings[x]) for x in range(len(embeddings))]
    sorted_cluster = sorted(cluster, key=lambda x: (-offset[x], info.iloc[x]["rq"]), reverse=True)
    return sorted_cluster


def cluster_coverage(cluster, info):
    return int(sum([info.iloc[c]["np"] for c in cluster]))


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Identify clusters from a 2D embedding"
)
@click.option("--similarity-threshold", "-s", type=int, default=10,
              help="minimum similarity percentile")
@click.option("--inflation", "-i", type=float, default=1.4,
              help="MCL cluster inflation")
@click.argument("embeddings", type=click.Path(exists=True))
@click.argument("sequence_info", type=click.Path(exists=True))
def cli_handler(similarity_threshold, inflation, embeddings, sequence_info):
    try:
        embedded = pd.read_csv(embeddings, sep="\t", header=None).values
        info = pd.read_csv(sequence_info, sep="\t")
    except pd.errors.EmptyDataError:
        print(json.dumps([]))
        return

    clusters = find_clusters(embedded, similarity_threshold, inflation)
    clusters = [sort_cluster_elements(c, embedded, info) for c in clusters]

    print(
        json.dumps(
            [
                {
                    "cluster": i,
                    "members": cluster,
                    "coverage": cluster_coverage(cluster, info),
                    "molecules": len(cluster),
                }
                for i, cluster in enumerate(clusters) if len(cluster) > 0
            ],
            sort_keys=True, indent=4, separators=(",", ": ")
        )
    )


if __name__ == '__main__':
    cli_handler()
