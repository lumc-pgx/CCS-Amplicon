#! /usr/bin/env python
from __future__ import print_function, division
import pandas as pd
import numpy as np
import markov_clustering as mc
import click

def similarity_from_dists(dists):
    max_dist = np.max(dists)
    return 1.0 - dists / max_dist


def euclidean_dist_squared(p1, p2):
    return np.power(p1[0] - p2[0], 2) + np.power(p1[1] - p2[1], 2)


def euclidean_dist(p1, p2):
    return np.sqrt(euclidean_dist_squared(p1, p2))


def distances_from_embeddings(embeddings):
    dists = []
    for p1 in embeddings:
        row = []
        for p2 in embeddings:
            dist = euclidean_dist(p1, p2)
            row.append(dist)
        dists.append(row)
    return dists


def find_clusters(embeddings, similarity_threshold, inflation):
    dists = distances_from_embeddings(embeddings)
    similarity = similarity_from_dists(dists)
    similarity[similarity < np.percentile(similarity, similarity_threshold)] = 0
    results = mc.run_mcl(similarity, inflation=1.4)
    clusters = sorted(mc.get_clusters(results), key=lambda x: len(x), reverse=True)
    return clusters


def sort_cluster_elements(cluster, embeddings):
    # prioritize by distance from cluster center
    cluster_center = [np.mean([embeddings[c][0] for c in cluster]), np.mean([embeddings[c][1] for c in cluster])]
    print(cluster_center)
    offset = [euclidean_dist(cluster_center, embeddings[x]) for x in range(len(embeddings))]
    
    sorted_cluster = sorted(cluster, key=lambda x: offset[x])
    cluster_offsets = [offset[c] for c in sorted_cluster]
    print(cluster_offsets)
    return sorted_cluster

@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Identify clusters from a 2D embedding"
)
@click.option("--similarity_threshold", "-s", type=int, default=80,
              help="minimum similarity percentile")
@click.option("--inflation", "-i", type=float, default=1.4,
              help="MCL cluster inflation")
@click.argument("embeddings", type=click.Path(exists=True))
def cli_handler(similarity_threshold, inflation, embeddings):
    embedded = pd.read_csv(embeddings, sep="\t", header=None).values
    clusters = find_clusters(embedded, similarity_threshold, inflation)
    clusters = [sort_cluster_elements(c, embedded) for c in clusters]
    for cluster in clusters:
        print(" ".join([str(c) for c in cluster]))


if __name__ == '__main__':
    cli_handler()
