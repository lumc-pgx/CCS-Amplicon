#! /usr/bin/env python
from __future__ import print_function, division
import pandas as pd
import numpy as np
import markov_clustering as mc
import click
import json

def similarity_from_dists(dists):
    max_dist = np.max(dists)
    return 1.0 - dists / max_dist


def euclidean_dist_squared(v1, v2):
    return np.sum([np.power(v1[i] - v2[i], 2) for i in range(len(v1))])


def euclidean_dist(v1, v2):
    return np.sqrt(euclidean_dist_squared(v1, v2))


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


def sort_cluster_elements(cluster, embeddings, info):
    # prioritize by distance from cluster center
    n_dim = len(embeddings[0])
    
    cluster_center = [
        np.mean([embeddings[c][i] for c in cluster for _ in range(info.iloc[c]["np"])])
        for i in range(n_dim)
    ]
    
    offset = [euclidean_dist(cluster_center, embeddings[x]) for x in range(len(embeddings))]
    sorted_cluster = sorted(cluster, key=lambda x: (-offset[x], info.iloc[x]["rq"]), reverse=True)
    return sorted_cluster


def cluster_coverage(cluster, info):
    return int(sum([info.iloc[c]["np"] for c in cluster]))


def get_coverage(clusters, info):
    return sum([sum([cluster_coverage(cluster, info)]) for cluster in clusters])


def filter_cluster_coverage(clusters, info, max_coverage):
    total_coverage = get_coverage(clusters, info)
    while total_coverage > max_coverage:
        clusters = [cluster[:-1] for cluster in clusters]
        total_coverage = get_coverage(clusters, info)
    return clusters


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Identify clusters from a 2D embedding"
)
@click.option("--similarity_threshold", "-s", type=int, default=80,
              help="minimum similarity percentile")
@click.option("--inflation", "-i", type=float, default=1.4,
              help="MCL cluster inflation")
@click.option("--coverage", "-c", type=int, default=500,
              help="limit the total coverage")
@click.argument("embeddings", type=click.Path(exists=True))
@click.argument("sequence_info", type=click.Path(exists=True))
def cli_handler(similarity_threshold, inflation, coverage, embeddings, sequence_info):
    embedded = pd.read_csv(embeddings, sep="\t", header=None).values
    info = pd.read_csv(sequence_info, sep="\t")
    clusters = find_clusters(embedded, similarity_threshold, inflation)
    clusters = [sort_cluster_elements(c, embedded, info) for c in clusters]
    clusters = filter_cluster_coverage(clusters, info, coverage)

    print(
        json.dumps(
            [
                {
                    "cluster": i,
                    "members": cluster,
                    "coverage": cluster_coverage(cluster, info),
                }
                for i, cluster in enumerate(clusters) if len(cluster) > 0
            ],
            sort_keys=True, indent=4, separators=(",", ": ")
        )
    )


if __name__ == '__main__':
    cli_handler()