#! /usr/bin/env python
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
    for cluster in clusters:
        print(" ".join([str(c) for c in cluster]))


if __name__ == '__main__':
    cli_handler()
