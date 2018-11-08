#! /usr/bin/env python
from __future__ import print_function
import sys
from sklearn.manifold import TSNE
import pandas as pd
import click


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Create 2D tSNE embedding from distance matrix"
)
@click.option("--iterations", "-i", type=int, default=5000,
              help="maximum number of tSNE iterations")
@click.option("--learning_rate", "-lr", type=int, default=50,
              help="tSNE learning rate")
@click.option("--seed", "-s", type=int, default=42,
              help="random seed to ensure reproducability")
@click.argument("distance_matrix", type=click.Path(exists=True))
def cli_handler(iterations, learning_rate, seed, distance_matrix):
    distances = pd.read_csv(distance_matrix, sep="\t", header=None)
    tsne=TSNE(n_components=2, metric="precomputed", 
          n_iter=iterations, learning_rate=learning_rate, 
          method="exact", verbose=2, random_state=seed)

    old_stdout = sys.stdout
    sys.stdout = sys.stderr
    embedded = tsne.fit_transform(distances)

    sys.stdout = old_stdout
    for e in embedded:
        print("{}\t{}".format(e[0], e[1]))


if __name__ == '__main__':
    cli_handler()
