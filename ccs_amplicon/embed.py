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
@click.option("--learning-rate", "-lr", type=int, default=50,
              help="tSNE learning rate")
@click.option("--method", "-m", type=click.Choice(["exact", "barnes_hut"]), default="exact",
              help="gradient calculation method. Can be one of ('exact', 'barnes_hut')")
@click.option("--seed", "-s", type=int, default=42,
              help="random seed to ensure reproducability")
@click.option("--dimensions", "-d", type=int, default=2,
              help="Number of dimensions for embedding")
@click.argument("distance_matrix", type=click.Path(exists=True))
def cli_handler(iterations, learning_rate, method, seed, dimensions, distance_matrix):
    try:
        distances = pd.read_csv(distance_matrix, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        sys.stderr.write("Empty distance matrix")
        return
    
    tsne=TSNE(n_components=dimensions, metric="precomputed",
          n_iter=iterations, learning_rate=learning_rate, 
          method=method, verbose=2, random_state=seed)

    old_stdout = sys.stdout
    sys.stdout = sys.stderr
    embedded = tsne.fit_transform(distances)

    sys.stdout = old_stdout
    for e in embedded:
        print("\t".join([str(i) for i in e]))


if __name__ == '__main__':
    cli_handler()
