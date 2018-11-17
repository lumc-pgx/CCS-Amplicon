#! /usr/bin/env python
import pandas as pd
import json

from bokeh.plotting import figure, ColumnDataSource
from bokeh.palettes import viridis
from bokeh.transform import factor_cmap
from bokeh.io import output_file, save
from bokeh.models.widgets import Div

import click


def make_plot(clusters, filtered_clusters, embeddings, info, title):
    data = embeddings.copy()
    
    embedded_clusters = [-1 for _ in range(data.shape[0])]

    for cluster in clusters:
        for member in cluster["members"]:
            embedded_clusters[member] = cluster["cluster"]

    data["cluster"] = embedded_clusters

    in_filtered = [0 for _ in range(data.shape[0])]
    for cluster in filtered_clusters:
        for member in cluster["members"]:
            in_filtered[member] = 1


    datasource1 = ColumnDataSource({
        "x": data[0], 
        "y": data[1],
        "passes": info["np"],
        "len": info["len"],
        "qual": info["rq"],
        "cluster": [str(c) for c in data["cluster"]],
        "filter": in_filtered
    })

    TOOLTIPS = [
        ("Passes", "@passes"),
        ("length", "@len"),
        ("qual", "@qual"),
        ("index", "$index")
    ]

    factors = [str(i) for i in set(data["cluster"])]

    p=figure(title=title, tooltips=TOOLTIPS)
    c = p.circle(x="x", y="y", source=datasource1, size="passes",
                 color=factor_cmap('cluster', palette=viridis(max(2, len(factors))), factors=factors),
                 alpha=0.5, line_alpha="filter", legend="cluster")

    return p


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Generate an html plot of the clustering"    
)
@click.option("--title", "-t", type=str, default="Clusters",
              help="Title for the plot")
@click.argument("cluster_file", type=click.Path(exists=True))
@click.argument("filtered_cluster_file", type=click.Path(exists=True))
@click.argument("embedding_file", type=click.Path(exists=True))
@click.argument("info_file", type=click.Path(exists=True))
@click.argument("plot_file", type=click.Path())
def cli_handler(title, cluster_file, filtered_cluster_file, embedding_file, info_file, plot_file):
    with open(cluster_file, "r") as infile:
        clusters = json.load(infile)

    with open(filtered_cluster_file, "r") as infile:
        filtered_clusters = json.load(infile)

    try:
        embeddings = pd.read_csv(embedding_file, sep="\t", header=None)
        bam_info = pd.read_csv(info_file, sep="\t")
    except pd.errors.EmptyDataError:
        pass

    output_file(plot_file, title="Cluster report", mode="inline")

    if len(clusters) > 0:
        figure = make_plot(clusters, filtered_clusters, embeddings, bam_info, title)
    else:
        figure = Div(text="No Data")

    save(figure) 


if __name__ == '__main__':
    cli_handler()
