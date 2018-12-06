#! /usr/bin/env python
from __future__ import print_function
import pandas as pd
import click
import json


def cluster_filter(clusters, threshold, max_size, info):
    sorted_clusters = sorted(clusters, key=lambda x: len(x["members"]), reverse=True)
    n_largest = len(sorted_clusters[0]["members"])
    cutoff = threshold * n_largest

    filtered_clusters = []
    for cluster in sorted_clusters:
        if len(cluster["members"]) >= cutoff:
            cluster["members"] = cluster["members"][:min(len(cluster["members"]), max_size)]
            cluster["coverage"] = int(sum([info.iloc[c]["np"] for c in cluster["members"]]))
            cluster["molecules"] = len(cluster["members"])
            filtered_clusters.append(cluster)

    return filtered_clusters


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Filter clusters"
)
@click.option("--inclusion-threshold", "-i", type=click.FloatRange(0.0, 1.0), default=0.2,
              help="fraction of molecules relative to the largest cluster required for a cluster to be included")
@click.option("--max-cluster-seqs", "-m", type=int, default=200,
              help="maximum number of CCS sequences to retain per-cluster")
@click.argument("clusters_json", type=click.Path(exists=True))
@click.argument("sequence_info", type=click.Path(exists=True))
def cli_handler(inclusion_threshold, max_cluster_seqs, clusters_json, sequence_info):
    with open(clusters_json, "r") as clusterfile:
        clusters = json.load(clusterfile)

    try:
        info = pd.read_csv(sequence_info, sep="\t")
    except pd.errors.EmptyDataError:
        assert len(clusters) == 0

    if len(clusters) > 0:
        filtered = cluster_filter(clusters, inclusion_threshold, max_cluster_seqs, info)
    else:
        filtered = clusters

    print(json.dumps(filtered, sort_keys=True, indent=4, separators=(",", ": ")))



if __name__ == '__main__':
    cli_handler()
