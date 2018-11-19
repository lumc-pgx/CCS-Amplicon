#! /usr/bin/env python
from __future__ import print_function
import sys
import click
import pysam
import pandas as pd
import json


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="summarize phasing and clustering"
)
@click.argument("info_file", type=click.Path(exists=True))
@click.argument("cluster_file", type=click.Path(exists=True))
@click.argument("tagged_bams", type=click.Path(exists=True), nargs=-1)
def cli_handler(info_file, cluster_file, tagged_bams):
    try:
        info = pd.read_csv(info_file, sep="\t")
    except pd.errors.EmptyDataError:
        return

    with open(cluster_file, "r") as infile:
        clusters = json.load(infile)

    assigned_clusters = [-1 for _ in range(info.shape[0])]
    for cluster in clusters:
        for c in cluster["members"]:
            assigned_clusters[c] = cluster["cluster"]
    info["cluster"] = assigned_clusters

    assigned_phase = [-1 for _ in range(info.shape[0])]
    for bamfile in tagged_bams:
        with pysam.AlignmentFile(bamfile, "rb", check_sq=False) as infile:
            for record in infile.fetch(until_eof=True):
                idx = info[info["id"] == record.query_name].index.item()
                try:
                    phase = record.get_tag("HP")
                except KeyError:
                    phase = 0
                assigned_phase[idx] = phase
    info["phase"] = assigned_phase

    info.to_csv(sys.stderr, sep="\t", index=False)

if __name__ == '__main__':
    cli_handler()
