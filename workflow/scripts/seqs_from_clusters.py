#! /usr/bin/env python
from __future__ import print_function
import click
from Bio import SeqIO
import json


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="generate fasta files for clusters"
)
@click.option("--prefix", "-p", type=str, default="",
              help="prefix for output file names")
@click.argument("sequence_fasta", type=click.Path(exists=True))
@click.argument("clusters_json", type=click.Path(exists=True))
def cli_handler(prefix, sequence_fasta, clusters_json):
    with open(clusters_json, "r") as infile:
        clusters = json.load(infile)

    if prefix != "":
        prefix = prefix + "."

    output_files = [open("{}cluster{}.fasta".format(prefix, cluster["cluster"]), "w") for cluster in clusters]

    for i, record in enumerate(SeqIO.parse(sequence_fasta, "fasta")):
        # output sequence to cluster sequence file
        for j, cluster in enumerate(clusters):
            if i in cluster["members"]:
                SeqIO.write(record, output_files[j], "fasta")
                break

    # cleanup
    for f in output_files:
        try:
            f.close()
        except:
            pass

if __name__ == '__main__':
    cli_handler()
