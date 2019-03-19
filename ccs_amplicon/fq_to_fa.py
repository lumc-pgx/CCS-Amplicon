#! /usr/bin/env python
from __future__ import print_function
import sys
from Bio import SeqIO
import click


@click.command(context_settings=dict(
    ignore_unknown_options=True,
),
    short_help="convert fastq to fasta"
)
@click.argument("input_fastq", type=click.Path(exists=True))
def cli_handler(input_fastq):
    with open(input_fastq, "r") as input_handle:
        seqs = SeqIO.parse(input_fastq, "fastq")
        SeqIO.write(seqs, sys.stdout, "fasta")


if __name__ == '__main__':
    cli_handler()
