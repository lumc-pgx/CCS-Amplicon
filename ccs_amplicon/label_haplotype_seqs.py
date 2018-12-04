#! /usr/bin/env python
from __future__ import print_function
import click
from Bio import SeqIO
import sys


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="rename the sequences created by LAA"
)
@click.option("--prefix", "-p", type=str, default="haplotype",
              help="prefix to use for the sequence id/name")
@click.argument("sequence_fastq", type=click.Path(exists=True))
def cli_handler(prefix, sequence_fastq):
    for record in SeqIO.parse(sequence_fastq, "fastq"):
        name = prefix
        record.id = name
        record.name = name
        record.description = name

        SeqIO.write(record, sys.stdout, "fastq")


if __name__ == '__main__':
    cli_handler()