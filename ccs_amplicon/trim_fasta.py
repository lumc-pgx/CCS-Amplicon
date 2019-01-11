#! /usr/bin/env python
from __future__ import print_function
import sys
from Bio import SeqIO
import click


def trim_ends(record, n):
    if n == 0:
        return record
    try:
        return record[n:-n]
    except IndexError:
        return None


def trim_seqs(seqs, n):
    for record in seqs:
        trimmed = trim_ends(record, n)
        if trimmed is not None:
            yield trimmed


@click.command(context_settings=dict(
    ignore_unknown_options=True,
),
    short_help="Trim nucleotides from start and end of sequences"
)
@click.option("--num-nucleotides", "-n", type=click.IntRange(0, None), default=20,
              help="number of nucleotides to remove from ends of sequence. "
                   "A value of 0 leaves the original sequence unaffected.")
@click.argument("input_fasta", type=click.Path(exists=True))
def cli_handler(num_nucleotides, input_fasta):
    seqs = SeqIO.parse(input_fasta, "fasta")
    trimmed = trim_seqs(seqs, num_nucleotides)
    SeqIO.write(trimmed, sys.stdout, "fasta")


if __name__ == '__main__':
    cli_handler()
