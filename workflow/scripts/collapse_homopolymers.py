#! /usr/bin/env python
from __future__ import print_function
import sys
from Bio import SeqIO, Seq
import click


def collapse_homopolymers(sequence, max_length=2):
    assert max_length > 0, "max_length must be greater than 0"
    collapsed = []

    def base_counter():
        return {"nuc": "", "count": 0}
                   
    current = base_counter()
    for base in sequence:
        if base != current["nuc"]:
            collapsed.append(current["nuc"] * min(max_length, current["count"]))
            current = base_counter()                         

        current["nuc"] = base
        current["count"] += 1

    collapsed.append(current["nuc"] * min(max_length, current["count"]))

    return "".join(collapsed)


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Hompolymer repeat collapser"
)
@click.option("--size", "-s", type=int, default=2,
              help="maximum length of homopolymer repeats")
@click.argument("input_fasta", type=click.Path(exists=True))
def cli_handler(input_fasta, size):
    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            collapsed = collapse_homopolymers(record.seq, size)
            record.seq = Seq.Seq(collapsed)
            SeqIO.write(record, sys.stdout, "fasta")


if __name__ == '__main__':
    cli_handler()
