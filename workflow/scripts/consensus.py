#! /usr/bin/env python
from __future__ import print_function
import click
from Bio import AlignIO
from Bio.Align import AlignInfo


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="generate consensus sequence from MSA"
)
@click.option("--fraction", "-f", type=float, default=0.51,
              help="call consensus for bases with frequency >= fraction at a position")
@click.option("--seqid", "-s", type=str, default="consensus",
              help="fasta sequence id to use for consensus sequence")
@click.argument("msa_fasta", type=click.Path(exists=True))
def cli_handler(fraction, seqid, msa_fasta):
    try:
        msa = AlignIO.read(msa_fasta, "fasta")
    except ValueError:
        return

    summary = AlignInfo.SummaryInfo(msa)
    consensus = str(summary.gap_consensus(fraction, require_multiple=True, ambiguous='')).replace("-", "")

    print(">{}".format(seqid))
    print(consensus)

if __name__ == '__main__':
    cli_handler()