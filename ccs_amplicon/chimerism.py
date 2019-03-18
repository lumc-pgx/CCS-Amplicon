#! /usr/bin/env python
from __future__ import print_function

from Bio import AlignIO
from Bio.Align import AlignInfo

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd
import subprocess
import os
import click


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Classify sequences which appear to be chimeric"
)
@click.argument("input_fasta", type=click.Path(exists=True))
@click.argument("phased_tsv", type=click.Path(exists=True))
def cli_handler(input_fasta, phased_tsv):
    # generate multiple sequence alignment from input sequences
    alignment = subprocess.run(
        ["mafft", "--auto", "--quiet", "--thread", "1", "--adjustdirection", input_fasta],
        stdout=subprocess.PIPE
    )

    # read the alignment
    try:
        msa = AlignIO.read(StringIO(alignment.stdout.decode(encoding="utf-8", errors="strict")), "fasta")
    except ValueError:
        msa = []

    if len(msa) < 3:
        return

    # generate alignment summary
    for alignment in msa:
        alignment.id = alignment.id.split("_R_")[-1]
        summary = AlignInfo.SummaryInfo(msa)

    # identify alignment portions which differ
    deltas = {
        s.id : { 
            "seq": [], # this will hold nucleotides
            "parent": [] # this will hold the hypothetical parent sequence for each nucleotide
        }
        for s in msa
    }

    # iterate over the columns of the msa
    for i in range(len(msa[0])):
        column = summary.get_column(i)
        # store the nucleotides for any columns which are not homogeneous
        if len([n for n in column if n == column[0]]) != len(column):
            for i, n in enumerate(column):
                deltas[msa[i].id]["seq"].append(n)

    # load the phasing info metadata
    info = pd.read_csv(phased_tsv, sep="\t")
    # group by cluster and phase then count sequences
    groups = info.groupby(["cluster", "phase"]).count()
        
    def molecule_count(seq_id, grouped_info):
        # get number of molecules for a given haplotype sequence
        fields = seq_id.split(".")
        cluster = int(fields[1].split("cluster")[-1])
        phase = int(fields[2].split("haplotype")[-1])
        return grouped_info.loc[(cluster, phase)]["np"]

    # list of sequence ids, sorted by molecule count
    ids = [s.id for s in msa]
    sorted_ids = sorted([s for s in deltas], key=lambda x: molecule_count(x, groups), reverse=True)

    def parent(pos, base):
        # determine the first sequence which supports the given base as the specified position
        for i, seq_id in enumerate(sorted_ids):
            parent_base = deltas[ids[ids.index(seq_id)]]["seq"][pos]
            if parent_base == base:
                return str(i)

        return "X" # no parent

    # determine parent sequence for each variant position
    for seq_id in deltas:
        deltas[seq_id]["parent"] = [parent(i,n) for i,n in enumerate(deltas[seq_id]["seq"])]
    
    # output the summary info
    seqid_size = max([len(a.id) for a in msa])
    
    format_string = "{{:{}}}\t{{:6d}}\t{{}}\t{{}}".format(seqid_size)
    
    for seq_id in sorted_ids:
        print(format_string.format(
            seq_id,
            molecule_count(seq_id, groups),
            "".join(deltas[seq_id]["seq"]),
            "".join(deltas[seq_id]["parent"])
        ))


if __name__ == '__main__':
    cli_handler()
