#! /usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
import numpy as np
import edlib
import click
from tqdm import tqdm

def get_direction(query, ref):
    """
    Determine the strand direction for query sequence by comparison with
    ref sequence
    """
    fwd = query
    rev = query.reverse_complement()
    d_fwd = edlib.align(str(fwd.seq), str(ref.seq), task="distance", mode="NW")["editDistance"]
    d_rev = edlib.align(str(rev.seq), str(ref.seq), task="distance", mode="NW")["editDistance"]
    return fwd if d_fwd < d_rev else rev


def distance_matrix(sequences):
    """
    Construct a distance matrix from pairwise alignments sequences
    """
    dists = np.array([np.array([0 for _ in range(len(sequences))]) for _ in range(len(sequences))])
    
    if dists.shape[0] == 0:
        return dists

    base_seq = sequences[0]
    adjusted_sequences = []
    for s in tqdm(sequences, desc="{:<10}".format("prescan")):
        adjusted_sequences.append(get_direction(s, base_seq))
    
    for i in tqdm(range(len(dists)), desc="{:<10}".format("align")):
        query = str(adjusted_sequences[i].seq)
    
        for j in range(i, len(dists[i])):
            if i != j:
                target = str(adjusted_sequences[j].seq)
                d = edlib.align(query, target, task="distance", mode="NW")["editDistance"]
                dists[i][j] = d
                dists[j][i] = d
    return dists


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Create a distance matrix for sequences in a fasta file"
)
@click.argument("input_fasta", type=click.Path(exists=True))
def cli_handler(input_fasta):
    with open(input_fasta, "r") as infile:
        sequences = [record for record in SeqIO.parse(infile, "fasta")]
   
    distances = distance_matrix(sequences)
    for d in distances:
        print("\t".join([str(x) for x in d]))


if __name__ == '__main__':
    cli_handler()
