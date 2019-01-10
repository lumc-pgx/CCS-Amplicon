#! /usr/bin/env python
from __future__ import print_function
import click
import pysam
import edlib
from Bio.Seq import reverse_complement

START_IDX = 0
END_IDX = 1

F_PRIMER = 1
F_RC_PRIMER = -1

R_PRIMER = 2
R_RC_PRIMER = -2

NULL_PRIMER = 0
NULL_SCORE = -1

def find_and_remove_primers(bam_handle, f, r, max_dist):
    for record in bam_handle.fetch(until_eof=True):
        quals = record.query_qualities
        seq = record.query_sequence
        seq_rc = reverse_complement(seq)

        primers_found = [NULL_PRIMER, NULL_PRIMER]
        primer_scores = [NULL_SCORE, NULL_SCORE]

        f_aln = edlib.align(f, seq, mode="SHW")
        f_rc_aln = edlib.align(f, seq_rc, mode="SHW")

        f_score = f_aln["editDistance"]
        f_rc_score = f_rc_aln["editDistance"]

        r_aln = edlib.align(r, seq, mode="SHW")
        r_rc_aln = edlib.align(r, seq_rc, mode="SHW")

        r_score = r_aln["editDistance"]
        r_rc_score = r_rc_aln["editDistance"]

        start, end = (0, 0)

        if f_score <= max_dist:
            primers_found[START_IDX] = F_PRIMER
            primer_scores[START_IDX] = f_score
            start = f_aln["locations"][-1][-1] + 1

        if f_rc_score <= max_dist:
            primers_found[END_IDX] = F_RC_PRIMER
            primer_scores[END_IDX] = f_rc_score
            end = -f_rc_aln["locations"][-1][-1] - 1

        if r_score <= max_dist:
            primers_found[START_IDX] = R_PRIMER
            primer_scores[START_IDX] = r_score
            start = r_aln["locations"][-1][-1] + 1

        if r_rc_score <= max_dist:
            primers_found[END_IDX] = R_RC_PRIMER
            primer_scores[END_IDX] = r_rc_score
            end = -r_rc_aln["locations"][-1][-1] - 1

        if start != 0:
            seq = seq[start:]
            quals = quals[start:]

        if end != 0:
            seq = seq[:end]
            quals = quals[:end]

        record.query_sequence = seq
        record.query_qualities = quals
        record.tags += [("pf", primers_found)]
        record.tags += [("pd", primer_scores)]

        yield record


@click.command(
    context_settings=dict(
        ignore_unknown_options=False,
    ),
    short_help="Identify and remove primers from sequences"
)
@click.argument(
    "input_bamfile",
    type=click.Path(
        exists=True,
        dir_okay=False,
        readable=True
    )
)
@click.argument(
    "output_bamfile",
    type=click.Path(
        exists=False,
        dir_okay=False,
        writable=True
    )
)
@click.option(
    "--forward", "-f",
    required=True,
    help="forward primer sequence"
)
@click.option(
    "--reverse", "-r",
    required=True,
    help="reverse primer sequence"
)
@click.option(
    "--max-dist", "-d",
    type=click.IntRange(0, None),
    default=5,
    help="maximum edit distance for a primer hit"
)
def cli_handler(input_bamfile, output_bamfile, forward, reverse, max_dist):

    with pysam.AlignmentFile(input_bamfile, "rb", check_sq=False) as input_handle:
        with pysam.AlignmentFile(output_bamfile, "wb", template=input_handle) as output_handle:
            records = find_and_remove_primers(input_handle, forward, reverse, max_dist)

            for r in records:
                output_handle.write(r)


if __name__ == '__main__':
    cli_handler()
