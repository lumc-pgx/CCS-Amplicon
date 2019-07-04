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

# degenerate nucleotides
IUPAC_EQUALITIES = [
    ("R", "A"), ("R", "G"),
    ("Y", "C"), ("Y", "T"),
    ("M", "A"), ("M", "C"),
    ("K", "G"), ("K", "T"),
    ("S", "C"), ("S", "G"),
    ("W", "A"), ("W", "T"),
    ("H", "A"), ("H", "C"), ("H", "T"),
    ("B", "C"), ("B", "G"), ("B", "T"),
    ("V", "A"), ("V", "C"), ("V", "G"),
    ("D", "A"), ("D", "G"), ("D", "T"),
    ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")
]


def align(query, target):
    """
    Align query sequence to target sequence using edlib in prefix mode.
    IUPAC codes for degenerate nucleotides are supported.
    :param query: The query sequence
    :type query: string
    :param target: The target sequence
    :type target: string
    :return: The alignment
    :rtype: edlib alignment dictionary
    """
    return edlib.align(query, target, mode="SHW", additionalEqualities=IUPAC_EQUALITIES)


def find_primers(f, r, target, max_dist):
    seq = target.upper()
    seq_rc = reverse_complement(seq)

    # initially, no primers found
    primers_found = [NULL_PRIMER, NULL_PRIMER]
    primer_scores = [NULL_SCORE, NULL_SCORE]

    # start and end positions for trimmed sequence
    start, end = (0, 0)

    if f:
        # align and score
        f_aln = align(f, seq)        # forward vs seq
        f_rc_aln = align(f, seq_rc)  # forward vs revcomp(seq)

        f_score = f_aln["editDistance"]
        f_rc_score = f_rc_aln["editDistance"]

        # forward primer match
        if f_score <= max_dist:
            primers_found[START_IDX] = F_PRIMER
            primer_scores[START_IDX] = f_score
            start = f_aln["locations"][-1][-1] + 1

        # forward primer revcomp match
        if f_rc_score <= max_dist:
            primers_found[END_IDX] = F_RC_PRIMER
            primer_scores[END_IDX] = f_rc_score
            end = -f_rc_aln["locations"][-1][-1] - 1

    if r:
        # align and score
        r_aln = align(r, seq)        # rev vs seq
        r_rc_aln = align(r, seq_rc)  # rev vs revcomp(seq)

        r_score = r_aln["editDistance"]
        r_rc_score = r_rc_aln["editDistance"]

        # reverse primer match
        if r_score <= max_dist:
            primers_found[START_IDX] = R_PRIMER
            primer_scores[START_IDX] = r_score
            start = r_aln["locations"][-1][-1] + 1

        # reverse primer revcomp match
        if r_rc_score <= max_dist:
            primers_found[END_IDX] = R_RC_PRIMER
            primer_scores[END_IDX] = r_rc_score
            end = -r_rc_aln["locations"][-1][-1] - 1

    return primers_found, primer_scores, start, end


def find_and_remove_primers(bam_handle, f, r, max_dist):
    f = f.upper()
    r = r.upper()

    # each bam record
    for record in bam_handle.fetch(until_eof=True):
        # qualities and sequence from bam record
        quals = record.query_qualities
        seq = record.query_sequence

        primers_found, primer_scores, start, end = find_primers(f, r, seq, max_dist)

        # remove primer from start of sequence.
        trimmed = record.query_sequence
        if start != 0:
            trimmed = trimmed[start:]
            quals = quals[start:]

        # remove primer from end of sequence
        if end != 0:
            trimmed = trimmed[:end]
            quals = quals[:end]

        # update bam record
        record.query_sequence = trimmed
        record.query_qualities = quals
        record.tags += [("pf", primers_found)]
        record.tags += [("pd", primer_scores)]
        record.tags += [("pp", (start, end))]

        yield record


@click.command(context_settings=dict(
    help_option_names=["-h", "--help"]
    ),
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
    type=click.STRING,
    default="",
    help="forward primer sequence"
)
@click.option(
    "--reverse", "-r",
    type=click.STRING,
    default="",
    help="reverse primer sequence"
)
@click.option(
    "--max-dist", "-d",
    type=click.IntRange(0, None),
    default=5,
    help="maximum edit distance for a primer hit"
)
def cli_handler(input_bamfile, output_bamfile, forward, reverse, max_dist):
    """Identify and remove primers from sequences"""
    with pysam.AlignmentFile(input_bamfile, "rb",
                             check_sq=False) as input_handle:
        with pysam.AlignmentFile(output_bamfile, "wb",
                                 template=input_handle) as output_handle:
            records = find_and_remove_primers(input_handle, forward,
                                              reverse, max_dist)
            for r in records:
                output_handle.write(r)


if __name__ == '__main__':
    cli_handler()
