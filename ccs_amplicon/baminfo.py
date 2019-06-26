#! /usr/bin/env python
from __future__ import print_function
import click
import pysam

def record_summary(record):
    """
    Summarize a single record
    
    :param record: The record to summarize
    :type record: A pysam.AlignedSegment object
    :return: Summary of specific fields from given record
    :rtype: dict
    """
    return {
        "id": record.query_name,
        "np": record.get_tag("np"),
        "rq": record.get_tag("rq"),
        "len": record.query_length
    }


def summarize(ccs_bam_path):
    """
    Summarize a bam file of CCS sequences.
    Summary information is written in tabular format to stdout.
    
    :param ccs_bam_path: Path to the ccs bam file
    :type ccs_bam_path: string
    """
    fields = ["id", "np", "rq", "len"]
    with pysam.AlignmentFile(ccs_bam_path, "rb", check_sq=False) as infile:
        print("\t".join(fields))
        for read in infile.fetch(until_eof=True):
            summary = record_summary(read)
            values = [str(summary[f]) for f in fields]
            print("\t".join(values))


@click.command(context_settings=dict(
    help_option_names=["-h", "--help"]
    ),
)
@click.argument("ccs_bamfile", type=click.Path(exists=True))
def cli_handler(ccs_bamfile):
    """
    Generate tabular summary information from a bam file containing
    pacbio CCS sequences
    """
    summarize(ccs_bamfile)


if __name__ == '__main__':
    cli_handler()
