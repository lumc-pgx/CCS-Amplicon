#! /usr/bin/env python
from __future__ import print_function
import click
import pysam

def record_summary(record):
    return {
        "id": record.query_name,
        "np": record.get_tag("np"),
        "rq": record.get_tag("rq"),
        "len": record.query_length
    }


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="summarize bam tag information"
)
@click.argument("ccs_bamfile", type=click.Path(exists=True))
def cli_handler(ccs_bamfile):
    fields = ["id", "np", "rq", "len"]
    with pysam.AlignmentFile(ccs_bamfile, "rb", check_sq=False) as infile:
        print("\t".join(fields))
        for read in infile.fetch(until_eof=True):
            summary = record_summary(read)
            values = [str(summary[f]) for f in fields]
            print("\t".join(values))


if __name__ == '__main__':
    cli_handler()
