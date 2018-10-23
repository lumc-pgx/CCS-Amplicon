#! /usr/bin/env python
from __future__ import print_function
import click
import pysam

@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="summarize bam tag information"
)
@click.argument("ccs_bamfile", type=click.Path(exists=True))
def cli_handler(ccs_bamfile):
    with pysam.AlignmentFile(ccs_bamfile, "rb", check_sq=False) as infile:
        print("\t".join(["id", "np", "rq", "len"]))
        for read in infile.fetch(until_eof=True):
            print(
                "\t".join([
                    str(x) for x in [ 
                        read.query_name,
                        read.get_tag("np"),
                        read.get_tag("rq"),
                        read.query_length
                    ]
                ])
            )


if __name__ == '__main__':
    cli_handler()
