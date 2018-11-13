#! /usr/bin/env python
from __future__ import print_function
import click
import pysam


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Generate haplotype whitelists"
)
@click.option("--prefix", "-p", type=str, default="",
              help="prefix for whitelist file names")
@click.option("--min-molecules", "-m", type=int, default=10,
              help="minimum number of molecules required for a haplotype")
@click.argument("bam_file", type=click.Path(exists=True))
def cli_handler(prefix, min_molecules, bam_file):
    tagged_reads = {}
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as infile:
        for r in infile.fetch(until_eof=True):
            qname = r.query_name
            molecule = "/".join(qname.split("/")[:-1])
            try:
                haplotype = r.get_tag("HP")
            except KeyError:
                haplotype = 0

            if haplotype not in tagged_reads:
                tagged_reads[haplotype] = []

            tagged_reads[haplotype].append(molecule)

    if prefix != "":
        prefix += "."

    tagged_reads = {haplotype: molecules for haplotype, molecules in tagged_reads.items() if len(molecules) >= min_molecules}

    if len(tagged_reads) == 0:
        with open("{}haplotype0.whitelist".format(prefix), "w") as whitelist:
            pass
    else:
        for haplotype in tagged_reads:
            with open("{}haplotype{}.whitelist".format(prefix, haplotype), "w") as whitelist:
                for molecule in tagged_reads[haplotype]:
                    print(molecule, file=whitelist)

if __name__ == '__main__':
    cli_handler()