#! /usr/bin/env python
from __future__ import print_function
import click
import pysam
import sys


@click.command(context_settings=dict(
    ignore_unknown_options=True,
    ),
    short_help="Ensure that single heterozygous variants are assigned a phase"
)
@click.argument("vcf_file", type=click.Path(exists=True))
def cli_handler(vcf_file):
    with pysam.VariantFile(vcf_file) as vcf_in:
        header = vcf_in.header
        records = list(r for r in vcf_in.fetch())

        if len(records) == 1:
            record = records[0]
            if not "PS" in header.formats:
                header.formats.add("PS", 1, "Integer", "Phase set identifier")

            sample = record.samples["sm"]
            genotype = sample["GT"]

            if not sample.phased and genotype[0] != genotype[1]: # heterozygous
                sample["PS"] = record.pos
                sample.phased = True

        with pysam.VariantFile(sys.stdout, "w", header=header) as vcf_out:
            for record in records:
                vcf_out.write(record)

if __name__ == '__main__':
    cli_handler()