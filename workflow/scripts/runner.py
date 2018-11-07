#! /usr/bin/env python
import click
from snakemake import snakemake

@click.command(context_settings=dict(
    ignore_unknown_options=False,
    ),
    short_help="CCS driven amplicon polishing and phasing for PacBio sequencing data"
)
@click.option("--directory", "-d", type=click.Path(file_okay=False), default=".",
              help="path to the directory which the program output will be written to."
                   "If directory is omitted, output will be written to the current directory")
@click.option("--min-ccs-length", type=click.IntRange(1, None), default=6000,
              help="minimum required CCS sequence length")
@click.option("--max-ccs-length", type=click.IntRange(1, None), default=7000,
              help="maximum allowed CCS sequence length")
@click.option("--min-ccs-passes", type=click.IntRange(1, None), default=1,
              help="minimum required CCS passes")
@click.option("--min-ccs-qual", type=click.FloatRange(0.0, 1.0), default=0.98,
              help="minimum required CCS sequence quality")
@click.option("--max-homopolymer", type=click.IntRange(1, None), default=2,
              help="homopolymer runs of length greater than max-homopolymer will be collapsed to this length")
@click.option("--tsne-iterations", type=click.IntRange(1, None), default=5000,
              help="number of iterations for tSNE")
@click.option("--tsne-rate", type=click.IntRange(1, None), default=50,
              help="tSNE learning rate")
@click.option("--cluster-percentile", type=click.FloatRange(0.0, 1.0), default=0.8,
              help="cluster members with a similarity above this percentile are considered to be connected."
                   "lower values result in fewer, larger clusters.")
@click.option("--cluster-inflation", type=click.FloatRange(0, None), default=1.4,
              help="markov clustering inflation parameter."
                   "Higer values result in more clusters")
@click.option("--max-cluster-passes", type=click.IntRange(1, None), default=500,
              help="Limit the number of sequences output by the cluster step")
@click.option("--consensus-fraction", type=click.FloatRange(0.0, 1.0), default=0.51,
              help="frequency of nucleotide at a given position required for rough consensus calling")
@click.option("--min-haplotype-molecules", type=click.IntRange(1, None), default=10,
              help="minimum number of molecules (CCS sequences) required for a haplotype")
@click.option("--min-variant-qual", type=click.IntRange(0, None), default=50,
              help="minimum variant qual score required for a variant to be used for phasing")
@click.argument("ccs_bam", type=click.Path(exists=True))
@click.argument("subreads_bam", type=click.Path(exists=True))
def cli_handler(directory, min_ccs_length, max_ccs_length, min_ccs_passes, min_ccs_qual, max_homopolymer,
                tsne_iterations, tsne_rate, cluster_percentile, cluster_inflation, max_cluster_passes,
                consensus_fraction, min_haplotype_molecules, min_variant_qual,
                ccs_bamfile, subreads_bamfile,):
    pass


if __name__ == "__main__":
    cli_handler()