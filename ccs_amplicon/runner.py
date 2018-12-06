#! /usr/bin/env python
import os
import subprocess
import click
from snakemake import snakemake

@click.command(context_settings=dict(
    ignore_unknown_options=False,
    ),
    short_help="CCS-driven amplicon phasing and polishing for targetted PacBio sequencing data"
)
@click.option("--directory", "-d", type=click.Path(file_okay=False), default=".",
              help="Path to the directory which the program output will be written to."
                   "If directory is omitted, output will be written to the current directory")
@click.option("--prefix", "-p", type=str, default="ccs_amplicon",
              help="Prefix to use for output file names")
@click.option("--profile", "-po", type=str, default="",
              help="The name of the snakemake profile to use when running the workflow."
                   "Use this profile to control how the workflow is run on your specific compute architecture."
                   "See https://snakemake.readthedocs.io/en/stable/executable.html?highlight=profile for details")
@click.option("--min-ccs-length", type=click.IntRange(1, None), default=6000,
              help="Minimum required CCS sequence length")
@click.option("--max-ccs-length", type=click.IntRange(1, None), default=7000,
              help="Maximum allowed CCS sequence length")
@click.option("--min-ccs-passes", type=click.IntRange(1, None), default=1,
              help="Minimum required CCS passes")
@click.option("--min-ccs-qual", type=click.FloatRange(0.0, 1.0), default=0.98,
              help="Minimum required CCS sequence quality")
@click.option("--max-homopolymer", type=click.IntRange(1, None), default=2,
              help="Homopolymer runs of length greater than max-homopolymer will be collapsed to this length")
@click.option("--trim-ends", type=int, default=5,
              help="remove n bases from start and ends of sequences before clustering")
@click.option("--tsne-iterations", type=click.IntRange(1, None), default=5000,
              help="Number of iterations for tSNE")
@click.option("--tsne-rate", type=click.IntRange(1, None), default=50,
              help="tSNE learning rate")
@click.option("--cluster-percentile", type=click.IntRange(0, 100), default=80,
              help="Cluster members with a similarity above this percentile are considered to be connected."
                   "Lower values result in fewer, larger clusters.")
@click.option("--cluster-inflation", type=click.FloatRange(0, None), default=1.4,
              help="Markov clustering inflation parameter."
                   "Higer values result in more clusters")
@click.option("--cluster-size-threshold", type=click.FloatRange(0.0, 1.0), default=0.2,
              help="fraction of molecules relative to the largest cluster required for a cluster to be included")
@click.option("--max-cluster-size", type=int, default=200,
              help="maximum number of molecules to retain per-cluster")
@click.option("--consensus-fraction", type=click.FloatRange(0.0, 1.0), default=0.51,
              help="Frequency of nucleotide at a given position required for rough consensus calling")
@click.option("--min-haplotype-molecules", type=click.IntRange(1, None), default=10,
              help="Minimum number of molecules (CCS sequences) required for a haplotype")
@click.option("--min-variant-qual", type=click.IntRange(0, None), default=50,
              help="Minimum variant qual score required for a variant to be used for phasing")
@click.argument("ccs_bam", type=click.Path(exists=True))
@click.argument("subreads_bam", type=click.Path(exists=True))
def cli_handler(directory, prefix, profile, min_ccs_length, max_ccs_length, min_ccs_passes, min_ccs_qual,
                max_homopolymer, trim_ends, tsne_iterations, tsne_rate, cluster_percentile, cluster_inflation,
                cluster_size_threshold, max_cluster_size, consensus_fraction, min_haplotype_molecules,
                min_variant_qual, ccs_bam, subreads_bam,):
    # dict of config values to pass to snakemake
    config = dict(
        PREFIX = prefix,
        MIN_READ_LENGTH = min_ccs_length,
        MAX_READ_LENGTH = max_ccs_length,
        MIN_PASSES = min_ccs_passes,
        MIN_QUAL = min_ccs_qual,
        MAX_HP_LEN = max_homopolymer,
        TRIM_ENDS = trim_ends,
        TSNE_ITERS = tsne_iterations,
        TSNE_LR = tsne_rate,
        CLUSTER_SIMILARITY_PERCENTILE = cluster_percentile,
        CLUSTER_INFLATION = cluster_inflation,
        CLUSTER_SIZE_THRESHOLD=cluster_size_threshold,
        MAX_CLUSTER_SIZE=max_cluster_size,
        CONSENSUS_FRACTION = consensus_fraction,
        MIN_HAPLOTYPE_MOLECULES = min_haplotype_molecules,
        MIN_VARIANT_QUAL = min_variant_qual,
        CCS_BAM = os.path.abspath(ccs_bam),
        SUBREADS_BAM = os.path.abspath(subreads_bam),
        PROFILE = profile,
    )

    config_items = ["{}={}".format(k, v) for k,v in config.items()]
    snakefile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "snakefiles/workflow.snake")

    snake_args = [
        "snakemake", "--rerun-incomplete", "--nolock",
        "--snakefile", snakefile,
        "--directory", directory,
        "--config"
    ] + config_items

    if profile != "":
        snake_args += ["--profile", profile]

    subprocess.run(snake_args)


if __name__ == "__main__":
    cli_handler()
