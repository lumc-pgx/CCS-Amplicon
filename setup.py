import os
from setuptools import setup

distmeta = {}
for line in open(os.path.join('ccs_amplicon', '__init__.py')):
    try:
        field, value = (x.strip() for x in line.split('='))
    except ValueError:
        continue
    if field == '__version_info__':
        value = value.strip('[]()')
        value = '.'.join(x.strip(' \'"') for x in value.split(','))
    else:
        value = value.strip('\'"')
    distmeta[field] = value

long_description = "See {}".format(distmeta["__homepage__"])

setup(
    name="ccs_amplicon",
    version=distmeta["__version_info__"],
    description="Haplotype polishing and phasing from targetted PacBio sequencing data.",
    long_description=long_description,
    author=distmeta["__author__"],
    author_email=distmeta["__contact__"],
    url=distmeta["__homepage__"],
    license="MIT",
    platforms=["linux"],
    packages=["ccs_amplicon"],
    include_package_data=True,
    python_requires='~=3.0',
    install_requires=[
        "numpy",
        "scikit-learn",
        "pandas",
        "bokeh",
        "pysam",
        "click>=7.0",
        "markov_clustering>=0.0.5.dev0",
        "networkx>=2.0",
        "pygraphviz>=1.4rc1",
        "biopython",
        "edlib",
        "tqdm",
        "snakemake>=5.0.0"
    ],
    entry_points={
        "console_scripts": [
            "ccs_amplicon.workflow = ccs_amplicon.runner:cli_handler",
            "ccs_amplicon.baminfo = ccs_amplicon.baminfo:cli_handler",
            "ccs_amplicon.trim_fasta = ccs_amplicon.trim_fasta:cli_handler",
            "ccs_amplicon.clusters = ccs_amplicon.cluster:cli_handler",
            "ccs_amplicon.collapse = ccs_amplicon.collapse_homopolymers:cli_handler",
            "ccs_amplicon.consensus = ccs_amplicon.consensus:cli_handler",
            "ccs_amplicon.distance = ccs_amplicon.distance_matrix:cli_handler",
            "ccs_amplicon.embed = ccs_amplicon.embed:cli_handler",
            "ccs_amplicon.label = ccs_amplicon.label_haplotype_seqs:cli_handler",
            "ccs_amplicon.sanitize_phase = ccs_amplicon.sanitize_phased_vcf:cli_handler",
            "ccs_amplicon.seqs_from_clusters = ccs_amplicon.seqs_from_clusters:cli_handler",
            "ccs_amplicon.whitelist = ccs_amplicon.whitelist:cli_handler",
            "ccs_amplicon.clusterplot = ccs_amplicon.cluster_view:cli_handler",
            "ccs_amplicon.filter_clusters = ccs_amplicon.filter_clusters:cli_handler",
            "ccs_amplicon.summarize_phasing = ccs_amplicon.phasing_summary:cli_handler",
            "ccs_amplicon.strip_primers = ccs_amplicon.strip_primers:cli_handler",
            "ccs_amplicon.tag_chimeras = ccs_amplicon.chimerism:cli_handler",
            "ccs_amplicon.fq2fa = ccs_amplicon.fq_to_fa:cli_handler"
        ]
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT License",
    ],
    keywords="bioinformatics"
)
