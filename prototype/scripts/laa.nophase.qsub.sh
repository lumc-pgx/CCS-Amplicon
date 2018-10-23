#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N LAA
#$ -l h_vmem=4G
#$ -pe BWA 1
#$ -cwd
#$ -j Y
#$ -V

PREFIX=$1
SUBREADS=$2
WHITELIST=$3

LAA=/usr/local/smrtlink/smrtcmds/bin/laa
${LAA} -n 1 --ignoreBc --noClustering --noPhasing --resultFile ${PREFIX}.laa.nophase.fastq --junkFile ${PREFIX}.laa.nophase.junk.fastq \
       --reportFile ${PREFIX}.laa.nophase.report.csv --inputReportFile ${PREFIX}.laa.nophase.input.csv \
       --subreadsReportPrefix ${PREFIX}.laa.nophase.subreads --whitelist ${WHITELIST} ${SUBREADS}
