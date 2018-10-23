#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N LAA
#$ -l h_vmem=4G
#$ -pe BWA 4
#$ -cwd
#$ -j Y
#$ -V

PREFIX=$1
SUBREADS=$2
WHITELIST=$3

LAA=/usr/local/smrtlink/smrtcmds/bin/laa
${LAA} -n 4 --resultFile ${PREFIX}.laa.fastq --junkFile ${PREFIX}.laa.nk.fastq \
       --reportFile ${PREFIX}.laa.report.csv --inputReportFile ${PREFIX}.laa.input.csv \
       --subreadsReportPrefix ${PREFIX}.laa.subreads --whitelist ${WHITELIST} ${SUBREADS}
