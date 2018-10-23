SAMPLE=$1

SMRTPATH=/usr/local/smrtlink/install/smrtlink-release_5.1.0.26412/bundles/smrttools/smrtcmds/bin/

${SMRTPATH}/pbindex ${SAMPLE}.bam
${SMRTPATH}/bam2fastq -o ${SAMPLE} ${SAMPLE}.bam
gunzip ${SAMPLE}.fastq.gz 
ngmlr -r ${SAMPLE}.consensus.fa -q ${SAMPLE}.fastq --rg-id foo --rg-sm sm --rg-lb lb > ${SAMPLE}.aligned.sam
samtools sort ${SAMPLE}.aligned.sam > ${SAMPLE}.aligned.bam
samtools index ${SAMPLE}.aligned.bam
bcftools mpileup -f ${SAMPLE}.consensus.fa ${SAMPLE}.aligned.bam | bcftools call -mv -Ob -o ${SAMPLE}.bcf
bcftools view -i '%QUAL>=50' ${SAMPLE}.bcf > ${SAMPLE}.vcf

#whatshap phase --indels --reference ${SAMPLE}.consensus.fa ${SAMPLE}.vcf ${SAMPLE}.aligned.bam -o ${SAMPLE}.phased.vcf
whatshap phase --indels ${SAMPLE}.vcf ${SAMPLE}.aligned.bam -o ${SAMPLE}.phased.vcf
bgzip ${SAMPLE}.phased.vcf
tabix ${SAMPLE}.phased.vcf.gz

#whatshap haplotag -o ${SAMPLE}.aligned.tagged.bam -r ${SAMPLE}.consensus.fa ${SAMPLE}.phased.vcf.gz ${SAMPLE}.aligned.bam
whatshap haplotag -o ${SAMPLE}.aligned.tagged.bam ${SAMPLE}.phased.vcf.gz ${SAMPLE}.aligned.bam
samtools index ${SAMPLE}.aligned.tagged.bam


# notes
# if the unphased vcf contains a single heterozygous variant, whatshap will not
# phase it. This can be overcome by manually adding the phasing information to the phased.vcf
# once this is done, haplotag will correctly assign reads based on the phase of the variant.
