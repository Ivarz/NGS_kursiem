#!/bin/bash

#Read mapping
tmap map3 -f ref/chr_merged.fa \
  -r reseq_reads.fastq \
  -i fastq \
  -R ID:SomeID \
  -R SM:Sample1 \
  -o 1 \
  -s reseq_reads_mapped.reheaded.bam

#Read sorting
samtools sort reseq_reads_mapped.reheaded.bam \
reseq_reads_mapped.reheaded.sorted

#SNP calling with bcftools
samtools mpileup -uf ref/chr_merged.fa \
  reseq_reads_mapped.reheaded.sorted.bam \
  | bcftools call -mv --skip-variants indels \
  > reseq_reads.bcftools_snps.vcf

#SNP annotating with snpEff
java -jar ~/programs/snpEff/snpEff.jar \
  ann \
  -c ~/programs/snpEff/snpEff.config \
  GRCh38.76 \
  reseq_reads.bcftools_snps.vcf \
  > reseq_reads.bcftools_snps.annotated.vcf

#filtering of missense and stop SNPs
java -jar ~/programs/snpEff/SnpSift.jar \
  filter \
  -f reseq_reads.bcftools_snps.annotated.vcf \
  "ANN[*].EFFECT='missense_variant'||ANN[*].EFFECT='stop_gained'"
  > reseq_reads.bcftools_snps.annotated.nonsyn_stop.vcf
