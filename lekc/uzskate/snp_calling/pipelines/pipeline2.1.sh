#!/bin/bash

#Setting sample name
sample_name=reseq_reads.fastq

#Read mapping
tmap map3 -f ref/chr_merged.fa \
  -r $sample_name \
  -i fastq \
  -R ID:SomeID \
  -R SM:Sample1 \
  -o 1 \
  -s $sample_name.mapped.reheaded.bam

#Read sorting
samtools sort $sample_name.mapped.reheaded.bam \
$sample_name.mapped.reheaded.sorted

#SNP calling with bcftools
samtools mpileup -uf ref/chr_merged.fa \
  $sample_name.mapped.reheaded.sorted.bam \
  | bcftools call --multiallelic-caller --variants-only --skip-variants indels \
  > $sample_name.bcftools_snps.vcf

#SNP annotating with snpEff
java -jar -Xmx3g ~/programs/snpEff/snpEff.jar \
  ann \
  -c ~/programs/snpEff/snpEff.config \
  GRCh38.76 \
  $sample_name.bcftools_snps.vcf \
  > $sample_name.bcftools_snps.annotated.vcf

java -jar ~/programs/snpEff/SnpSift.jar \
  annotate \
  -id \
  ~/data/day1/snpEffData/00-All.vcf \
  reseq_reads.bcftools_snps.annotated.vcf \
  > reseq_reads.bcftools_snps.annotated.id.vcf

#filtering of missense and stop SNPs
java -jar -Xmx3g ~/programs/snpEff/SnpSift.jar \
  filter \
  -f $sample_name.bcftools_snps.annotated.id.vcf \
  "(ANN[*].EFFECT='missense_variant')||(ANN[*].EFFECT='stop_gained')" \
  > $sample_name.bcftools_snps.annotated.id.nonsyn_stop.vcf
