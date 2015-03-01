#!/bin/bash

#Setting sample name
sample_name=$1

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
  | bcftools call -mv --skip-variants indels \
  > $sample_name.bcftools_snps.vcf

#SNP annotating with snpEff
java -jar ~/programs/snpEff/snpEff.jar \
  ann \
  -c ~/programs/snpEff/snpEff.config \
  GRCh38.76 \
  $sample_name.bcftools_snps.vcf \
  > $sample_name.bcftools_snps.annotated.vcf

#filtering of missense and stop SNPs
java -jar ~/programs/snpEff/SnpSift.jar \
  filter \
  -f $sample_name.bcftools_snps.annotated.vcf \
  "ANN[*].EFFECT='missense_variant'||ANN[*].EFFECT='stop_gained'"
  > $sample_name.bcftools_snps.annotated.nonsyn_stop.vcf
