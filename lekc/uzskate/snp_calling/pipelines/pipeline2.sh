#!/bin/bash

#Setting sample name
sample_name=reseq_reads

#Read mapping
tmap map3 -f ref/chr_merged.fa \
  -r $sample_name.fastq \
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

#SNP annotating with snpEff and SnpSift
java -jar ~/programs/snpEff/snpEff.jar \
  ann \
  -c ~/programs/snpEff/snpEff.config \
  GRCh38.76 \
  $sample_name.bcftools_snps.vcf \
  > $sample_name.bcftools_snps.annotated.vcf

java -jar ~/programs/snpEff/SnpSift.jar \
  annotate \
  -id \
  ~/data/day1/snpEffData/00-All.vcf \
  $sample_name.bcftools_snps.annotated.vcf \
  > $sample_name.bcftools_snps.annotated.id.vcf

java -jar ~/programs/snpEff/SnpSift.jar \
  dbnsfp \
  $sample_name.bcftools_snps.annotated.id.vcf \
  > $sample_name.bcftools_snps.annotated.id.dbnsfp.vcf

#filtering of missense and stop SNPs
java -jar ~/programs/snpEff/SnpSift.jar \
  filter \
  -f $sample_name.bcftools_snps.annotated.id.dbnsfp.vcf \
  "(ANN[*].EFFECT='missense_variant')||(ANN[*].EFFECT='stop_gained')" \
  > $sample_name.bcftools_snps.annotated.id.dbnsfp.nonsyn_stop.vcf

#Table creation
java -jar ~/programs/snpEff/SnpSift.jar
	extractFields \
	$sample_name.bcftools_snps.annotated.id.dbnsfp.nonsyn_stop.vcf \
	CHROM POS ID REF ALT dbNSFP_SIFT_pred GEN[*].GT \
	> $sample_name.bcftools_snps.annotated.id.dbnsfp.nonsyn_stop.txt
