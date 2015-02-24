#Read mapping
tmap map3 -f reference.fasta \
-r reads.fastq \
-i fastq \
-o 2 \
-s mapped_reads.bam

#Header modifying
samtools view -H mapped_reads.bam \
| sed 's/SM:NOSM/SM:Sample1/' \
| samtools reheader - mapped_reads.bam \
> mapped_reads.reheaded.bam

#Bam file sorting
samtools sort mapped_reads.reheaded.bam mapped_reads.reheaded.sorted
#Variant detection
samtools mpileup -uf reference.fasta mapped_reads.reheaded.sorted.bam \
| bcftools call -cv \
| vcfutils.pl varFilter > variants.samtools.vcf

