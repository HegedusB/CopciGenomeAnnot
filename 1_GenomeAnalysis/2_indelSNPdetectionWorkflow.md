# raw reads qc
mkdir raw_reads_QC
fastqc SRR10162428_1.fastq.gz SRR10162428_2.fastq.gz -t 10 -o raw_reads_QC

# trim
fastp -i SRR10162428_1.fastq.gz -I SRR10162428_2.fastq.gz -o SRR10162428_1_t.fastq.gz -O SRR10162428_2_t.fastq.gz --thread 10

# trimmed reads qc
mkdir timmed_reads_QC
fastqc SRR10162428_1_t.fastq.gz SRR10162428_2_t.fastq.gz -t 10 -o timmed_reads_QC

# index genome 

## CopciAB JGI new
bwa index CopciAB_new_jgi_20220113.fasta
samtools faidx CopciAB_new_jgi_20220113.fasta
java -jar /SSD/software/picard/picard.jar CreateSequenceDictionary REFERENCE=CopciAB_new_jgi_20220113.fasta OUTPUT=CopciAB_new_jgi_20220113.dict

## CopciAB JGI old
bwa index Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.fasta
samtools faidx Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.fasta
java -jar /SSD/software/picard/picard.jar CreateSequenceDictionary REFERENCE=Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.fasta OUTPUT=Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.dict

## CopciAB Xie
cp Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fsa Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fasta
bwa index Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fasta
samtools faidx Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fasta
java -jar /SSD/software/picard/picard.jar CreateSequenceDictionary REFERENCE=Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fasta OUTPUT=Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.dict

# map reads to reference

## CopciAB JGI new
bwa mem -t 20 CopciAB_new_jgi_20220113.fasta  ../genome_reseq_reads/SRR10162428_1_t.fastq.gz ../genome_reseq_reads/SRR10162428_2_t.fastq.gz | samtools view -Sb - > CopciAB_new_jgi_20220113.bam
samtools sort -@ 20 CopciAB_new_jgi_20220113.bam -o CopciAB_new_jgi_20220113_s.bam
samtools index -@ 20 CopciAB_new_jgi_20220113_s.bam

## CopciAB JGI old
bwa mem -t 20 Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.fasta  ../genome_reseq_reads/SRR10162428_1_t.fastq.gz ../genome_reseq_reads/SRR10162428_2_t.fastq.gz | samtools view -Sb - > Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.bam
samtools sort -@ 20 Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.bam -o Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_s.bam
samtools index -@ 20 Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_s.bam

## CopciAB Xie
bwa mem -t 20 Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fasta  ../genome_reseq_reads/SRR10162428_1_t.fastq.gz ../genome_reseq_reads/SRR10162428_2_t.fastq.gz | samtools view -Sb - > Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.bam
samtools sort -@ 20 Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.bam -o Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_s.bam
samtools index -@ 20 Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_s.bam

# gatk reindexing

## CopciAB JGI new
java -Xms1g -Xmx2000m -jar /SSD/software/picard/picard.jar AddOrReplaceReadGroups I=CopciAB_new_jgi_20220113_s.bam O=CopciAB_new_jgi_20220113_s_rg.bam SO=coordinate RGID=S205.1 RGSM=S205.1 RGLB=Lib1 RGPU=Barcode RGPL=ILLUMINA
samtools index CopciAB_new_jgi_20220113_s_rg.bam

## CopciAB JGI old
java -Xms1g -Xmx2000m -jar /SSD/software/picard/picard.jar AddOrReplaceReadGroups I=Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_s.bam O=Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_s_rg.bam SO=coordinate RGID=S205.1 RGSM=S205.1 RGLB=Lib1 RGPU=Barcode RGPL=ILLUMINA
samtools index Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_s_rg.bam

## CopciAB Xie
java -Xms1g -Xmx2000m -jar /SSD/software/picard/picard.jar AddOrReplaceReadGroups I=Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_s.bam O=Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_s_rg.bam SO=coordinate RGID=S205.1 RGSM=S205.1 RGLB=Lib1 RGPU=Barcode RGPL=ILLUMINA
samtools index Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_s_rg.bam

# perform GATK variant call, only INDEL

## CopciAB JGI new
java -jar /home/balintb/GenomeAnalysisTK3.jar -glm BOTH -T UnifiedGenotyper -R CopciAB_new_jgi_20220113.fasta -I CopciAB_new_jgi_20220113_s_rg.bam  -o CopciAB_new_jgi_20220113.vcf

java -jar /home/balintb/GenomeAnalysisTK3.jar -glm INDEL -T UnifiedGenotyper -R CopciAB_new_jgi_20220113.fasta -I CopciAB_new_jgi_20220113_s_rg.bam  -o CopciAB_new_jgi_20220113_indel.vcf

java -jar /home/balintb/GenomeAnalysisTK3.jar -glm SNP -T UnifiedGenotyper -R CopciAB_new_jgi_20220113.fasta -I CopciAB_new_jgi_20220113_s_rg.bam  -o CopciAB_new_jgi_20220113_SNP.vcf

## CopciAB JGI old
java -jar /home/balintb/GenomeAnalysisTK3.jar -glm BOTH -T UnifiedGenotyper -R Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.fasta -I Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_s_rg.bam  -o Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.vcf

java -jar /home/balintb/GenomeAnalysisTK3.jar -glm INDEL -T UnifiedGenotyper -R Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.fasta -I Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_s_rg.bam  -o Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_indel.vcf

java -jar /home/balintb/GenomeAnalysisTK3.jar -glm SNP -T UnifiedGenotyper -R Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.fasta -I Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_s_rg.bam  -o Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked_SNP.vcf

## CopciAB Xie
java -jar /home/balintb/GenomeAnalysisTK3.jar -glm BOTH -T UnifiedGenotyper -R Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fasta -I Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_s_rg.bam  -o Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.vcf

java -jar /home/balintb/GenomeAnalysisTK3.jar -glm INDEL -T UnifiedGenotyper -R Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fasta -I Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_s_rg.bam  -o Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_indel.vcf

java -jar /home/balintb/GenomeAnalysisTK3.jar -glm SNP -T UnifiedGenotyper -R Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fasta -I Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_s_rg.bam  -o Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_SNP.vcf




