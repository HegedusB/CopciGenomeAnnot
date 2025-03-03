# basecalling (ONT Guppy basecalling software version 3.2.4+d9ed22f)
guppy_basecaller --flowcell FLO-MIN106 --kit SQK-LSK109 --barcode_kits EXP-NBD104 -r --input_path ./raw_reads/ --save_path ./basecalled_demultiplexed_reads/ --num_callers 7 --cpu_threads_per_caller 10 --qscore_filtering

# pychopper (Tool to identify, orient and rescue full-length cDNA reads)
-NB01
cdna_classifier.py -b ../cdna_barcodes_teloprime_NB01.fas -r NB01_report.pdf -S NB01_report.stats -A NB01_report.scores -u NB01_report_unclass.fastq -w NB01_report_rescue.fastq -t 40 -m edlib -c primer_config.txt ../barcode01_all.fastq NB01_report_full_length.fastq

-NB02
cdna_classifier.py -b ../cdna_barcodes_teloprime_NB02.fas -r NB02_report.pdf -S NB02_report.stats -A NB02_report.scores -u NB02_report_unclass.fastq -w NB02_report_rescue.fastq -t 40 -m edlib -c primer_config.txt ../barcode02_all.fastq NB02_report_full_length.fastq

-NB03
cdna_classifier.py -b ../cdna_barcodes_teloprime_NB03.fas -r NB03_report.pdf -S NB03_report.stats -A NB03_report.scores -u NB03_report_unclass.fastq -w NB03_report_rescue.fastq -t 40 -m edlib -c primer_config.txt ../barcode03_all.fastq NB03_report_full_length.fastq

-NB04
cdna_classifier.py -b ../cdna_barcodes_teloprime_NB04.fas -r NB04_report.pdf -S NB04_report.stats -A NB04_report.scores -u NB04_report_unclass.fastq -w NB04_report_rescue.fastq -t 10 -m edlib -c primer_config.txt ../barcode04_all.fastq NB04_report_full_length.fastq

-NB05
cdna_classifier.py -b ../cdna_barcodes_teloprime_NB05.fas -r NB05_report.pdf -S NB05_report.stats -A NB05_report.scores -u NB05_report_unclass.fastq -w NB05_report_rescue.fastq -t 10 -m edlib -c primer_config.txt ../barcode05_all.fastq NB05_report_full_length.fastq

# Illumina reads subsampling for ONT polishing

all_R1.fastq = "It contains all the R1 read pairs of the samples in the GSE125200 dataset"

reformat.sh in=all_R1.fastq out=subsampled.fastq samplerate=0.1

ORNA -input subsampled_10p.fastq -output subsampled_10p_norm.fastq -base 1.7 -kmer 21 -nb-cores 20 -type fasta 2> run.info


# ONT polishing using lordec
singularity run -B /SSD1:/SSD1 -B /SSD2:/SSD2 docker://hbotond/lordec:latest

readlink -f subsampled_10p_norm.fastq.fa > meta_for_lordec.tsv

-NB01/VM
lordec-correct -i ../all_pychopper_pass_reads/length_filtered_NB01.fastq -2 meta_for_lordec.tsv -k 19 -s 3 -T 70 -o length_filtered_NB01_lordec.fasta

-NB02/P1
lordec-correct -i ../all_pychopper_pass_reads/length_filtered_NB02.fastq -2 meta_for_lordec.tsv -k 19 -s 3 -T 70 -o length_filtered_NB02_lordec.fasta

-NB03/P2
lordec-correct -i ../all_pychopper_pass_reads/length_filtered_NB03.fastq -2 meta_for_lordec.tsv -k 19 -s 3 -T 70 -o length_filtered_NB03_lordec.fasta

-NB04/YC
lordec-correct -i ../all_pychopper_pass_reads/length_filtered_NB04.fastq -2 meta_for_lordec.tsv -k 19 -s 3 -T 70 -o length_filtered_NB04_lordec.fasta

-NB05/YS
lordec-correct -i ../all_pychopper_pass_reads/length_filtered_NB05.fastq -2 meta_for_lordec.tsv -k 19 -s 3 -T 70 -o length_filtered_NB05_lordec.fasta

