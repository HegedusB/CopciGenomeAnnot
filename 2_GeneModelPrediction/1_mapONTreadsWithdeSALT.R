library(stringr)


# index ref genome
# deSALT index ../ref_genome/Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.fasta refgen_index/

refgenome_index_path <- "./refgen_index"
fasta_path <- list.files("../lordec_out/", "length_filtered_NB.*_lordec.fasta", full.names = TRUE)

annot_path <- "RoundNrefAnnotFile.gtf"

# ref annotációs fájl elokészítése a deSALT számára
prepannot_command <- str_c("/SSD/software/deSALT/install/deSALT/src/Annotation_Load.py ",annot_path, " r4.annot")
print(prepannot_command)
system(command = prepannot_command)
file.exists("r4.annot")

# fasta_path <- fasta_path[c(3,5)] # ha csak a korábbi hibás illesztéseket akarom megismételni

nCPU <- 48

for(i in fasta_path) {
  out_tag <- str_extract(i, pattern = "NB[0-9]+")
  deSALT_command <- str_c("deSALT aln -t ", nCPU, " -I 1000 -a 6 -k 22 -l 14 -s 2 -x ccs -O6,24 -M4 -G r4.annot -o ", out_tag, "_lordec_desalt_a6k22l14s2.sam ", refgenome_index_path, " ", i)
  # print(deSALT_command)
  system(command = deSALT_command)
}

# sam to sorted bam

sam_path <- list.files(".", pattern = "_lordec_desalt_a6k22l14s2.sam$")
# sam_path <- sam_path[c(3,5)] # ha csak a korábbi hibás illesztéseket akarom megismételni

for(i in sam_path) {
  out_tag <- str_remove(i, pattern = "\\.sam$")
  command <- str_c("samtools view -Sb -@", nCPU, " ", i ," | samtools sort -@", nCPU, " -o ", out_tag, "_s.bam")
  # print(command)
  system(command = command)
}

# index

bam_path <- list.files(".", pattern = "_s.bam$")
# bam_path <- bam_path[c(5,8)] # ha csak a korábbi hibás illesztéseket akarom megismételni

for(i in bam_path) {
  command <- str_c("samtools index -@", nCPU, " ", i)
  system(command = command)
}

# remove reads from sam with long soft clipps (soft clipp filter = sclipf)

samclip_path <- "/work/SSD1_move/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_illumina/additional_apps/samclip"
refgenome_path <- "../ref_genome/Copci_AmutBmut1_AssemblyScaffolds_Repeatmasked.fasta"

for(i in sam_path) {
  out_tag <- str_remove(i, pattern = "\\.sam$")
  command <- str_c(samclip_path, " --ref ", refgenome_path, " --max 25 ", i, " > ", out_tag, "_sclipf.sam")
  # print(command)
  system(command = command)
}

# sam to sorted bam

sclipf_sam_path <- list.files(".", pattern = "_lordec_desalt_a6k22l14s2_sclipf.sam$")
# sclipf_sam_path <- sclipf_sam_path[c(3,5)] # ha csak a korábbi hibás illesztéseket akarom megismételni

for(i in sclipf_sam_path) {
  out_tag <- str_remove(i, pattern = "\\.sam$")
  command <- str_c("samtools view -Sb -@", nCPU, " ", i ," | samtools sort -@", nCPU, " -o ", out_tag, "_s.bam")
  # print(command)
  system(command = command)
}

# index

sclipf_bam_path <- list.files(".", pattern = "_sclipf_s.bam$")
# sclipf_bam_path <- sclipf_bam_path[c(3,5)] # ha csak a korábbi hibás illesztéseket akarom megismételni

for(i in sclipf_bam_path) {
  command <- str_c("samtools index -@", nCPU, " ", i)
  system(command = command)
}


