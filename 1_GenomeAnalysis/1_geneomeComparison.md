# unify genome seq format 

```R
library(Biostrings)
library(readr)

cons_scaf <- read_tsv("consensus_scaffold_names.tsv")

# Okayama
refO <- readDNAStringSet("Coprinopsis_cinerea.masked.fasta")
refO13 <- refO[grepl(names(refO), pattern="Chr")]
rev_index <- names(refO13) %in% c("Chr_1","Chr_4", "Chr_5", "Chr_6", "Chr_7", "Chr_9")
refO13[rev_index] <- reverseComplement(refO13[rev_index])
names(refO13) <- cons_scaf$Consensus_name[match(names(refO13), cons_scaf$Okayama)]
writeXStringSet(refO13, "Coprinopsis_cinerea.masked_C13_revComp.fasta")

# JGI V2
refJGI <- readDNAStringSet("CopciAB_new_jgi_20220113.fasta")
refJGI13 <- refJGI[1:13]
writeXStringSet(refJGI13, "CopciAB_new_jgi_20220113_C13.fasta")

# ONT/Xie
refXie <- readDNAStringSet("Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fasta")
refXie13 <- refXie[1:13]
names(refXie13) <- cons_scaf$Consensus_name[match(names(refXie13), cons_scaf$Xie)]
writeXStringSet(refXie13, "Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_C13.fasta")

```

# genome comparison --------------------------------------------------------

# nucmer
# comp AB
nucmer --maxmatch -c 500 -b 500 -l 100 --prefix A_B -t 10 Coprinopsis_cinerea.masked_C13_revComp.fasta CopciAB_new_jgi_20220113_C13.fasta
delta-filter -m -i 90 -l 100 A_B.delta > A_B_m_i90_l100.delta; 
show-coords -THrd A_B_m_i90_l100.delta > A_B_m_i90_l100.coords;

# comp BC
nucmer --maxmatch -c 500 -b 500 -l 100 --prefix B_C -t 10 CopciAB_new_jgi_20220113_C13.fasta Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_C13.fasta
delta-filter -m -i 90 -l 100 B_C.delta > B_C_m_i90_l100.delta; 
show-coords -THrd B_C_m_i90_l100.delta > B_C_m_i90_l100.coords;

# plot --------------------------------------------------
syri -c A_B_m_i90_l100.coords -d A_B_m_i90_l100.delta -r Coprinopsis_cinerea.masked_C13_revComp.fasta -q CopciAB_new_jgi_20220113_C13.fasta --prefix Am_Bm
syri -c B_C_m_i90_l100.coords -d B_C_m_i90_l100.delta -r CopciAB_new_jgi_20220113_C13.fasta -q Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs_C13.fasta --prefix Bm_Cm

plotsr --sr Am_Bmsyri.out --sr Bm_Cmsyri.out --genomes genomes.txt --chrord scaffold_order.txt --cfg base.cfg -o output_plotX5.pdf -W 16 -H 5



