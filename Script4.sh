awk 'BEGIN {OFS="\t"} {print $1, $2-1, $3}' FinDMR.txt > Fin-DMRs.bed

awk 'NR>1' Fin-DMRs.bed > Fin-DMRs_fixed.bed

awk '{print $1, $2, $3}' OFS="\t" Fin-DMRs_fixed.bed > Fin-DMRs_final.bed

awk '$2 >= 0' Fin-DMRs_final.bed > Fin-DMRs_cleaned.bed


bedtools intersect -a Fin-DMRs_cleaned.bed -b genes.bed -wa -wb > Fin_Gene_Overlap.txt
