#run modkit
for barcode in 13 14 15; do
    modkit pileup --cpg \
    --ref ~/reference_Eel/GCF_013347855.1_fAngAng1.pri_genomic.fna \
    --ignore h \
    --log-filepath ~/Eel-Data/log_barcode${barcode}.txt \
    --combine-strands /lustre/nobackup/WUR/ABGC/lynde002/Eel-Data/merged_pool19_barcode${barcode}_sorted.bam \
    /lustre/nobackup/WUR/ABGC/lynde002/Eel-Data/merged_pool19_barcode${barcode}_modkit.txt
done
