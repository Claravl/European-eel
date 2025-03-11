#running basecaller
dorado-0.7.1-linux-x64/bin/dorado basecaller sup,5mCG_5hmCG /lustre/scratch/WUR/ABGC/gunda005/eel_nanopore_data/20240228_Aanguilla_pool_${POOL}/loading_${LOADING}/*/pod5_pass/barcode${
BARCODE}/  --reference ../reference_genome/GCF_013347855.1_fAngAng1.pri_genomic.fna --recursive --device cuda:0 --kit-name SQK-NBD114-
24 >Pool${POOL}_loading${LOADING}_barcode${BARCODE}_met_D071.bam
