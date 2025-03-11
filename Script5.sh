# Extract matching lines
awk 'NR==FNR {genes[$1]; next} $2 in genes {print $4}' foundgenes.txt OMA_file.tsv | sort | uniq > zebrafish_orthologs.txt

echo "Conversion complete. Zebrafish genes saved"
