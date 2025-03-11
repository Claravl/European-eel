# Load libraries
library(methylKit)
library(data.table)
library(GenomicRanges)

# Define file names and corresponding sample IDs
files <- c(
  "merged_pool17_barcode18_modkit.txt",
  "merged_pool19_barcode13_modkit.txt",
  "merged_pool19_barcode14_modkit.txt",
  "merged_pool19_barcode15_modkit.txt",
  "merged_pool11_barcode11_modkit.txt",
  "merged_pool14_barcode05_modkit.txt",
  "merged_pool11_barcode10_modkit.txt",
  "merged_pool11_barcode09_modkit.txt"
)

sample_ids <- c(
  "17-bc18", "19-bc13", "19-bc14", "19-bc15",
  "11-bc11", "14-bc05", "11-bc10", "11-bc09"
)

#define treatment vector
treatment_vector <- c(1,1,1,1,0,0,0,0)

# Load the GRanges object with predefined regions
gene_annotations <- readRDS("gene_annotations.rds")  

TSS_regions <- promoters(gene_annotations, upstream = 2000, downstream = 500)

# empty list to store methylation objects
meth_objects <- list()

# Process each file and create methylation objects
for (i in seq_along(files)) {
  # Read and modify the modkit file
  modkit_data <- fread(files[i], header = FALSE)
  modkit_data[, V11 := V11 / 100]  # Convert percentages to fractions
  modified_file <- paste0("modified_", files[i])
  fwrite(modkit_data, modified_file, sep = "\t", col.names = FALSE)  # Save modified file

  # Create a methylation object
  meth_obj <- methRead(
    modified_file,
    sample.id = sample_ids[i],
    assembly = "eel_genome",  # genome assembly name
    header = FALSE,
    context = "CpG",
    resolution = "base",      # Base-level resolution
    pipeline = list(
      fraction = TRUE,          # Use methylation fraction
      chr.col = 1,              # Chromosome column
      start.col = 2,            # Start position
      end.col = 3,              # End position
      coverage.col = 10,        # Coverage column
      strand.col = 6,           # Strand column
      freqC.col = 11            # Methylation fraction column
    )
  )

  # Append the methylation object to the list
  meth_objects[[i]] <- meth_obj
}

# Combine the individual methylRaw objects into a methylRawList
meth_raw_list <- methylRawList(meth_objects, treatment = treatment_vector)

# Filter the data 
filtered_meth <- filterByCoverage(meth_raw_list, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.5)

# Combine filtered data
meth_combined <- unite(filtered_meth)

# Perform PCA
pdf("PCA.pdf")
PCASamples(meth_combined)  
dev.off()

# Aggregate methylation data over predefined regions
meth_regions <- regionCounts(
  obj = meth_combined,
  regions = TSS_regions,    # Use predefined GRanges regions
  min.coverage = 5,    # Minimum coverage 
  context = "CpG"
)

# Omit empty values
meth_regions <- na.omit(meth_regions)

# Perform differential methylation analysis
myDiff <- calculateDiffMeth(meth_regions)

# Get all differentially methylated regions with a 5% difference and q-value < 0.01
myDiff25p <- getMethylDiff(myDiff, difference = 5, qvalue = 0.01)

# Save the results for later review
write.table(
  myDiff25p,
  file = "FinDMR2.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Save the R object for further use in R
saveRDS(myDiff25p, file = "myDiff25p_results.rds")

# Output a message
cat("Differential methylation analysis complete. Significant results saved.\n")

# Compute the distribution of differentially methylated regions per chromosome
diff_per_chr <- diffMethPerChr(
  myDiff,
  plot = FALSE,        
  qvalue.cutoff = 0.01, # Use q-value cutoff of 0.01
  meth.cutoff = 5      # Use methylation difference cutoff of 5%
)

# Save the results to a tab-delimited text file
write.table(
  diff_per_chr,
  file = "FinChrDMR2.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Save the R object for further analysis
saveRDS(diff_per_chr, file = "DiffMethPerChr_Results.rds")

# Output a message
cat("Differential methylation per chromosome results saved.\n")
                                                                 
