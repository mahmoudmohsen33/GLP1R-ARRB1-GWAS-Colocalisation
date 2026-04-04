# ================================
# 01_read_magic_vcf.R
# Read and parse MAGIC GWAS VCF files
# ================================

library(data.table)
library(stringr)

# Function: read MAGIC VCF
read_magic_vcf <- function(vcf_file){

  dt <- fread(
    vcf_file,
    sep = "\t",
    header = TRUE,
    skip = "#CHROM",
    data.table = FALSE,
    quote = ""
  )

  # Rename columns
  names(dt)[1:5] <- c("chrom","pos","rsid","ref","alt")

  # Extract FORMAT structure
  format_vec <- strsplit(dt$FORMAT[1], ":")[[1]]

  sample_mat <- str_split_fixed(dt[[10]], ":", length(format_vec))
  colnames(sample_mat) <- format_vec

  df <- cbind(
    dt[, c("chrom","pos","rsid","ref","alt")],
    as.data.frame(sample_mat)
  )

  # Convert numeric fields
  num_fields <- intersect(c("ES","SE","LP","AF","SS"), colnames(df))
  for (fld in num_fields){
    df[[fld]] <- suppressWarnings(as.numeric(df[[fld]]))
  }

  # Compute p-values
  if ("LP" %in% colnames(df)){
    df$pval <- 10^(-df$LP)
  } else {
    df$pval <- NA
  }

  return(df)
}
