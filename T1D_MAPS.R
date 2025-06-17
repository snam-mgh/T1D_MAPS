library(data.table)
library(tools)
library(R.utils)

# Function to extract HLA rows from VCF
extract_HLA_rows <- function(vcf_file, output_file) {
  df <- fread(vcf_file)
  df <- data.frame(df)
  HLA_df <- df[grepl('HLA', df$ID), ]
  names(HLA_df)[10:ncol(HLA_df)] <- sub('.','', names(HLA_df)[10:ncol(HLA_df)])
  names(HLA_df)[1] <- 'CHROM'
  write.table(HLA_df, file=output_file, row.names=FALSE, quote=FALSE)
  return(output_file)
}

# Function to extract GT and DS
extract_GT_DS <- function(hla_file, working_dir, gt_file, ds_file) {
  df <- fread(hla_file)
  df <- data.frame(df)
  names(df)[10:ncol(df)] <- sub('.','', names(df)[10:ncol(df)])
  HLA_list <-  read.csv(paste0(working_dir, 'HLA_allele_list_for_T1D_PRS.txt'), head=F)
  df <- df[df$ID %in% HLA_list$V1, ]
  GT_df <- DS_df <- df
  GT_df$FORMAT <- 'GT'
  DS_df$FORMAT <- 'DS'
  
  for (i in 1:nrow(df)) {
    for (j in 10:ncol(df)) {
      split_vals <- strsplit(as.character(df[i, j]), ":", fixed=TRUE)[[1]]
      GT_df[i, j] <- split_vals[1]
      DS_df[i, j] <- split_vals[2]
    }
  }
  
  write.table(GT_df, file=gt_file, row.names=FALSE, quote=FALSE)
  write.table(DS_df, file=ds_file, row.names=FALSE, quote=FALSE)
}

# Function to create haplotypes
generate_haplotypes <- function(gt_file, working_dir, output_file) {
  part1_df <- fread(gt_file)
  transpose_df <- transpose(part1_df)
  transpose_df <- data.frame(transpose_df)
  transpose_df <- cbind(names(part1_df), transpose_df)
  names(transpose_df) <- c('ID', transpose_df[3, 2:ncol(transpose_df)])
  transpose_df <- transpose_df[10:nrow(transpose_df), ]
  
  mat_chr_df <- transpose_df
  pat_chr_df <- transpose_df
  for (i in 2:ncol(transpose_df)) {
    curr_allele <- data.frame(do.call('rbind', strsplit(as.character(transpose_df[, i]), '|', fixed=TRUE)))
    mat_chr_df[, i] <- as.numeric(curr_allele[, 1])
    pat_chr_df[, i] <- as.numeric(curr_allele[, 2])
  }
  
  haplotype_table <- read.table(paste0(working_dir, 'HLA_haplotype_table.txt'), head=T)
  
  mat_hap_df <- data.frame(ID = mat_chr_df[, 1])
  pat_hap_df <- data.frame(ID = pat_chr_df[, 1])
  
  for (i in 1:nrow(haplotype_table)) {
    DRB1 <- haplotype_table[i, 'DRB1_code']
    DQA1 <- haplotype_table[i, 'DQA1_code']
    DQB1 <- haplotype_table[i, 'DQB1_code']
    
    mat_hap_df[, i + 1] <- ifelse(mat_chr_df[, DRB1] == 1 & mat_chr_df[, DQA1] == 1 & mat_chr_df[, DQB1] == 1, 1, 0)
    pat_hap_df[, i + 1] <- ifelse(pat_chr_df[, DRB1] == 1 & pat_chr_df[, DQA1] == 1 & pat_chr_df[, DQB1] == 1, 1, 0)
  }
  
  total_df <- data.frame(ID = mat_chr_df[, 1])
  total_df[, 2:(nrow(haplotype_table) + 1)] <- mat_hap_df[, -1] + pat_hap_df[, -1]
  names(total_df)[2:(nrow(haplotype_table) + 1)] <- paste0('Hap', 1:nrow(haplotype_table))
  
  write.table(total_df, file=output_file, row.names=FALSE, quote=FALSE)
}

# Function to score PRS
calculate_PRS <- function(haplotype_file, working_dir, output_file) {
  HLA_df <- read.csv(haplotype_file, sep='')
  haplotype_table <- read.table(paste0(working_dir, 'HLA_haplotype_table.txt'), head=T)
  weights <- haplotype_table$prs_weights
  
  num_haps <- sum(grepl("^Hap[0-9]+$", names(HLA_df)))
  if (length(weights) != num_haps) {
    stop("Number of weights does not match number of Hap columns in HLA_df.")
  }
  
  hap_columns <- paste0("Hap", 1:num_haps)
  hap_matrix <- as.matrix(HLA_df[, hap_columns])
  
  HLA_df$T1D_PRS_HLA <- hap_matrix %*% weights
  write.table(HLA_df[, c('ID', 'T1D_PRS_HLA')], file=output_file, sep='\t', row.names=FALSE, quote=FALSE)
}

combine_HLA_nonHLA_scores <- function(hla_score_file, nonhla_score_file, output_file) {
  hla_prs <- fread(hla_score_file)
  nonhla_prs <- fread(nonhla_score_file)
  hla_prs[, ID := gsub("\\.", "-", ID)]
  merged <- merge(hla_prs, nonhla_prs, by.x = "ID", by.y = "IID")
  
  names(merged)[names(merged) == "T1D_PRS_HLA"] <- "HLA_score"
  score_col <- grep("^SCORE", names(merged), value = TRUE) 
  if (length(score_col) != 1) {paste(stop("Unable to uniquely identify score column in PLINK output. Found: ", paste(score_col, collapse=", ")))
  }  
  names(merged)[names(merged) == score_col] <- "nonHLA_score"
    
  merged$Total_score <- merged$HLA_score * 0.79003 + merged$nonHLA_score * 8.86562 + 20
  fwrite(merged[, .(ID, HLA_score, nonHLA_score, Total_score)], file = output_file, sep = "\t", quote = FALSE)
}

# --------- Argument Handling ---------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5 || length(args) == 6 || length(args) > 7) {
  stop("Usage: Rscript T1D_MAPS.R <vcf.gz> <BFILE> <genome_build (hg19|hg38)> <path_to_PLINK> <working_dir> [genetic_prob_file] [t1d_grs2_file]")}

input_vcf <- args[1]
input_plink <- args[2]
genome_build <- args[3]
path_to_PLINK <- args[4]
working_dir <- args[5]
genetic_prob_file <- if (length(args) >= 6 && args[6] != "NA") args[6] else NA
t1d_grs2_file <- if (length(args) == 7 && args[7] != "NA") args[7] else NA

if (!grepl("/$", working_dir)) {
  working_dir <- paste0(working_dir, "/")}    
    
if (xor(is.na(genetic_prob_file), is.na(t1d_grs2_file))) {
  stop("Both genetic_prob_file and t1d_grs2_file must be provided together, or both set to NA.")}

bfile_base <- basename(input_plink)

# --------- Temp and Output Filenames ---------
hla_vcf <- sub("\\.gz$", "", input_vcf)
hla_gt <- paste0(working_dir, "T1DMAPS_", bfile_base, ".GT")
hla_ds <- paste0(working_dir, "T1DMAPS_", bfile_base, ".DS")
hap_out <- paste0(working_dir, "T1DMAPS_", bfile_base, ".haplotypes")
prs_out <- paste0(working_dir, "T1DMAPS_", bfile_base, ".haplotypes_scored.txt")
final_score_file <- paste0(working_dir, "T1DMAPS_", bfile_base, ".total_score.txt")

# --------- Run T1D MAPS Pipeline ---------
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

log_msg("Starting extraction of HLA rows from VCF file...\n")
hla_vcf_file <- extract_HLA_rows(input_vcf, hla_vcf)
log_msg("Completed extraction of HLA rows.\n")
log_msg("Extracting genotype (GT) and dosage (DS) data...\n")
extract_GT_DS(hla_vcf_file, working_dir, hla_gt, hla_ds)
log_msg("Completed extraction of GT and DS data.\n")
log_msg("Generating HLA haplotypes...\n")
generate_haplotypes(hla_gt, working_dir, hap_out)
log_msg("HLA haplotype generation completed.\n")
log_msg("Calculating HLA polygenic risk scores (PRS)...\n")
calculate_PRS(hap_out, working_dir, prs_out)
log_msg("HLA PRS calculation completed.\n")


# --------- PLINK Non-HLA PRS Scoring ---------
score_col_chr <- if (genome_build == "hg38") 6 else if (genome_build == "hg19") 3 else stop("Invalid genome build. Use 'hg19' or 'hg38'.")

gunzip(paste0(working_dir, "PRScs_EUR_phi1e-05_weights.txt.gz"),overwrite = FALSE,remove = FALSE)
scorefile <- paste0(working_dir, "PRScs_EUR_phi1e-05_weights.txt")
outfile <- paste0(working_dir, "T1DMAPS_", bfile_base, ".txt")

# Ensure .bim file SNP IDs are in chr:pos format
bim_file <- paste0(input_plink, ".bim")
bim <- fread(bim_file, header = FALSE)
bim[, V2 := sub('^"?chr"?', '', V2)]                     # Remove 'chr' or "chr"
bim[, V2 := sub('^(\\d+):(\\d+):.*', '\\1:\\2', V2)]     # Keep only chrom:pos if alleles exist
bim[, V2 := sub('^(\\d+):(\\d+)$', '\\1:\\2', V2)]       # If already in chrom:pos, leave it
bim[, V1 := sub('^"?chr"?', '', V1)] # In case "chr" in front
bim[, V2 := ifelse(!grepl('^\\d+:\\d+$', V2), paste0(V1, ":", V4), V2)] # If V2 still doesn't match pattern, replace with chrom:pos using V1 and V4
fwrite(bim, bim_file, sep = "\t", quote = FALSE, col.names = FALSE)

plink_cmd <- paste(
  shQuote(path_to_PLINK),
  "--bfile", shQuote(input_plink),
  "--score", shQuote(scorefile), score_col_chr, "8 10 sum",
  "--out", shQuote(outfile)
)

log_msg("Running PLINK non-HLA scoring step...\n")
if (.Platform$OS.type == "windows") {
  shell(plink_cmd)
} else {
  system(plink_cmd)
}
log_msg("PLINK non-HLA scoring step completed.\n")

# --------- Combine HLA and non-HLA Scores ---------
log_msg("Combining HLA and non-HLA scores...\n")
if (!grepl("\\.profile$", outfile)) {
  outfile2 <- paste0(outfile, ".profile")
} else {
  outfile2 <- outfile
}
combine_HLA_nonHLA_scores(prs_out, outfile2, final_score_file)
log_msg("HLA and non-HLA score combination completed.\n")

# --------- Calculate T1D_MAPS2_HLA and T1D_MAPS2 ---------
if (!is.na(genetic_prob_file) && !is.na(t1d_grs2_file)) {
  log_msg("Calculating T1D_MAPS2_HLA and T1D_MAPS2 scores...\n")
  hla_data <- fread(t1d_grs2_file)
  ancestry_data <- fread(genetic_prob_file)

  ancestry_data[, IID := as.character(IID)]
  hla_data[, IID := as.character(IID)]
  
  merged_data <- merge(hla_data, ancestry_data[, .(IID, prob_EUR)], by = "IID", all = FALSE)
  
  if (any(is.na(merged_data$prob_EUR))) {
    stop("Some IDs in HLA file do not have matching ancestry probabilities.")}
  
  final_scores <- fread(paste0(working_dir, "T1DMAPS_", bfile_base, ".total_score.txt"))
  final_scores[, match_key := sub(".*-", "", ID)]
  merged_data[, match_key := sub(".*-", "", IID)]
  merged_data <- merge(merged_data, final_scores[, .(match_key, HLA_score, nonHLA_score, Total_score)], by = "match_key")
  
  merged_data[, T1D_GRS2_HLA := HLA_DRDQ + HLA_Class_1 + HLA_Class_2]
  merged_data[, T1D_MAPS2_HLA_score := prob_EUR * T1D_GRS2_HLA + (1 - prob_EUR) * (HLA_score+8.6226)]
  
  merged_data[, T1D_MAPS2_score := 0.82424 * merged_data$T1D_MAPS2_HLA_score + 8.94820 * merged_data$nonHLA_score + 13]
  
  merged_clean <- merged_data[, .(IID, T1D_MAPS2_HLA_score = T1D_MAPS2_HLA_score,
                                  T1D_MAPS2_score = T1D_MAPS2_score)]
  final_scores <- final_scores[, .(match_key, ID, T1D_MAPS_HLA_score = HLA_score,
                                   T1D_MAPS_nonHLA_score = nonHLA_score, 
                                   T1D_MAPS_Total_score = Total_score)]
  final_output <- merge(final_scores, merged_clean, by.x = "match_key", by.y = "IID", all.x = TRUE)
  final_output[, match_key := NULL]
    
  fwrite(final_output, paste0(working_dir, "T1DMAPS2_", bfile_base, ".total_score.txt"), sep = "\t")
  log_msg("Final T1D_MAPS2 scores written to file.\n")
}