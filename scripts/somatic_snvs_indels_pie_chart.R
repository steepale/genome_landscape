#===============================================================================
#
#         FILE: /Users/Alec/Documents/Bioinformatics/MDV_Project/genome_landscape/scripts/somatic_snvs_indels_pie_chart.R
#
#  DESCRIPTION: The generation of a pie chart of all the somatic SNVs and INDELs in a example tumor
#                
# REQUIREMENTS:  R
#        NOTES:  ---
#       AUTHOR:  Alec Steep, steepale@msu.edu
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#        			         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2017.05.02
#     REVISION:  
#===============================================================================

# Permanent PROJECT DIRECTORY (Alec's MacBook)
setwd('/Users/Alec/Documents/Bioinformatics/MDV_Project/genome_landscape')


setwd('/Users/Alec/Documents/Bioinformatics/MDV_Project/DNASeq_Analysis/somatic_snp_candidate_determination')

S10 <- read.delim("./data/somatic_snps/918-3_S10_somatic_snp_annotated_1kb1.vcf", header=TRUE)
S10T <- transform(S10, QUAL = as.numeric(QUAL))
library(dplyr)
S10 <- tbl_df(S10T)
NM <- regexpr("NUM_SMMETHODS=[^;]+", S10$INFO, perl=TRUE)
callers <- regmatches(S10$INFO, NM)
callers2 <- sub(pattern="NUM_SMMETHODS=", replacement="", callers)
S10[grep("NUM_SMMETHODS=", S10$INFO), "CALLERS"] <- callers2
S10T <- transform(S10, CALLERS = as.numeric(CALLERS))
AN <- (regexpr("ANN=[^;]+", S10T$INFO, perl=TRUE))
annotations <- regmatches(S10T$INFO, AN)
annotations2 <- strsplit(annotations, "\\|")

S10T$ALT_ALLELE <- sapply(annotations2, function(x) x[1])
S10T$ANNOTATION <- sapply(annotations2, function(x) x[2])
S10T$IMPACT <- sapply(annotations2, function(x) x[3])
S10T$GENE <- sapply(annotations2, function(x) x[4])
S10T$FEATURE_TYPE <- sapply(annotations2, function(x) x[6])
S10T$TRANSCRIPT_BIOTYPE <- sapply(annotations2, function(x) x[8])
S10T$HGVS.c <- sapply(annotations2, function(x) x[10])
S10T$HGVS.p <- sapply(annotations2, function(x) x[11])
S10T$AA_POS_LENGTH <- sapply(annotations2, function(x) x[14])

S10 <- arrange(S10T, CHROM, POS)
S10ft <- filter(S10, CALLERS >=2, QUAL >=30)

Total_mut_S10 <- filter(S10ft, CALLERS >= 3, ANNOTATION != "5_prime_UTR_premature_start_codon_gain_variant", ANNOTATION != "missense_variant&splice_region_variant", ANNOTATION != "non_coding_exon_variant", ANNOTATION != "splice_acceptor_variant&intron_variant", ANNOTATION != "splice_donor_variant&intron_variant", ANNOTATION != "splice_region_variant&intron_variant", SAMPLE == "918-3_S10")
ggplot(Total_mut_S10, aes(x=factor(1), fill = ANNOTATION)) + geom_bar(width = 1) + coord_polar(theta = "y")