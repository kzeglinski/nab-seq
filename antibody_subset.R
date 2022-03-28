args <- commandArgs(trailingOnly = TRUE)
path_to_paf_file <- args[1]
output_read_names <- args[2]

# read in the paf file
ab_alignment_paf <- read_tsv(path_to_paf_file, col_names = FALSE)
ab_alignment_paf <- ab_alignment_paf[, c(1, 6)] # only need the 1st and 6th columns

# how many heavy chain transcripts?
# btw we need to use unique() because one Ab read can align to multiple
# references: the heavy constant region has CH1, CH2, CH3 etc and all of these
# are separate reference sequences in IMGT. Using unique() means we don't count
# the same read twice
length(unique(ab_alignment_paf[grepl("IGH", ab_alignment_paf$X6) ,]$X1))

# heavy chain isotypes?
length(unique(ab_alignment_paf[grepl("IGHD", ab_alignment_paf$X6) ,]$X1))
length(unique(ab_alignment_paf[grepl("IGHM", ab_alignment_paf$X6) ,]$X1))

# what about light chain?
# kappa
length(unique(ab_alignment_paf[grepl("IGK", ab_alignment_paf$X6) ,]$X1))

# lambda
length(unique(ab_alignment_paf[grepl("IGL", ab_alignment_paf$X6) ,]$X1))

# make a list of all the Ab reads
# using the paf file previously read in
all_ab_reads <- as.data.frame(unique(ab_alignment_paf$X1))
write_tsv(all_ab_reads, output_read_names, col_names = FALSE)