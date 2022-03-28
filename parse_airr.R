args <- commandArgs(trailingOnly = TRUE)
path_to_constant_calls <- args[1]
path_to_variable_calls <- args[2]
path_to_pre_consensus_gene_calls <- args[3]
number_of_consensus <- args[4]
consensus_directory <- args[5]

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# process the AIRR format files
##### barcode 1 #####
# read in AIRR .tsv
per_read_calls <- read_tsv(path_to_variable_calls)

# read in constant calls from minimap2
constant_calls <- read_tsv(path_to_constant_calls, col_names = FALSE)
constant_calls <- constant_calls[, c(1, 2, 6)] # query name, query length, target name
colnames(constant_calls) <- c("read_name", "read_length", "c_call")
constant_calls <- distinct(constant_calls, read_name, .keep_all = TRUE)

# make the read names the same (and strip a lot of the crap out so its shorter)
constant_calls$read_name <- str_replace(constant_calls$read_name, "_.*", "")
per_read_calls$sequence_id <- str_replace(per_read_calls$sequence_id, "_.*", "")

# combine constant calls with AIRR data
per_read_calls <- left_join(per_read_calls, constant_calls, by = c("sequence_id" = "read_name"))

# fix constant gene name (make it just like IGKC*02 etc)
per_read_calls$c_call <- str_replace(per_read_calls$c_call, ".*?\\|", "")
per_read_calls$c_call <- str_replace(per_read_calls$c_call, "\\|.*", "")

# add n column for the collapsing
per_read_calls$n <- rep(1, nrow(per_read_calls))

# collapse by gene calls 
per_read_calls %>%
  group_by(v_call, d_call, j_call, c_call) %>%
  summarise(count = sum(n), .groups = "keep", 
            reads = paste(sequence_id, collapse = "\n")) -> grouped_clones

# write out the pre-consensus gene call tables
write_tsv(grouped_clones, path_to_pre_consensus_gene_calls)

# write out the names of reads in each group (for subsetting the fasta)
# but only the top 25 bc we're only taking the consensus of those
grouped_clones_H <- filter(grouped_clones, str_detect(v_call, "H"))
grouped_clones_L <- filter(grouped_clones, str_detect(v_call, "H", 
                                                            negate = TRUE))
grouped_clones_H <- arrange(grouped_clones_H, desc(count))
grouped_clones_H <- grouped_clones_H[1:number_of_consensus,]
grouped_clones_H$clone_id <- paste0("H", seq_len(nrow(grouped_clones_H)))

grouped_clones_L <- arrange(grouped_clones_L, desc(count))
grouped_clones_L <- grouped_clones_L[1:number_of_consensus,]
grouped_clones_L$clone_id <- paste0("L", seq_len(nrow(grouped_clones_L)))

# then export for subsetting to fasta
for (i in seq_along(grouped_clones_H$clone_id)) {
  this_clone_id <- grouped_clones_H$clone_id[i]
  write_lines(x = grouped_clones_H$reads[i],
              file = paste0(consensus_directory, this_clone_id, ".txt"))
}

for (i in seq_along(grouped_clones_L$clone_id)) {
  this_clone_id <- grouped_clones_L$clone_id[i]
  write_lines(x = grouped_clones_L$reads[i],
              file = paste0(consensus_directory, this_clone_id, ".txt"))
}


##### choose a starting copy #####
# starting copy is just the longest one from each group

# barcode 1
# making a data frame where one read is one row and it also says which group that
# read is in
grouped_clones_L_sep <- separate_rows(grouped_clones_L, reads, sep = "\n")
grouped_clones_H_sep <- separate_rows(grouped_clones_H, reads, sep = "\n")
groups_keys <- bind_rows(separate_rows(grouped_clones_L, reads, sep = "\n"), separate_rows(grouped_clones_H, reads, sep = "\n"))
groups_keys <- groups_keys[, c("reads", "clone_id")]

# adding this to the per_read_calls
per_read_calls_with_groups <- left_join(per_read_calls, groups_keys, by = c("sequence_id" = "reads"))
per_read_calls_with_groups <- per_read_calls_with_groups[!is.na(per_read_calls_with_groups$clone_id), c("sequence_id", "read_length", "clone_id")]

# select the read for each group that is the longest (most likely to contain a
# full constant region)
per_read_calls_with_groups %>% 
  group_by(clone_id) %>% 
  slice(which.max(read_length)) -> longest_reads

##### write out starting copies #####
for (i in seq_along(longest_reads$sequence_id)) {
  this_clone_id <- longest_reads$clone_id[i]
  write_lines(x = longest_reads$sequence_id[i],
              file = paste0(consensus_directory, this_clone_id, "_starting_point_name.txt"))
}