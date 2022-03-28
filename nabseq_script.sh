#!/bin/bash

while getopts nf:tor:c:h flag
do
    case "${flag}" in
        n) sample_name=${OPTARG};;
        f) reads_fastq=${OPTARG};;
		t) num_threads=${OPTARG};;
		o) out_dir=${OPTARG};;
		r) references=${OPTARG};;
		c) constant_references=${OPTARG};;
		h) how_many_consensus=${OPTARG};;
    esac
done

# default values in case not provided
# sample name
if [[ -z $sample_name ]]
then
sample_name="nabseq_sample"
fi

# output directory
if [[ -z $out_dir ]]
then
out_dir=$(pwd)
fi

# number of threads
if [[ -z $num_threads ]]
then
num_threads=2
fi

# how many consensus sequences to output per chain
if [[ -z $how_many_consensus ]]
then
how_many_consensus=1
fi

# alignment to imgt references
initial_ref_alignment_paf=${out_dir}/${sample_name}_initial_ref_alignment.paf
minimap2 -x map-ont -n $num_threads $references $reads_fastq > $initial_ref_alignment_paf

# subsetting the antibody reads
antibody_read_names=${out_dir}/${sample_name}_antibody_read_names.txt
antibody_reads=${out_dir}/${sample_name}_antibody_reads.fastq
Rscript --vanilla antibody_subset.R --args $initial_ref_alignment_paf $antibody_read_names
seqkit grep -f $antibody_read_names $reads_fastq > $antibody_reads

# trimming of polyA tails
trimmed_antibody_reads=${out_dir}/${sample_name}_trimmed_antibody_reads.fastq
trimmed_antibody_reads_fasta=${out_dir}/${sample_name}_trimmed_antibody_reads.fasta
cutadapt -a "A{100}" -g "T{100}" -n 2 -o $trimmed_antibody_reads $antibody_reads
seqkit fq2fa $trimmed_antibody_reads -o $trimmed_antibody_reads_fasta

# annotating constant regions
constant_calls_pre_consensus=${out_dir}/${sample_name}_pre_consensus_constant_calls.paf
minimap2 --secondary=no -x map-ont -o $constant_calls_pre_consensus $constant_references $trimmed_antibody_reads

# igblast annotation of variable regions
# edit this command as necessary based on your desired igblast parameters
# alternatively, if command line igblast is giving you trouble you can always generate
# an AIRR formatted file using the igblast web server, place it in the output directory
# and name as follows: SAMPLENAME_igblast_pre_consensus.tsv
igblast_pre_consensus_output=${out_dir}/${sample_name}_igblast_pre_consensus.tsv
igblastn -germline_db_V database/rat_V -germline_db_J database/rat_J -germline_db_D database/rat_D \
-organism rat -query $trimmed_antibody_reads_fasta -auxiliary_data optional_file/rat_gl.aux -show_translation \
-num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -outfmt 19 > $igblast_pre_consensus_output


# parsing the AIRR file in R
consensus_directory=${out_dir}/${sample_name}_consensus_dir/
pre_consensus_gene_call_table=${out_dir}/${sample_name}_pre_consensus_gene_calls.tsv
Rscript --vanilla parse_airr.R --args $constant_calls_pre_consensus $igblast_pre_consensus_output \
$pre_consensus_gene_call_table $how_many_consensus $consensus_directory

# making the consensus 
# heavy chain
for i in $( seq 0 $how_many_consensus )
do
	clone_name=${sample_name}_H${i}
	echo $clone_name
	read_name_file=${out_dir}${clone_name}.txt
	starting_point_name_file=${out_dir}${clone_name}_starting_point_name.txt
	clone_reads_file=${out_dir}${clone_name}_reads.fastq
	clone_reads_file_renamed=${out_dir}${clone_name}_reads_renamed.fastq
	clone_reads_file_renamed_dup_fixed=${out_dir}${clone_name}_reads_renamed_dup_fixed.fastq
	clone_reads_file_renamed_dup_fixed_len_fixed=${out_dir}${clone_name}_reads_renamed_dup_fixed_len_fixed.fastq
	polishing_starting_point=${out_dir}${clone_name}_starting_point.fastq
	polishing_starting_point_renamed=${out_dir}${clone_name}_starting_point_renamed.fastq
	all_to_start_overlaps=${out_dir}${clone_name}_overlaps.sam
	racon_consensus=${out_dir}${clone_name}_racon_consensus.fasta
	medaka_consensus_folder=${out_dir}/${clone_name}/

	
	# test if we have a full-length starting point. for now we are only proceeding with full-length starting copies (something to change in future???)
	if [ -f "$starting_point_name_file" ]; then
    	echo "found starting point"
	else 
    	continue
	fi

	# extract reads from fastq
	seqkit grep -n -r -f $read_name_file $trimmed_antibody_reads -o $clone_reads_file
	seqkit replace -p "\_.*" -r "" $clone_reads_file > $clone_reads_file_renamed
	seqkit rename $clone_reads_file_renamed > $clone_reads_file_renamed_dup_fixed
	seqkit seq -m 1 $clone_reads_file_renamed_dup_fixed > $clone_reads_file_renamed_dup_fixed_len_fixed

	# extract starting copy
	seqkit grep -n -r -f $starting_point_name_file $trimmed_antibody_reads -o $polishing_starting_point
	seqkit replace -p ".*" -r $clone_name $polishing_starting_point > $polishing_starting_point_renamed
	
	# grab the first entry to use for polishing !!!!!EDIT 9/9/21: not doing this anymore
	#seqkit head -n 1 $clone_reads_file_renamed > $polishing_starting_point
	#seqkit replace -p ".*" -r $clone_name $polishing_starting_point > $polishing_starting_point_renamed

	# run racon
	minimap2 $polishing_starting_point_renamed $clone_reads_file_renamed_dup_fixed_len_fixed --secondary=no -ax map-ont -o $all_to_start_overlaps
	racon -w 5000 -t 4 -u -g -8 -x -6 -m 8 --no-trimming $clone_reads_file_renamed_dup_fixed_len_fixed $all_to_start_overlaps $polishing_starting_point_renamed > $racon_consensus 

	# run medaka
	mkdir $medaka_consensus_folder
	medaka_consensus -i $clone_reads_file_renamed_dup_fixed_len_fixed -d $racon_consensus -o $medaka_consensus_folder -t 4 -m r941_min_sup_g507
	mv ${medaka_consensus_folder}/consensus.fasta ${out_dir}/${clone_name}_medaka_consensus.fasta
	rm -rf $medaka_consensus_folder
done

# light chain
for i in $( seq 0 $how_many_consensus )
do
	clone_name=${sample_name}_L${i}
	echo $clone_name
	read_name_file=${out_dir}${clone_name}.txt
	starting_point_name_file=${out_dir}${clone_name}_starting_point_name.txt
	clone_reads_file=${out_dir}${clone_name}_reads.fastq
	clone_reads_file_renamed=${out_dir}${clone_name}_reads_renamed.fastq
	clone_reads_file_renamed_dup_fixed=${out_dir}${clone_name}_reads_renamed_dup_fixed.fastq
	clone_reads_file_renamed_dup_fixed_len_fixed=${out_dir}${clone_name}_reads_renamed_dup_fixed_len_fixed.fastq
	polishing_starting_point=${out_dir}${clone_name}_starting_point.fastq
	polishing_starting_point_renamed=${out_dir}${clone_name}_starting_point_renamed.fastq
	all_to_start_overlaps=${out_dir}${clone_name}_overlaps.sam
	racon_consensus=${out_dir}${clone_name}_racon_consensus.fasta
	medaka_consensus_folder=${out_dir}/${clone_name}/

	
	# test if we have a full-length starting point. for now we are only proceeding with full-length starting copies (something to change in future???)
	if [ -f "$starting_point_name_file" ]; then
    	echo "found starting point"
	else 
    	continue
	fi

	# extract reads from fastq
	seqkit grep -n -r -f $read_name_file $trimmed_antibody_reads -o $clone_reads_file
	seqkit replace -p "\_.*" -r "" $clone_reads_file > $clone_reads_file_renamed
	seqkit rename $clone_reads_file_renamed > $clone_reads_file_renamed_dup_fixed
	seqkit seq -m 1 $clone_reads_file_renamed_dup_fixed > $clone_reads_file_renamed_dup_fixed_len_fixed

	# extract starting copy
	seqkit grep -n -r -f $starting_point_name_file $trimmed_antibody_reads -o $polishing_starting_point
	seqkit replace -p ".*" -r $clone_name $polishing_starting_point > $polishing_starting_point_renamed
	
	# grab the first entry to use for polishing !!!!!EDIT 9/9/21: not doing this anymore
	#seqkit head -n 1 $clone_reads_file_renamed > $polishing_starting_point
	#seqkit replace -p ".*" -r $clone_name $polishing_starting_point > $polishing_starting_point_renamed

	# run racon
	minimap2 $polishing_starting_point_renamed $clone_reads_file_renamed_dup_fixed_len_fixed --secondary=no -ax map-ont -o $all_to_start_overlaps
	racon -w 5000 -t 4 -u -g -8 -x -6 -m 8 --no-trimming $clone_reads_file_renamed_dup_fixed_len_fixed $all_to_start_overlaps $polishing_starting_point_renamed > $racon_consensus 

	# run medaka
	mkdir $medaka_consensus_folder
	medaka_consensus -i $clone_reads_file_renamed_dup_fixed_len_fixed -d $racon_consensus -o $medaka_consensus_folder -t 4 -m r941_min_sup_g507
	mv ${medaka_consensus_folder}/consensus.fasta ${out_dir}/${clone_name}_medaka_consensus.fasta
	rm -rf $medaka_consensus_folder
done

# annotation of final corrected consensus sequences
constant_calls_post_consensus=${out_dir}/${sample_name}_post_consensus_constant_calls.paf
minimap2 --secondary=no -x map-ont -o $constant_calls_post_consensus $constant_references $consensus_antibody_reads
igblast_post_consensus_output=${out_dir}/${sample_name}_igblast_post_consensus.tsv
igblastn -germline_db_V database/rat_V -germline_db_J database/rat_J -germline_db_D database/rat_D \
-organism rat -query $consensus_antibody_reads -auxiliary_data optional_file/rat_gl.aux -show_translation \
-num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -outfmt 19 > $igblast_post_consensus_output
