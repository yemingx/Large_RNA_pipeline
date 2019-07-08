#!/bin/bash
set -e
set -u
set -o pipefail
#"""hisat2,stringtie, pipeline"""
#"""author = Yeming 07-08-19"""

thread='8'
GTF="${project_dir}/hisat2Genome/mm10/genes.gtf"
hisat_index="${project_dir}/hisat2Genome/mm10/genome.fa"
trim="${project_dir}/reference/trim/trimmomatic-0.33.jar"
adaptor_SE="${project_dir}/reference/trim/TruSeq3-SE.fa"
adaptor_PE="${project_dir}/reference/trim/TruSeq3-PE.fa"

project_name=$PWD
project_dir=${PWD%/*} 


# create a Bash array from the first col for sample names, because the sample names are the same
# for R1 and R2, we use uniq
cd raw_fastq
seqfile=($(ls AA*_L*))
echo ${seqfile[@]}

##organize samples to creat project, trim everything after _SS?
samplename=(${seqfile[@]%_L*})
samplename=($(printf "%q\n" "${samplename[@]}"|sort -u))
echo ${samplename[@]}
cd ..

##trim illumina adaptor sequence
if test ! -d "./trimmed_fastq"; then
	mkdir ${project_name}/trimmed_fastq
fi

for i in ${samplename[@]}; do
	echo ${i}
	#cut forward adaptor read1
	Rfile=(${i}_L000_R{1,2}.fq)
	echo ${Rfile[@]}
	if [ -f "${project_name}/raw_fastq/${Rfile[0]}" ] && [ ! -f "${project_name}/trimmed_fastq/${Rfile[0]}_trim.fq" ] && [ -f "${project_name}/raw_fastq/${Rfile[1]}" ] && [ ! -f "${project_name}/trimmed_fastq/${Rfile[1]}_trim.fq" ]; then ##test if the input file exists and the output file does not exist
		echo "starting trimming for pe reads for ${i}"
		echo "java -jar ${trim} PE -phred33 ${project_name}/raw_fastq/${Rfile[0]} ${project_name}/raw_fastq/${Rfile[1]} ${project_name}/trimmed_fastq/${Rfile[0]}_trim.fq ${project_name}/trimmed_fastq/${i}_forward_unpaired.fq \
		${project_name}/trimmed_fastq/${Rfile[1]}_trim.fq ${project_name}/trimmed_fastq/${i}_reverse_unpaired.fq ILLUMINACLIP:${adaptor_PE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
		java -jar ${trim} PE -phred33 ${project_name}/raw_fastq/${Rfile[0]} ${project_name}/raw_fastq/${Rfile[1]} ${project_name}/trimmed_fastq/${Rfile[0]}_trim.fq ${project_name}/trimmed_fastq/${i}_forward_unpaired.fq \
		${project_name}/trimmed_fastq/${Rfile[1]}_trim.fq ${project_name}/trimmed_fastq/${i}_reverse_unpaired.fq ILLUMINACLIP:${adaptor_PE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	fi
	if [ -f "${project_name}/raw_fastq/${Rfile[0]}" ] && [ ! -f "${project_name}/raw_fastq/${Rfile[1]}" ] && [ ! -f "${project_name}/trimmed_fastq/${Rfile[0]}_trim.fq" ]; then
		echo "startig trimming for se reads for ${i}"
		echo "java -jar ${trim} SE -phred33 ${project_name}/raw_fastq/${Rfile[0]} ${project_name}/trimmed_fastq/${Rfile[0]}_trim.fq ILLUMINACLIP:${adaptor_SE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
		java -jar ${trim} SE -phred33 ${project_name}/raw_fastq/${Rfile[0]} ${project_name}/trimmed_fastq/${Rfile[0]}_trim.fq ILLUMINACLIP:${adaptor_SE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	fi


	T_file=(${i}_L000_R{1,2}.fq_trim.fq)
	echo ${T_file[@]}
	#fastqc#
	if test ! -d "./fastqc_out"; then
		mkdir ./fastqc_out
	fi
	if [ ! -f "./fastqc_out/${i}_L000_R1.fq_trim_fastqc.html" ] && [ -f "${project_name}/trimmed_fastq/${T_file[1]}" ] && [ -f "${project_name}/trimmed_fastq/${T_file[0]}" ]; then
		echo "fastqc ${project_name}/trimmed_fastq/${T_file[0]} -o ./fastqc_out"	
		fastqc ${project_name}/trimmed_fastq/${T_file[0]} -o ./fastqc_out
	fi
	if [ ! -f "./fastqc_out/${i}_L000_R2.fq_trim_fastqc.html" ] && [ -f "${project_name}/trimmed_fastq/${T_file[1]}" ] && [ -f "${project_name}/trimmed_fastq/${T_file[0]}" ]; then
		echo "fastqc ${project_name}/trimmed_fastq/${T_file[1]} -o ./fastqc_out"
		fastqc ${project_name}/trimmed_fastq/${T_file[1]} -o ./fastqc_out
	fi
	if  [ ! -f "./fastqc_out/${i}_L000_R1.fq_trim_fastqc.html" ] && [ ! -f "${project_name}/trimmed_fastq/${T_file[1]}" ] && [ -f "${project_name}/trimmed_fastq/${T_file[0]}" ]; then
		echo "fastqc ${project_name}/trimmed_fastq/${T_file[0]} -o ./fastqc_out"
		fastqc ${project_name}/trimmed_fastq/${T_file[0]} -o ./fastqc_out
	fi
	####Hisat2 alignemnt and converting to bam output####	
	if test ! -d "./hisat2_out"; then
		mkdir ./hisat2_out
	fi
	if [ ! -f "./hisat2_out/${i}.bam" ] && [ -f "${project_name}/trimmed_fastq/${T_file[1]}" ] && [ -f "${project_name}/trimmed_fastq/${T_file[0]}" ]; then
		echo "starting hisat alignment for ${i} PE"
		echo "hisat2 -p ${thread} --dta -x $hisat_index -1 ${project_name}/trimmed_fastq/${T_file[0]} -2 ${project_name}/trimmed_fastq/${T_file[1]} 2> ./hisat2_out/${i}_hisat_alignment.txt | samtools sort -@ ${thread} -o ./hisat2_out/${i}.bam"
		hisat2 -p ${thread} -N 1 -L 20 -i S,1,0.5 -D 25 -R 5 --pen-noncansplice 12 --mp 1,0 --sp 3,0 --dta -x $hisat_index -1 ${project_name}/trimmed_fastq/${T_file[0]} -2 ${project_name}/trimmed_fastq/${T_file[1]} 2> ./hisat2_out/${i}_hisat_alignment.txt | samtools sort -@ ${thread} -o ./hisat2_out/${i}.bam
	fi

	if [ ! -f "./hisat2_out/${i}.bam" ] && [ ! -f "${project_name}/trimmed_fastq/${T_file[1]}" ] && [ -f "${project_name}/trimmed_fastq/${T_file[0]}" ]; then
		echo "starting hisat alignment for ${i} SE"
		echo "hisat2 -p ${thread} --dta -x $hisat_index -U ${project_name}/trimmed_fastq/${T_file[0]} 2> ./hisat2_out/${i}_hisat_alignment.txt| samtools sort -@ ${thread} -o ./hisat2_out/${i}.bam"
		hisat2 -p ${thread} -N 1 -L 20 -i S,1,0.5 -D 25 -R 5 --pen-noncansplice 12 --mp 1,0 --sp 3,0 --dta -x $hisat_index -U ${project_name}/trimmed_fastq/${T_file[0]} 2> ./hisat2_out/${i}_hisat_alignment.txt | samtools sort -@ ${thread} -o ./hisat2_out/${i}.bam
	fi
	###Stringtie for making gtf files
 	if test ! -d "./stringtie_out"; then
 		mkdir ./stringtie_out
 	fi
 	if [ ! -f "./stringtie_out/${i}.gtf" ]; then
 		echo "stringtie -p ${thread} -G $GTF -o ./stringtie_out/${i}.gtf -l ${i} ./hisat2_out/${i}.bam"
		stringtie -p ${thread} -G $GTF -o ./stringtie_out/${i}.gtf -l ${i} ./hisat2_out/${i}.bam
	fi
	t=($(find ./stringtie_out/${i}.gtf))
	echo ${t[@]} >> ./stringtie_out/merge_list.txt 
done

###featureCounts###
if test ! -d './featureCounts_out'; then
	mkdir ./featureCounts_out
fi
g=($(find ./hisat2_out/*.bam))
echo ${g[@]}	
if [ ! -f "./featureCounts_out/gene_counts.txt" ]; then
	echo "featureCounts -a $GTF -g gene_id -o ./featureCounts_out/gene_counts.txt ${g[@]}"
	featureCounts -a $GTF -g gene_id -o ./featureCounts_out/gene_counts.txt ${g[@]}
fi
if [ ! -f "./featureCounts_out/transcript_counts.txt" ]; then
	echo "featureCounts -a $GTF -g transcript_id -o ./featureCounts_out/transcript_counts.txt ${g[@]}"
	featureCounts -a $GTF -g transcript_id -o ./featureCounts_out/transcript_counts.txt ${g[@]}
fi

###Merge transcripts from all samples###
if [ ! -f "./stringtie_out/stringtie_merged.gtf" ]; then
	echo "stringtie --merge -p ${thread} -G $GTF -o ./stringtie_out/stringtie_merged.gtf ./stringtie_out/merge_list.txt"
	stringtie --merge -p ${thread} -G $GTF -o ./stringtie_out/stringtie_merged.gtf ./stringtie_out/merge_list.txt 
fi

###Compare the stringtie assemblied transcript with the refrence annotation###
if test ! -d "./gffcompare"
		then mkdir ./gffcompare
		echo "start gff compare"
		gffcompare -r $GTF -o ./gffcompare/merged ./stringtie_out/stringtie_merged.gtf
fi