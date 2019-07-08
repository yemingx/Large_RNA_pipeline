# Large_RNA_pipeline

Large RNA data analysis from illumina raw fastq to tabular data. 
This pipeline is designed to be deployed on local computing with super user permission.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

System:
Ubuntu or macOS

Software dependencies:
Trimmomatic v0.33
FastQC v0.11.8
HISAT2 v2.1.0
StringTie v1.3.4d
featureCounts v1.6.3

Reference index:
Prebuilt hisat2 index (genome.fa)
GTF file
Trimmomatic adaptor files according to the sequencing platform (TruSeq3-PE.fa)

### Fastq nomination

The fastq files must be renamed according to the following policy.
AA[animal]-FF[fraction]-TT[tissue]-SS[replicate number]_L[lane number]_R[read number]fq

For example
```
AAwt-FFround-TTtestis-SS1_L001_R1.fq

```

AA: animal. AAwt means the sample is wild type.
FF: fraction or function. FFround means round spermatid.
TT: tissue. TTtestis means the sample is the testis tissue.
SS: replicate number.
L: lane number. If the sample is sequenced twice in different lane, "L001" and "L002" should be assigned.
R: read number. For paired-end sequencing, the fastq files should contain R1 and R2. And R1 only for single-end.


The arguments should be numbers or lower case. Arguments after AA, FF, and TT are necessary but not required.

For example
```
AAwt-FFround-TT-SS1_L001_R1.fq

```
This is also legal as long as the user can identify each sample based on the nomination.


### File structure

Files should be organized as described below.

Top level:  folder for all large RNA projects

Secondary: specific large RNA project folder
           hisat2 genome index folder
           reference folder

The renamed fastq files must be put in the folder "raw_fastq" under the specific large RNA project folder.


For example:
```
#specific large RNA project folder
/Users/yemingx/Yan_lab/lRNA_project/spermatogenesis/raw_fastq:
AAko-FFround-TT-SS1_L001_R1.fq AAko-FFround-TT-SS2_L001_R1.fq AAko-FFround-TT-SS3_L001_R1.fq 
AAwt-FFround-TT-SS1_L001_R1.fq AAwt-FFround-TT-SS2_L001_R1.fq AAwt-FFround-TT-SS3_L001_R1.fq 

#hisat2 genome index folder
/Users/yemingx/Yan_lab/lRNA_project/hisat2Genome/mm10:
genes.gtf	genome.fa.1.ht2	genome.fa.3.ht2	genome.fa.5.ht2	genome.fa.7.ht2
genome.fa.2.ht2	genome.fa.4.ht2	genome.fa.6.ht2	genome.fa.8.ht2

#other reference folder
/Users/yemingx/Yan_lab/lRNA_project/reference/trim:
NexteraPE-PE.fa		TruSeq2-SE.fa		TruSeq3-SE.fa
TruSeq.fa		TruSeq3-PE-2.fa		trimmomatic-0.33.jar
TruSeq2-PE.fa		TruSeq3-PE.fa

```

In the example, 
"/Users/yemingx/Yan_lab/lRNA_project/" is the folder for all large RNA projects.
"/Users/yemingx/Yan_lab/lRNA_project/spermatogenesis" is the specific large RNA project folder, which contains the renamed fastq files under the folder named "raw_fastq".
"/Users/yemingx/Yan_lab/lRNA_project/hisat2Genome" is the hisat2 genome index folder, which contains mm10 gtf and prebuilt hisat2 index.
"/Users/yemingx/Yan_lab/lRNA_project/reference" is the reference folder, which contains the trimmomatic adaptor sequences.


### Setting function options

## Usage

1. Get the project folders, genomic index, gtf, and trimmomatic adaptor reference prepared as described in the File structure section.

2. Rename the sample fastq files and copy them to the "raw_fastq" folder, which is under the specfic large RNA project folder.

3. Copy the script under the specific large RNA project folder.
In the example, copy the script to this directory:
```
/Users/yemingx/Yan_lab/lRNA_project/spermatogenesis

```
4. Check the arguments in the script.

5. Make the script executable. Run the pipeline script while generating the log file.

```
sudo chmod 755 ./fastq_hisat2_stringtie_07_08_19.sh
./fastq_hisat2_stringtie_07_08_19.sh 2>&1 | tee make.log

```
## Arguments

Thread number, reference and trimmomatic adaptor sequence should be preset to each project before running the script. The arguments are at the first few lines of the script.

```
thread='8'
GTF="${project_dir}/hisat2Genome/mm10/genes.gtf"
hisat_index="${project_dir}/hisat2Genome/mm10/genome.fa"
trim="${project_dir}/reference/trim/trimmomatic-0.33.jar"
adaptor_SE="${project_dir}/reference/trim/TruSeq3-SE.fa"
adaptor_PE="${project_dir}/reference/trim/TruSeq3-PE.fa"
```

## Output
The script will trim the raw fastq file, generate outputs from fastqc, hisat2, stringtie, and featureCounts in each folder accordingly.
The hisat2 alignment parameters are adjusted for the better sensitivity according the following paper.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5792058/


## Built With

* Bash 3.2.57(1) - macOS Bash used
* Bash 4.4.19(1) - Ubuntu Bash used

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Yeming Xie** - *Chong Tang* - *Tian Yu* - [Yanlab](http://www.weiyanlab.com/home.html)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thanks Dr. Wei Yan lab for the access to multiple RNA-Seq projects that allows this script to be reshaped again and again.
