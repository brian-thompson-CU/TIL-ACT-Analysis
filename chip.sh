#! /bin/bash
# Author: Brian Thompson 
# Date: August 11 2022
# About: One complete code for processing the TIL Epigenome ChIP Data

#Setting the number of threads to use
THREADS=8

#Setting root folder
ROOT_FOLDER="/Volumes/2TB_SSD/TIL_Epigenome_Chip"
FASTQ_FOLDER_Group1="/Volumes/2TB_SSD/TIL_Epigenome_Chip/group_1_fastq"
FASTQ_FOLDER_Group2="/Volumes/2TB_SSD/TIL_Epigenome_Chip/group_2_fastq"
INDEX_FILE="/Volumes/2TB_SSD/TIL_Epigenome_Chip/Index/GRCh38_noalt_as"

#Creating log file
touch ${ROOT_FOLDER}/log.txt
echo "Analysis started on $(date +%y-%m-%d\ %H:%M:%S)" "\n" >> ${ROOT_FOLDER}/log.txt
echo Root folder location: $ROOT_FOLDER >> ${ROOT_FOLDER}/log.txt
echo Fastq folder location: $FASTQ_FOLDER_Group1 >> ${ROOT_FOLDER}/log.txt
echo Fastq folder location: $FASTQ_FOLDER_Group2 >> ${ROOT_FOLDER}/log.txt
echo Index file location: $INDEX_FILE >> ${ROOT_FOLDER}/log.txt
echo Genome used as index: "GRCh38_noalt" >> ${ROOT_FOLDER}/log.txt


#Zipping files to save space, if they are not already compressed
for f in $(ls $FASTQ_FOLDER_Group1/*.fastq); do
pigz -p ${THREADS} $f 
done

mkdir ${ROOT_FOLDER}/fastqc
for f in $(ls $FASTQ_FOLDER_Group1/*.gz); do 
fastqc -o ${ROOT_FOLDER}/fastqc -t ${THREADS} $f
done

for f in $(ls $FASTQ_FOLDER_Group2/*.gz); do 
fastqc -o ${ROOT_FOLDER}/fastqc -t ${THREADS} $f
done

#Trimming reads, although they all look to be of pretty good quality
#Reads are trimmed to remove adatper sequence and low quality bases\
#that may impair alignment
#First make directory for trimmed files, keep groups separate
mkdir ${ROOT_FOLDER}/trimmed_group1
for f in $(ls $FASTQ_FOLDER_Group1/*.gz); do
trim_galore --fastqc --gzip --quality 20 --cores 8 --output_dir ${ROOT_FOLDER}/trimmed_group1 $f
mkdir ${ROOT_FOLDER}/group1_trim_results
mv ${ROOT_FOLDER}/trimmed_group1/*.txt ${ROOT_FOLDER}/group1_trim_results
done

#Trimming group 2
mkdir ${ROOT_FOLDER}/trimmed_group2
for f in $(ls $FASTQ_FOLDER_Group2/*.gz); do
trim_galore --fastqc --gzip --quality 20 --cores 8 --output_dir ${ROOT_FOLDER}/trimmed_group2 $f
mkdir ${ROOT_FOLDER}/group2_trim_results
mv ${ROOT_FOLDER}/trimmed_group2/*.txt ${ROOT_FOLDER}/group2_trim_results
done

#Now I need to align the reads to the human genome with bowtie2
mkdir ${ROOT_FOLDER}/aligned_group1
mkdir ${ROOT_FOLDER}/aligned_group1_stats
echo $ROOT_FOLDER
for f in $(ls ${ROOT_FOLDER}/trimmed_group1/*.fq.gz); do
echo $f
 sam_file="${f/$ROOT_FOLDER/}"
  echo $sam_file
 sam_file="${sam_file/*_group1/}"
 echo $sam_file
  sam_file="${sam_file///}"
 echo $sam_file
 sam_file="${sam_file/.fq.gz/.sam}"
 echo $sam_file
 touch ${ROOT_FOLDER}/aligned_group1_stats/$sam_file.txt
BOWTIE="(bowtie2 -p $THREADS -x $INDEX_FILE -U $f -S ${ROOT_FOLDER}/aligned_group1/$sam_file) 2> ${ROOT_FOLDER}/aligned_group1_stats/$sam_file.txt"
    echo $BOWTIE 
    echo $BOWTIE >> ${ROOT_FOLDER}/log.txt
    eval $BOWTIE
done


#Aligning group 2
mkdir ${ROOT_FOLDER}/aligned_group2
mkdir ${ROOT_FOLDER}/aligned_group2_stats
echo $ROOT_FOLDER
for f in $(ls ${ROOT_FOLDER}/trimmed_group2/*.fq.gz); do
echo $f
 sam_file2="${f/$ROOT_FOLDER/}"
  echo $sam_file2
 sam_file2="${sam_file2/*_group2/}" 
 echo $sam_file2
  sam_file2="${sam_file2///}"
 echo $sam_file2
 sam_file2="${sam_file2/.fq.gz/.sam}" 
 echo "sam file: $sam_file2"
 touch ${ROOT_FOLDER}/aligned_group2_stats/$sam_file.txt
BOWTIE="(bowtie2 -p $THREADS -x $INDEX_FILE -U $f -S ${ROOT_FOLDER}/aligned_group2/$sam_file2) >& ${ROOT_FOLDER}/aligned_group2_stats/$sam_file2.txt" 
    echo $BOWTIE 
    echo $BOWTIE >> ${ROOT_FOLDER}/log.txt
    eval $BOWTIE #could have put stderr redirect here
done

#Converting sam to bam
for f in $(ls ${ROOT_FOLDER}/aligned_group1/*.sam); do
echo $f
bam_file="${f/$ROOT_FOLDER/}"
echo $bam_file
bam_file="${bam_file/*_group1/}"
echo $bam_file
bam_file="${bam_file/.sam/.bam}" 
 echo $bam_file
 bam_file="${bam_file///}"
 echo $bam_file
samtools view -S -b $f > ${ROOT_FOLDER}/aligned_group1/$bam_file
done

#Cleaning bam reads
#The bam reads must be cleaned before peak calling
#First check line counts
mkdir ${ROOT_FOLDER}/raw_line_counts_group1
for f in $(ls ${ROOT_FOLDER}/aligned_group1/*.sam); do
echo $f
bam_file="${f/$ROOT_FOLDER/}"
echo $bam_file
bam_file="${bam_file/*_group1/}" 
 echo $bam_file
 bam_file="${bam_file///}"
 echo $bam_file
samtools view -c $f > ${ROOT_FOLDER}/raw_line_counts_group1/$bam_file.txt
done

#First sort the bam files with samtools
mkdir ${ROOT_FOLDER}/sorted_group1
for f in $(ls ${ROOT_FOLDER}/aligned_group1/*.bam); do
echo $f
bam_file="${f/$ROOT_FOLDER/}"
echo $bam_file
bam_file="${bam_file/*_group1/}" 
 echo $bam_file
 bam_file="${bam_file///}"
 echo $bam_file
Samtools sort -@ ${THREADS} $f -o ${ROOT_FOLDER}/sorted_group1/$bam_file.sorted.bam
done

#Index the sorted files
#for f in $(ls ${ROOT_FOLDER}/sorted_group1/*.bam); do
echo $f
bam_file="${f/$ROOT_FOLDER/}"
echo $bam_file
bam_file="${bam_file/*_group1/}" 
echo $bam_file
 bam_file="${bam_file///}"
 echo $bam_file
Samtools index -@ ${THREADS} $f 
done

#Removing multimapping reads
mkdir ${ROOT_FOLDER}/filtered_group1
for f in $(ls ${ROOT_FOLDER}/sorted_group1/*.bam); do
echo $f
bam_file="${f/$ROOT_FOLDER/}"
echo $bam_file
bam_file="${bam_file/*_group1/}" 
 echo $bam_file
 bam_file="${bam_file///}"
 echo $bam_file
Samtools view $f -@ ${THREADS} -q 20 -o ${ROOT_FOLDER}/filtered_group1/$bam_file.filtered.bam
done

#Removing duplicate reads
mkdir ${ROOT_FOLDER}/dedup_group1
for f in $(ls ${ROOT_FOLDER}/filtered_group1/*.bam); do
echo $f
bam_file="${f/$ROOT_FOLDER/}"
echo $bam_file
bam_file="${bam_file/*_group1/}" 
 echo $bam_file
 bam_file="${bam_file///}"
 echo $bam_file
picard MarkDuplicates INPUT=$f OUTPUT=${ROOT_FOLDER}/dedup_group1/$bam_file.marked.bam METRICS_FILE=${ROOT_FOLDER}/dedup_group1/$bam_file.sorted.metrics \
         REMOVE_DUPLICATES=TRUE CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
done

#Index the sorted files
for f in $(ls ${ROOT_FOLDER}/dedup_group1/*.bam); do
echo $f
bam_file="${f/$ROOT_FOLDER/}"
echo $bam_file
bam_file="${bam_file/*_group1/}" 
echo $bam_file
 bam_file="${bam_file///}"
 echo $bam_file
Samtools index -@ ${THREADS} $f 
done

#Final line count 
mkdir ${ROOT_FOLDER}/final_line_counts_group1
for f in $(ls ${ROOT_FOLDER}/dedup_group1/*.bam); do
echo $f
bam_file="${f/$ROOT_FOLDER/}"
echo $bam_file
bam_file="${bam_file/*_group1/}" 
 echo $bam_file
 bam_file="${bam_file///}"
 echo $bam_file
samtools view -c $f > ${ROOT_FOLDER}/final_line_counts_group1/$bam_file.txt
done

#Calling peaks 
mkdir ${ROOT_FOLDER}/peaks_group1

#Input file is the ChIP Control file
GROUP1_INPUT=${ROOT_FOLDER}/dedup_group1/13*.bam
echo Group 1 Input file is: $GROUP1_INPUT >> ${ROOT_FOLDER}/log.txt
for f in $(ls in ${ROOT_FOLDER}/dedup_group1/*.bam); do
bam_file="${f/$ROOT_FOLDER/}"
echo $bam_file
bam_file="${bam_file/*_group1/}" 
 echo $bam_file
 bam_file="${bam_file///}"
 echo $bam_file
macs2 callpeak -t $f -c ${GROUP1_INPUT} -f BAM -g 2.7e+9 --broad -n $bam_file --outdir ${ROOT_FOLDER}/peaks_group1
done 

#####Repeating cleanup and later steps for group 2
#Converting sam to bam
for f in $(ls ${ROOT_FOLDER}/aligned_group2/*.sam); do
echo $f
bam_file2="${f/$ROOT_FOLDER/}"
echo $bam_file2
bam_file2="${bam_file2/*_group2/}"
echo $bam_file2
bam_file2="${bam_file2/.sam/.bam}" 
 echo $bam_file2
 bam_file2="${bam_file2///}"
 echo $bam_file2
samtools view -S -b $f > ${ROOT_FOLDER}/aligned_group2/$bam_file2
done

#Cleaning bam reads
#The bam reads must be cleaned before peak calling
#First check line counts
mkdir ${ROOT_FOLDER}/raw_line_counts_group2
for f in $(ls ${ROOT_FOLDER}/aligned_group2/*.bam); do
echo $f
bam_file2="${f/$ROOT_FOLDER/}"
echo $bam_file2
bam_file2="${bam_file2/*_group2/}" 
 echo $bam_file2
 bam_file2="${bam_file2///}"
 echo $bam_file2
samtools view -c $f > ${ROOT_FOLDER}/raw_line_counts_group2/$bam_file2.txt
done

#First sort the bam files with samtools
mkdir ${ROOT_FOLDER}/sorted_group2
for f in $(ls ${ROOT_FOLDER}/aligned_group2/*.bam); do
echo $f
bam_file2="${f/$ROOT_FOLDER/}"
echo $bam_file2
bam_file2="${bam_file2/*_group2/}" 
 echo $bam_file2
 bam_file2="${bam_file2///}"
 echo $bam_file2
Samtools sort -@ ${THREADS} $f -o ${ROOT_FOLDER}/sorted_group2/$bam_file2.sorted.bam
done

#Index the sorted files
for f in $(ls ${ROOT_FOLDER}/sorted_group2/*.bam); do
echo $f
bam_file2="${f/$ROOT_FOLDER/}"
echo $bam_file2
bam_file2="${bam_file2/*_group2/}" 
 echo $bam_file2
 bam_file2="${bam_file2///}"
 echo $bam_file2
Samtools index -@ ${THREADS} $f 
done

#Removing multimapping reads
mkdir ${ROOT_FOLDER}/filtered_group2
for f in $(ls ${ROOT_FOLDER}/sorted_group2/*.bam); do
echo $f
bam_file2="${f/$ROOT_FOLDER/}"
echo $bam_file2
bam_file2="${bam_file2/*_group2/}" 
 echo $bam_file2
 bam_file2="${bam_file2///}"
 echo $bam_file2
Samtools view $f -@ ${THREADS} -q 20 -o ${ROOT_FOLDER}/filtered_group2/$bam_file2.filtered.bam
done

mkdir ${ROOT_FOLDER}/dedup_group2
for f in $(ls ${ROOT_FOLDER}/filtered_group2/*.bam); do
echo $f
bam_file2="${f/$ROOT_FOLDER/}"
echo $bam_file2
bam_file2="${bam_file2/*_group2/}" 
 echo $bam_file2
 bam_file2="${bam_file2///}"
 echo $bam_file2
picard MarkDuplicates INPUT=$f OUTPUT=${ROOT_FOLDER}/dedup_group2/$bam_file2.marked.bam METRICS_FILE=${ROOT_FOLDER}/dedup_group2/$bam_file2.sorted.metrics \
         REMOVE_DUPLICATES=TRUE CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
done

#Index the sorted files
for f in $(ls ${ROOT_FOLDER}/dedup_group2/*.bam); do
echo $f
bam_file2="${f/$ROOT_FOLDER/}"
echo $bam_file2
bam_file2="${bam_file2/*_group2/}" 
 echo $bam_file2
 bam_file2="${bam_file2///}"
 echo $bam_file2
Samtools index -@ ${THREADS} $f 
done

Final line count 
mkdir ${ROOT_FOLDER}/final_line_counts_group2
for f in $(ls ${ROOT_FOLDER}/dedup_group2/*.bam); do
echo $f
bam_file2="${f/$ROOT_FOLDER/}"
echo $bam_file2
bam_file2="${bam_file2/*_group2/}" 
 echo $bam_file2
 bam_file2="${bam_file2///}"
 echo $bam_file2
samtools view -c $f > ${ROOT_FOLDER}/final_line_counts_group2/$bam_file2.txt
done

#Calling peaks 
mkdir ${ROOT_FOLDER}/peaks_group2

#Input file is the ChIP Control file
GROUP2_INPUT=${ROOT_FOLDER}/dedup_group2/17*.bam
echo Group 2 Input file is: $GROUP2_INPUT >> ${ROOT_FOLDER}/log.txt
for f in $(ls in ${ROOT_FOLDER}/dedup_group2/*.bam); do
echo $f
bam_file2="${f/$ROOT_FOLDER/}"
echo $bam_file2
bam_file2="${bam_file2/*_group2/}" 
 echo $bam_file2
 bam_file2="${bam_file2///}"
 echo $bam_file2
macs2 callpeak -t $f -c ${GROUP2_INPUT} -f BAM -g 2.7e+9 --broad -n $bam_file2 --outdir ${ROOT_FOLDER}/peaks_group2
done 




