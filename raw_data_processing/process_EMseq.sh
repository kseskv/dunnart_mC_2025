#!/bin/bash

#activate mamba
eval "$(conda shell.bash hook)"
source /home/913/aa3618/miniforge3/etc/profile.d/mamba.sh #replace with your own mamba path

###--- run fastqc (pre-trimming)
module load fastqc/0.11.7

INPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/fastq/ANG15762/"
OUTPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/fastq/ANG15762/fastqc/"
mkdir -p ${OUTPUT_DIR}
cd ${INPUT_DIR}

for read in *.fastq*; do
fastqc -o ${OUTPUT_DIR} -f fastq $read
done

###--- trim reads (paired end)
mamba activate trimmomatic-0.39

FASTQ_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/fastq/ANG15762/"
ADAPTER_DIR="/home/913/aa3618/miniforge3/envs/trimmomatic-0.39/share/trimmomatic-0.39-2/adapters/"
cd ${FASTQ_DIR}

for read_1 in *_R1_001.fastq.gz; do #check the name of fastq file read_1
echo $read_1
read_2=$(echo $read_1 | sed s/R1/R2/) # check the read_1 name
echo $read_2
ext_paired='_paired.fastq'
ext_unpaired='_unpaired.fastq'
output_paired_1=$(echo $read_1 | rev | cut -c 10- | rev)$ext_paired
output_unpaired_1=$(echo $read_1 | rev | cut -c 10- | rev)$ext_unpaired
output_paired_2=$(echo $read_2 | rev | cut -c 10- | rev)$ext_paired
output_unpaired_2=$(echo $read_2 | rev | cut -c 10- | rev)$ext_unpaired
CMD="trimmomatic PE -threads 24 $read_1 $read_2 $output_paired_1 $output_unpaired_1 $output_paired_2 $output_unpaired_2 \
ILLUMINACLIP:${ADAPTER_DIR}TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50"
echo $CMD && eval $CMD
done

###--- run fastqc (post-trimming)
INPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/fastq/ANG15762/"
OUTPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/fastq/ANG15762/fastqc/"
mkdir -p ${OUTPUT_DIR}
cd ${INPUT_DIR}

for read in *_paired.fastq; do
fastqc -o ${OUTPUT_DIR} -f fastq $read
done

###--- index genome
WALT="/g/data/vn68/aa3618/tools/walt/bin/"
FASTA_DIR="/g/data/vn68/aa3618/genomes/Sminthopsis_crassicaudata/T2T_ONT_genome/"
INDEX_DIR="/g/data/vn68/aa3618/genomes/Sminthopsis_crassicaudata/T2T_ONT_genome/WALTIndex_lambda_pUC19/"
mkdir -p ${INDEX_DIR}

${WALT}makedb -c ${FASTA_DIR}"dunnart_male_T2T_21082025_lambda_pUC19.fa" -o ${INDEX_DIR}"genome_lambda_pUC19.dbindex"

###--- map reads
OUTPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/walt/"
mkdir -p ${OUTPUT_DIR}
cd ${FASTQ_DIR}

for read_1 in *_R1*_001_paired.fastq; do
echo $read_1
read_2=$(echo $read_1 | sed 's/_R1_/_R2_/')
echo $read_2
ext='.sam'
output_sam=$(echo $read_1 | rev | cut -c 7- | rev)$ext # Removes n last characters from the string to prepare output file name
echo $output_sam
CMD="${WALT}walt -i ${INDEX_DIR}"genome_lambda_pUC19.dbindex" \
-1 $read_1 -2 $read_2 \
-sam -m 10 -o ${OUTPUT_DIR}$output_sam \
-t 24 -N 10000000 -L 2000"
echo $CMD && eval $CMD
done

###--- convert sam to bam
module load samtools/1.19
INPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/walt/"
cd ${INPUT_DIR}

for input_sam in *.sam; do
echo $input_sam
ext_bam='.bam'
output_bam=$(echo $input_sam | rev | cut -c 5- | rev)$ext_bam
echo $output_bam
CMD="samtools view -Sb $input_sam > $output_bam"
echo $CMD && eval $CMD
done

###--- sort bam
for input_bam in *.bam; do
echo $input_bam
ext_sorted_bam='_sorted.bam'
output_sorted_bam=$(echo $input_bam | rev | cut -c 5- | rev)$ext_sorted_bam 
echo $output_sorted_bam
CMD="samtools sort -o $output_sorted_bam -@ 4 $input_bam;
samtools index $output_sorted_bam"
echo $CMD && eval $CMD
done

###--- deduplicate bam files using picard
mamba activate picard-3.1.1
module load sambamba/0.8.1

for input_bam in *_sorted.bam; do
echo $input_bam
ext_dedup_bam='_dedup.bam'
output_dedup_bam=$(echo $input_bam | rev | cut -c 5- | rev)$ext_dedup_bam 
echo $output_dedup_bam
ext_dup_txt='_marked_dup_metrics.txt'
output_dup_txt=$(echo $input_bam | rev | cut -c 5- | rev)$ext_dup_txt 
echo $output_dup_txt
CMD="picard MarkDuplicates REMOVE_DUPLICATES=true \
I=$input_bam \
O=$output_dedup_bam \
M=$output_dup_txt"
echo $CMD && eval $CMD
done

###--- sort bam
for input_bam in *_dedup.bam; do
echo $input_bam
ext_sorted_bam='_sorted.bam'
output_sorted_bam=$(echo $input_bam | rev | cut -c 5- | rev)$ext_sorted_bam 
echo $output_sorted_bam
CMD="samtools sort -o $output_sorted_bam -@ 4 $input_bam;
samtools index $output_sorted_bam"
echo $CMD && eval $CMD
done

###--- identify reads with >=3 non-converted cytosines
module load samtools/1.19

OUTPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/walt/"
cd ${OUTPUT_DIR}

for input_bam in *_dedup_sorted.bam; do
echo $input_bam
output_filter_list="${input_bam%.bam}_filter_list"
echo $output_filter_list
CMD="samtools view -f2 $input_bam | \
awk '{if (\$2==163 || \$2==83) {print \$1, gsub(/TG|AG|GG/,\"\",\$10)} \
else if (\$2==147 || \$2==99) {print \$1, gsub(/CC|CA|CT/,\"\",\$10)}}' | \
awk '{if (\$2>3) {print \$1}}' > $output_filter_list"
echo $CMD && eval $CMD
done

###--- filter reads from bam file
mamba activate picard-3.1.1
OUTPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/walt/"
cd ${OUTPUT_DIR}

###--- filter reads from bam file
mamba activate picard-3.1.1

for input_bam in *_dedup_sorted.bam; do
echo $input_bam
ext='_filter.bam'
ext2='_filter_list'
output_bam=$(echo $input_bam | rev | cut -c 5- | rev)$ext 
filter_list=$(echo $input_bam | rev | cut -c 5- | rev)$ext2
echo $output_bam
echo $filter_list
CMD="picard -Xmx100g FilterSamReads I=$input_bam O=$output_bam READ_List_File=$filter_list Filter=excludeReadList"
echo $CMD && eval $CMD
done

###--- sort bam
module load samtools/1.19

for input_bam in *filter.bam; do
echo $input_bam
ext_sorted_bam='_sorted.bam'
output_sorted_bam=$(echo $input_bam | rev | cut -c 5- | rev)$ext_sorted_bam 
echo $output_sorted_bam
CMD="samtools sort -o $output_sorted_bam -@ 4 $input_bam;
samtools index $output_sorted_bam"
echo $CMD && eval $CMD
done


###--- MethylDackel mbias
mamba activate /g/data/ih05/aa3618/miniforge3/envs/methyldackel-0.6.1

FASTA_DIR="/g/data/vn68/aa3618/genomes/Sminthopsis_crassicaudata/T2T_ONT_genome/"
INPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/walt/"
OUTPUT_DIR="/g/data/vn68/aa3618/dunnart_reprogramming/EMseq/methyldackel/"
mkdir -p ${OUTPUT_DIR}

for input_bam in *_dedup_sorted_filter_sorted.bam; do
echo $input_bam
ext='.svg'
output=$(echo $input_bam | rev | cut -c 5- | rev)$ext
echo $output
CMD="MethylDackel mbias \
${FASTA_DIR}"dunnart_male_T2T_21082025_lambda_pUC19.fa" \
$input_bam \
${OUTPUT_DIR}$output"
echo $CMD && eval $CMD
done

###--- MethylDackel CpG methylation calling 
for input_bam in *_dedup_sorted_filter_sorted.bam; do
echo $input_bam
output=$(echo $input_bam | rev | cut -c 32- | rev)
echo $output
CMD="MethylDackel extract \
${FASTA_DIR}"dunnart_male_T2T_21082025_lambda_pUC19.fa" \
$input_bam \
-o ${OUTPUT_DIR}$output \
--mergeContext \
--minOppositeDepth 5 --maxVariantFrac 0.5 \
--OT 20,130,20,130 --OB 20,130,20,130"
echo $CMD && eval $CMD
done
