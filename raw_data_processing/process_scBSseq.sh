#!/bin/bash

# forloops are written so that the command only runs for input files when an output file doesn't already exist
# this is necessary because of the time constraints of jobs on NCI
# if you require >48 hours, make sure that you don't have any incomplete outputs before rerunning this job

###--- trim reads
TRIMGALORE="/g/data/vn68/aa3618/tools/TrimGalore-0.6.10/"
CUTADAPT="/g/data/vn68/aa3618/tools/cutadapt/bin/cutadapt"
export PYTHONPATH=/g/data/vn68/aa3618/tools/cutadapt/

INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/fastq/"
OUTPUT_DIR="${INPUT_DIR}/trimmed/"
mkdir -p "${OUTPUT_DIR}"
cd "${INPUT_DIR}"

MAX_JOBS=5
JOB_COUNT=0

for read_1 in *_R1_001.fastq.gz; do
    echo "Processing $read_1"
    read_2=$(echo "$read_1" | sed 's/_R1_/_R2_/')
    echo "Pair: $read_2"

    sample=$(basename "$read_1" | sed 's/_R1_001.fastq.gz//')
    output_1="${OUTPUT_DIR}${sample}_R1_001_val_1.fq.gz"
    output_2="${OUTPUT_DIR}${sample}_R2_001_val_2.fq.gz"

    if [[ -f "$output_1" && -f "$output_2" ]]; then
        echo "Output files for ${sample} exist, skipping..."
        continue
    fi

    (
    echo "Starting trimming for ${sample}..."
    ${TRIMGALORE}trim_galore --cores 4 \
        --path_to_cutadapt ${CUTADAPT} \
        --paired --output_dir ${OUTPUT_DIR} \
        --clip_r1 6 --clip_r2 6 \
        "$read_1" "$read_2"
    echo "Finished trimming for ${sample}"
    ) &

    ((JOB_COUNT++))
    if [[ ${JOB_COUNT} -ge ${MAX_JOBS} ]]; then
        wait
        JOB_COUNT=0
    fi
done

wait
echo "All trimming jobs complete."

###--- fastqc
module load fastqc/0.12.1

# pre trimming
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/fastq/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/fastq/fastqc/"
mkdir -p ${OUTPUT_DIR}
cd ${INPUT_DIR}

for read in *.fastq*; do
fastqc -o ${OUTPUT_DIR} -f fastq $read
done

# post trimming
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/fastq/trimmed/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/fastq/fastqc/trimmed/"
mkdir -p ${OUTPUT_DIR}
cd ${INPUT_DIR}

for read in *.fq.gz; do
fastqc -o ${OUTPUT_DIR} -f fastq $read
done

###--- bismark genome indexing
module load bowtie2/2.3.5.1 
module load samtools/1.22 
BISMARK="/g/data/vn68/aa3618/tools/Bismark-0.25.1/"
GENOME="/g/data/vn68/aa3618/genomes/Sminthopsis_crassicaudata/T2T_ONT_genome/lambda/"

# run bismark indexing
${BISMARK}bismark_genome_preparation ${GENOME}

###--- map reads with bismark
GENOME_DIR="/g/data/vn68/aa3618/genomes/Sminthopsis_crassicaudata/T2T_ONT_genome/lambda/"
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/fastq/trimmed/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/"
mkdir -p "${OUTPUT_DIR}"
cd "${INPUT_DIR}"

for read_1 in Dunart-d9-P3*_L001_R1_001_val_1.fq.gz; do
    echo "$read_1"
    read_2="${read_1/_R1_001_val_1/_R2_001_val_2}"
    echo "$read_2"
    read_basename=$(basename "$read_1")
    read_basename_no_ext="${read_basename%.fq.gz}"
    output_bam="${OUTPUT_DIR}/${read_basename_no_ext}_bismark_bt2_pe.bam"
    if [ ! -f "$output_bam" ]; then
        CMD="${BISMARK}bismark --non_directional --unmapped --multicore 2 \
        -o \"${OUTPUT_DIR}\" --genome \"${GENOME_DIR}\" \
        -1 \"$read_1\" -2 \"$read_2\""
        echo "$CMD"
        eval "$CMD"
    else
        echo "Skipping $read_1 – output file already exists: $output_bam"
    fi
done

###--- map unmapped reads using single end alignment
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/merged_SE_fastq/"
cd ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

# combine read 1 and read 2
for sample in Dunart-d9-P3-*L001_R1_001_val_1.fq.gz_unmapped_reads_1.fq.gz; do
    base=${sample%%_R1_001_val_1.fq.gz_unmapped_reads_1.fq.gz}
    r1="${base}_R1_001_val_1.fq.gz_unmapped_reads_1.fq.gz"
    r2="${base}_R2_001_val_2.fq.gz_unmapped_reads_2.fq.gz"
    output="${base}_remainder.fq.gz"
    cat "$r1" "$r2" > ${OUTPUT_DIR}"$output"
done

# map
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/merged_SE_fastq/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/merged_SE_bam/"
cd ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

for read in Dunart-d9-P3*L001*_remainder.fq.gz; do
    echo "$read"
    read_basename=$(basename "$read")
    output_bam="${OUTPUT_DIR}/${read_basename%.fq.gz}_bismark_bt2.bam"
    if [ ! -f "$output_bam" ]; then
        CMD="${BISMARK}bismark --non_directional --multicore 2 \
        -o \"${OUTPUT_DIR}\" --genome \"${GENOME_DIR}\" \
        \"$read\""
        echo "$CMD"
        eval "$CMD"
    else
        echo "Skipping $read – output file already exists: $output_bam"
    fi
done

###--- deduplicate
# PAIRED END
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/dedup/"
mkdir -p "${OUTPUT_DIR}"
cd "${INPUT_DIR}"

for input_L001 in Dunart-d9-P3*_R1_001_val_1_bismark_bt2_pe.bam; do
    [ -e "$input_L001" ] || continue
    echo "$input_L001"
    input_L002="${input_L001/_L001_/_L002_}"
    if [ ! -f "$input_L002" ]; then
        echo "Warning: Paired file not found for $input_L001 (expected $input_L002). Skipping."
        continue
    fi
    # extract just the sample prefix
    # assumes format: NAME-EXTRA-SAMPLE-<rest>
    # splits on "-" and joins first 4 fields
    sample_prefix=$(basename "$input_L001" | cut -d'-' -f1-4)
    sample_id="${sample_prefix}_PE"
    output="${OUTPUT_DIR}/${sample_id}.multiple.deduplicated.bam"
    if [ ! -f "$output" ]; then
        echo "Running deduplication on $input_L001 and $input_L002"
        "${BISMARK}deduplicate_bismark" --bam --multiple \
            --output_dir "${OUTPUT_DIR}" -o "${sample_id}" \
            "$input_L001" "$input_L002"
    else
        echo "Skipping $input_L001 and $input_L002 – output file already exists: $output"
    fi
done

# SINGLE END
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/merged_SE_bam/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/dedup/"
mkdir -p "${OUTPUT_DIR}"
cd "${INPUT_DIR}"

for input_L001 in Dunart-d9-P3*_remainder_bismark_bt2.bam; do
    [ -e "$input_L001" ] || continue
    echo "$input_L001"
    input_L002="${input_L001/_L001_/_L002_}"
    if [ ! -f "$input_L002" ]; then
        echo "Warning: Paired file not found for $input_L001 (expected $input_L002). Skipping."
        continue
    fi
    # extract just the sample prefix
    # assumes format: NAME-EXTRA-SAMPLE-<rest>
    # splits on "-" and joins first 4 fields
    sample_prefix=$(basename "$input_L001" | cut -d'-' -f1-4)
    sample_id="${sample_prefix}_SE"
    output="${OUTPUT_DIR}/${sample_id}.multiple.deduplicated.bam"
    if [ ! -f "$output" ]; then
        echo "Running deduplication on $input_L001 and $input_L002"
        "${BISMARK}deduplicate_bismark" --bam --multiple \
            --output_dir "${OUTPUT_DIR}" -o "${sample_id}" \
            "$input_L001" "$input_L002"
    else
        echo "Skipping $input_L001 and $input_L002 – output file already exists: $output"
    fi
done

###--- bismark methylation extractor
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/dedup/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_extractor/"
mkdir -p "${OUTPUT_DIR}"
cd "${INPUT_DIR}"

for input in Dunart-d9-P3*.deduplicated.bam; do
    sample_name=$(basename "$input" .deduplicated.bam)
    output_cov="${OUTPUT_DIR}/${sample_name}.bismark.cov.gz"
    if [ ! -f "$output_cov" ]; then
        echo "Running methylation extraction on $input"
        # determine if BAM is paired-end or single-end
        if [[ "$input" == *"_PE."* ]]; then
            mode="--paired-end"
        else
            mode="--single-end"
        fi
        "${BISMARK}bismark_methylation_extractor" --gzip --bedGraph \
            "$mode" \
            --output_dir "${OUTPUT_DIR}" \
            "$input"
    else
        echo "Skipping $input – output file already exists: $output_cov"
    fi
done

###--- get CpG methylation calls in bedGraph format
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_extractor/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_bedGraph/"
mkdir -p "${OUTPUT_DIR}"
cd "${INPUT_DIR}"

# group all relevant input files
files=(CpG_CTOB_*d9-P3*.deduplicated.txt.gz CpG_CTOT_*d9-P3*.deduplicated.txt.gz CpG_OB_*d9-P3*.deduplicated.txt.gz CpG_OT_*d9-P3*.deduplicated.txt.gz)

# extract unique sample names and group matching files
declare -A sample_files

for file in "${files[@]}"; do
    # extract sample name (without prefix and PE/SE info)
    sample=$(echo "$file" | sed -E 's/CpG_(CTOB|CTOT|OB|OT)_//; s/_(PE|SE)\.multiple\.deduplicated\.txt\.gz//')
    sample_files["$sample"]+="$file "
done

# process each sample
for sample in "${!sample_files[@]}"; do
    input_list=${sample_files[$sample]}
    echo "Processing sample: $sample"
    echo "Input files: $input_list"

    "${BISMARK}bismark2bedGraph" \
        --dir "${OUTPUT_DIR}" \
        -o "${sample}.bedGraph" \
        $input_list
done
