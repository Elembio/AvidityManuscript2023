set -euo pipefail
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#bash ~/git/wholegenomesequencing_bash/bash/run_fastq_highAT.sh -f s3://elembio-usw2-s3-runs/blajoie/data/workflows/fastq/20210322_P1-04_SI330-test/ -o s3://elembio-usw2-s3-runs/blajoie/projects/highAT/20210322_P1-04_SI330-test/ -t 0.9 -w ~/scratch/

maybe_download_prefix(){
    if [[ $1 == s3://* ]]
    then
        dirname=$(dirname ${1})
        aws s3 ls ${1} | awk '{print $4}' | xargs -I {} aws s3 sync --force-glacier-transfer --only-show-errors --exclude "*" --include {} ${dirname}/ ${2} 2> err
        if [ -s "err" ]
        then
            echo -e "\n--- ERROR downoading $1 ---\n" >&2
            cat err >&2
            exit 1
        fi
        echo $2
    else
        echo $1
    fi
}

maybe_download_dir(){
    if [[ $1 == s3://* ]]
    then
        aws s3 sync --force-glacier-transfer --only-show-errors ${1} ${2} 2> err
        if [ -s "err" ]
        then
            echo -e "\n--- ERROR downoading $1 ---\n" >&2
            cat err >&2
            exit 1
        fi
        echo $2
    else
        echo $1
    fi
}

locate_single_file(){    
    files=$(ls ${1}/*.${2} 2> /dev/null | wc -l)
    if (( $files > 1 ))
    then
        echo -e "\nERROR: found multiple files with exten-${2}" >&2
        ls ${1}/*.${2} >&2
        exit 1
    fi
    file=$(ls ${1}/*.${2} 2> /dev/null)
    echo $file
}

get_seeded_random(){
    seed="$1"
    openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
}

HOMOPOLYMER=false
while getopts r:b:g:i:o:hw: opt; do
    case $opt in
        r) BAM_URL="$OPTARG" ;;
        b) BED_URL="$OPTARG" ;;
        g) GENOME_URL="$OPTARG" ;;
        i) ID="$OPTARG" ;;
        o) OUTPUT_URL="$OPTARG" ;;
        h) HOMOPOLYMER=true;;
        w) WORK_ROOT="$OPTARG" ;;
        *) 
            echo "Invalid option -$OPTARG" >&2
            exit 1
    esac
done

: ${BAM_URL:?Missing -r}
: ${BED_URL:?Missing -b}
: ${GENOME_URL:?Missing -g}
: ${OUTPUT_URL:?Missing -o}
: ${ID:?Missing -i}
THREADS=${THREADS:-"4"}
CONCURRENCY=${CONCURRENCY:-"6"}
RATE=${RATE:-"0.1"}
SLOP=${SLOP:-"0"}
MIN_INTERVAL=${MIN_INTERVAL:-"1"}
MIN_PRE=${MIN_PRE:-"5"}
MIN_POST=${MIN_POST:-"5"}
PRE_CONTEXT=${PRE_CONTEXT:-"0"}
POST_CONTEXT=${POST_CONTEXT:-"0"}
WORK_ROOT=${WORK_ROOT:-"/"}
# add trailing slash if missing
[[ "${WORK_ROOT}" != */ ]] && WORK_ROOT="${WORK_ROOT}/"

now=$(date)
echo -e "\nStarting stack-reads-by-interval - ${now}"

echo -e "\n--- Inputs ---"
echo -e "SCRIPT_DIR": ${SCRIPT_DIR}""
echo "BAM_URL: ${BAM_URL}"
echo "BED_URL: ${BED_URL}"
echo "GENOME_URL: ${GENOME_URL}"
echo "OUTPUT_URL: ${OUTPUT_URL}"
echo "ID: ${ID}"
echo "THREADS: ${THREADS}"
echo "RATE: ${RATE}"
echo "SLOP: ${SLOP}"
echo "MIN_INTERVAL: ${MIN_INTERVAL}"
echo "MIN_PRE: ${MIN_PRE}"
echo "MIN_POST: ${MIN_POST}"
echo "PRE_CONTEXT: ${PRE_CONTEXT}"
echo "POST_CONTEXT: ${POST_CONTEXT}"
echo "HOMOPOLYMER: ${HOMOPOLYMER}"
echo "WORK_ROOT: ${WORK_ROOT}"

# get bam
echo -e "\nDownloading BAM: ${WORK_ROOT}input/bam"
bam_dir=$(maybe_download_prefix ${BAM_URL} ${WORK_ROOT}"input/bam")
echo -e "Downloading BED: ${WORK_ROOT}input/bed"
bed_dir=$(maybe_download_prefix ${BED_URL} ${WORK_ROOT}"input/bed")
echo -e "Downloading GENOME: ${WORK_ROOT}input/genome"
genome_dir=$(maybe_download_prefix ${GENOME_URL} ${WORK_ROOT}"input/genome")

# check dir contents
echo -e "\nDirectory contents - ${bam_dir}"
ls -ld ${bam_dir}/*
echo -e "\nDirectory contents - ${bed_dir}"
ls -ld ${bed_dir}/*
echo -e "\nDirectory contents - ${genome_dir}"
ls -ld ${genome_dir}/*

# create output dir
output_dir=${WORK_ROOT}"output" 
rm -rf ${output_dir}
mkdir -p ${output_dir}

logfile=${output_dir}/run.log
exec > >(tee $logfile)
exec 2>&1

# run stack-reads-by-interval
bam_file=$(locate_single_file ${bam_dir} "bam")
bed_file=$(locate_single_file ${bed_dir} "bed.gz")
genome_file=$(locate_single_file ${genome_dir} "genome")

echo -e "\n--- Starting stack-reads-by-interval ---"
echo -e "bam: ${bam_file}"
echo -e "bed: ${bed_file}"
echo -e "genome: ${genome_file}"

HOMOPOLYMER_OPTION=""
if [ "$HOMOPOLYMER" == "true" ]
then
    HOMOPOLYMER_OPTION="--homopolymer"
fi

genome_size=`cat ${genome_file} | awk '{ sum += $2 } END { print sum }'`
bed_size=`zcat ${bed_file} | awk '{ sum += $3-$2 } END { print sum }'`
num_intervals=`zcat ${bed_file} | wc -l`
echo -e "genome_size: ${genome_size}"
echo -e "bed_size: ${bed_size}"
echo -e "num_intervals: ${num_intervals}"

max_intervals=1000000
if (( $(echo "$num_intervals > ${max_intervals}" | bc -l) ))
then
    echo -e "\nBED [${num_intervals} > ${max_intervals}] is too large ... randomly downsampling\n"
    downsampled_bed_file=`basename ${bed_file} | sed s/.bed.gz//`
    downsampled_bed_file="${downsampled_bed_file}_downsampled.bed.gz"
    echo -e "downsampled_bed: ${downsampled_bed_file}"
    zcat ${bed_file} | grep -v "#" | shuf --random-source=<(get_seeded_random 42) -n ${max_intervals} | bedtools sort -g ${genome_file} | bgzip > ${downsampled_bed_file}
    bed_file=${downsampled_bed_file}
    num_intervals=`zcat ${bed_file} | wc -l`
    echo -e "num_intervals (modified): ${num_intervals}"
fi

overlapped_bam_file=`basename ${bam_file} | sed s/.bam//`
overlapped_bam_file="${overlapped_bam_file}_overlap.bam"

# intersect bam with bed, produce overlapped (smaller) bam file
echo -e "\noverlapping bam with bed [samtools view -hb -L ${bed_file} -@${THREADS} ${bam_file} -o ${overlapped_bam_file} && samtools index -@${THREADS} ${overlapped_bam_file}]"
samtools view -hb -L ${bed_file} -@${THREADS} ${bam_file} -o ${overlapped_bam_file} && samtools index -@${THREADS} ${overlapped_bam_file}
echo -e "overlapped_bam_file: ${overlapped_bam_file}"

bam_size=`samtools idxstats ${overlapped_bam_file} | awk '{ sum += $3 } END { print sum }'`
echo -e "bam_size: ${bam_size}"

if (( $(echo "$bam_size > 30000000" | bc -l) ))
then
    RATE=`echo "scale=10; 30000000/${bam_size}" | bc`
    echo "RATE-OVERRIDE: ${RATE}"
fi

# remove any existing job files
rm -f READ_STACK.jobs

# create read_stack.py commands per chr
while read p; 
    do chr=`echo $p | awk '{print $1}'`; 
    #python3 python/stack_reads_by_interval.py
    echo "python3 ${SCRIPT_DIR}/../python/stack_reads_by_interval.py --bam ${overlapped_bam_file} --bed ${bed_file} --chr ${chr} --rate ${RATE} --slop ${SLOP} --min-interval ${MIN_INTERVAL} --min-pre ${MIN_PRE} --min-post ${MIN_POST} --pre-context ${PRE_CONTEXT} --post-context ${POST_CONTEXT} --output-prefix ${output_dir}/${ID}__${chr} ${HOMOPOLYMER_OPTION}" >> READ_STACK.jobs
done < ${genome_file}

# print all read_stack jobs
cat READ_STACK.jobs

# run all read_stack jobs
cat READ_STACK.jobs | parallel -u -j ${CONCURRENCY}

# combine all chr.offset-error.tsv
cat ${output_dir}/${ID}__chr*.offset-error.tsv | sort -k1,1 -k2,2n > ${output_dir}/${ID}.offset-error.tsv
rm ${output_dir}/${ID}__chr*.offset-error.tsv

# combine all chr.interval-error.tsv
cat ${output_dir}/*__chr*.interval-error.tsv | sort -k1,1n > ${output_dir}/${ID}.interval-error.tsv
rm ${output_dir}/*__chr*.interval-error.tsv

echo -e "\n--- Completed stack-reads-by-interval ---"

now=$(date)
echo -e "\nCompleted stack-reads-by-interval - ${now}"



