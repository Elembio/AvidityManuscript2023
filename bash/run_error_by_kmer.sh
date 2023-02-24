set -euo pipefail
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

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

while getopts r:b:g:k:i:o:w: opt; do
    case $opt in
        r) BAM_URL="$OPTARG" ;;
        b) BED_URL="$OPTARG" ;;
        g) GENOME_URL="$OPTARG" ;;
        k) KMER_LENGTH="$OPTARG" ;;
        i) ID="$OPTARG" ;;
        o) OUTPUT_URL="$OPTARG" ;;
        w) WORK_ROOT="$OPTARG" ;;
        *) 
            echo "Invalid option -$OPTARG" >&2
            exit 1
    esac
done

: ${BAM_URL:?Missing -r}
: ${GENOME_URL:?Missing -g}
: ${KMER_LENGTH:?Missing -k}
: ${OUTPUT_URL:?Missing -o}
: ${ID:?Missing -i}
BED_URL=${BED_URL:-""}
SAMPLE=${SAMPLE:-""}
THREADS=${THREADS:-"4"}
CONCURRENCY=${CONCURRENCY:-"6"}
RATE=${RATE:-"0.1"}
WORK_ROOT=${WORK_ROOT:-"/"}
# add trailing slash if missing
[[ "${WORK_ROOT}" != */ ]] && WORK_ROOT="${WORK_ROOT}/"

now=$(date)
echo -e "\nStarting compute-error-by-kmer - ${now}"

echo -e "\n--- Inputs ---"
echo "SCRIPT_DIR: ${SCRIPT_DIR}"
echo "BAM_URL: ${BAM_URL}"
echo "BED_URL: ${BED_URL}"
echo "KMER_LENGTH: ${KMER_LENGTH}"
echo "GENOME_URL: ${GENOME_URL}"
echo "SAMPLE: ${SAMPLE}"
echo "OUTPUT_URL: ${OUTPUT_URL}"
echo "ID: ${ID}"
echo "THREADS: ${THREADS}"
echo "RATE: ${RATE}"
echo "WORK_ROOT: ${WORK_ROOT}"

# get bam
echo -e "\nDownloading BAM: ${WORK_ROOT}input/bam"
bam_dir=$(maybe_download_prefix ${BAM_URL} ${WORK_ROOT}"input/bam")
echo -e "Downloading BED: ${WORK_ROOT}input/bed"
bed_dir=$(maybe_download_prefix ${BED_URL} ${WORK_ROOT}"input/bed")
echo -e "Downloading GENOME: ${WORK_ROOT}input/genome"
genome_dir=$(maybe_download_prefix ${GENOME_URL} ${WORK_ROOT}"input/genome")

# create output dir
output_dir=${WORK_ROOT}"output" 
rm -rf ${output_dir}
mkdir -p ${output_dir}

logfile=${output_dir}/run.log
exec > >(tee $logfile)
exec 2>&1

# run compute-error-by-kmer
bam_file=$(locate_single_file ${bam_dir} "bam")
bed_file=$(locate_single_file ${bed_dir} "bed.gz")
genome_file=$(locate_single_file ${genome_dir} "genome")

echo -e "\n--- Starting compute-error-by-kmer ---"
echo -e "kmerlength: ${KMER_LENGTH}"
echo -e "bam: ${bam_file}"
echo -e "bed: ${bed_file}"
echo -e "genome: ${genome_file}"

genome_size=`cat ${genome_file} | awk '{ sum += $2 } END { print sum }'`
echo -e "genome_size: ${genome_size}"

if [[ ! -z $bed_file ]]    
then
    bed_size=`zcat ${bed_file} | awk '{ sum += $3-$2 } END { print sum }'`
    num_intervals=`zcat ${bed_file} | wc -l`
    echo -e "bed_size: ${bed_size}"
    echo -e "num_intervals: ${num_intervals}"

    overlapped_bam_file=`basename ${bam_file} | sed s/.bam//`
    overlapped_bam_file="${overlapped_bam_file}_overlap.bam"

    # intersect bam with bed, produce overlapped (smaller) bam file
    echo -e "\noverlapping bam with bed [samtools view -hb -L ${bed_file} -@${THREADS} ${bam_file} -o ${overlapped_bam_file} && samtools index -@${THREADS} ${overlapped_bam_file}]"
    samtools view -hb -L ${bed_file} -@${THREADS} ${bam_file} -o ${overlapped_bam_file} && samtools index -@${THREADS} ${overlapped_bam_file}
    echo -e "overlapped_bam_file: ${overlapped_bam_file}"
else
    overlapped_bam_file=${bam_file}
fi

bam_size=`samtools idxstats ${overlapped_bam_file} | awk '{ sum += $3 } END { print sum }'`
echo -e "bam_size: ${bam_size}"

# remove any existing job files
rm -f KMER_ERROR.jobs

# create read_stack.py commands per chr
while read p; 
    do chr=`echo $p | awk '{print $1}'`; 
    #python3 python/stack_reads_by_interval.py
    echo "python3 ${SCRIPT_DIR}/../python/compute_error_by_kmer.py --bam ${overlapped_bam_file} --chr ${chr} --kmer_length ${KMER_LENGTH} --rate ${RATE} --output-prefix ${output_dir}/${ID}__${chr} --sample_name ${SAMPLE}"  >> KMER_ERROR.jobs
done < ${genome_file}

# print all read_stack jobs
cat KMER_ERROR.jobs

# run all read_stack jobs
cat KMER_ERROR.jobs | parallel -u -j ${CONCURRENCY}

# combine all chr.offset-error.tsv
echo "cat ${output_dir}/${ID}__chr*.kmer-error.tsv | sort -k1,1 -k2,2n | python3 ${SCRIPT_DIR}/../python/combine_error_by_kmer.py --kmer_length ${KMER_LENGTH} --output-prefix ${output_dir}/${ID}"

cat ${output_dir}/${ID}__chr*.kmer-error.tsv | sort -k1,1 -k2,2n | python3 ${SCRIPT_DIR}/../python/combine_error_by_kmer.py --kmer_length ${KMER_LENGTH} --output-prefix ${output_dir}/${ID}
rm ${output_dir}/${ID}__chr*.kmer-error.tsv

echo -e "\n--- Completed compute-error-by-kmer ---"

now=$(date)
echo -e "\nCompleted compute-error-by-kmer - ${now}"
