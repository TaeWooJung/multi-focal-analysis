normal=$1
tumor=$2
output_dir=$3

reference=/home/bioit/mlarmuse/Documents/Data/Almac_project/CopenHagen_oligo/scripts_variant_calling/human_g1k_v37_decoy.fasta

start_time=$(date)


#if [ ! -f "${normal}.pileup" ]; then
#    samtools mpileup -f $reference -r ${normal} > ${normal}.pileup
#fi

#if [ ! -f "${normal}.pileup" ]; then
#    samtools mpileup -f $reference -r ${tumor} > ${tumor}.pileup
#fi

#samtools mpileup -f $reference ${normal} > ${normal}.pileup
#samtools mpileup -f $reference ${tumor} > ${tumor}.pileup

echo "Starting the variant calling: "

varscan somatic ${normal}.pileup ${tumor}.pileup ${output_dir} --min-var-freq 0

# rm normal_tumor_${sample_name}.pileup

echo 'Performing the filtering of the variants: '
echo ${output_file}
varscan somaticFilter ${output_dir}.snp --min-var-freq 0 --indel-file ${output_dir}.indel --output-file ${output_dir}.somaticFilter

echo 'This code started: '
echo ${start_time}
echo 'This code terminated: '
date
