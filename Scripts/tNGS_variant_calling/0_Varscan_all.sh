
sample_ids=(ID4 ID5 ID6 ID9 ID10 ID12)
path=/home/bioit/mlarmuse/Documents/Data/Almac_project/CopenHagen_oligo/tNGS_bams
write_dir=/home/bioit/mlarmuse/Documents/Data/Almac_project/CopenHagen_oligo/Strelka_output

for s in ${sample_ids[@]}
do
echo #####################################################
echo $s
sample_path=${path}/$s
echo $sample_path

normal_sample=$(ls ${sample_path}/*${s}-N*.bam)
echo $normal_sample

samtools index $normal_sample

tumor_samples=$(ls ${sample_path}/*${s}-T*.bam)

for ts in ${tumor_samples[@]}
do

echo $ts
sample_name=$(echo $ts | cut -f 5 -d\-)

echo $sample_name

outdir=${write_dir}/${sample_name}

bash callVarscanFastPaired.sh $normal_sample $ts $outdir

done

done
