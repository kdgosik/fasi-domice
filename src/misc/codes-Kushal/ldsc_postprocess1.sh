annot_cell=/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Miao_DO_propQTL
output_cell_pre=/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/Miao_DO_propQTL/
baseline_version=baseline_Epi_hg38
output_cell=$output_cell_pre/$baseline_version

sumstats_taskfile=/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt
#sumstats_taskfile=/n/groups/price/kushal/singlecellLDSC/data/traits_bio.txt

IFS="
"

module load gcc/10.2.0
module load R

flag=0
index_in_results=1 ## which annotation to choose from the .results file in case of multiple annotations


for step in `cat $sumstats_taskfile | awk '{print $1}' | sort | uniq`;
do
sumstats_file=`echo $step | awk '{print $1}'`
echo $sumstats_file
sumstats_file2=${sumstats_file%.sumstats}

counter1=0
for step2 in `ls $output_cell | awk '{print $1}' | sort | uniq`;
do
    annot_name=`echo $step2 | awk '{print $1}'`
    if [ ! -f $output_cell/$annot_name/${sumstats_file2}_ldsc_postprocess.txt ]
    then
	counter1=$(($counter1+1))
    fi
done

if (( $counter1 > 0 ))
then
    echo $sumstats_file2
    cmd="Rscript ldsc_postprocess.R  $annot_cell $output_cell $sumstats_file $flag $index_in_results"
    bsub -W 270 -R "rusage[mem=20]" -e ldsc_post.err -o ldc_post.out -n 1 "$cmd"
    #sbatch --time=40:00 --mem=20000 --output=ldsc_post.out --error=ldsc_post.err -p short -c 1 --wrap="$cmd"
fi
done
