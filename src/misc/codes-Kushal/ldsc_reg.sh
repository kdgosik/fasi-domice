annot_cell=/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Miao_DO_propQTL
baseline_cell=/data/deyk/kushal/ENCODE_Flagship2023/ANNOTATIONS/Baselines
baseline_version=baseline_Epi_hg38
ldsc_path=/data/deyk/kushal/ldsc/
weights_path=/data/deyk/kushal/extras/1000G_EUR_Phase3_hg38/weights
freq_path=/data/deyk/kushal/extras/1000G_EUR_Phase3_hg38/plink_files
#sumstats_cell=/n/groups/price/ldsc/sumstats_formatted
sumstats_cell=/data/deyk/kushal/ENCODE_Flagship2023/sumstats
output_cell_pre=/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/Miao_DO_propQTL/
IFS="
"

#sumstats_taskfile=/n/groups/price/kushal/singlecellLDSC/data/traits_bio.txt
sumstats_taskfile=/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt  
annot_taskfile=/data/deyk/kushal/Miao_DO/data/miao2.txt

module load conda2
source activate ldsc

if [ ! -d $output_cell_pre ]
then
    mkdir $output_cell_pre
fi

output_cell=$output_cell_pre/$baseline_version

if [ ! -d $output_cell ]
then
    mkdir $output_cell
fi

echo $output_cell
for line in `cat $annot_taskfile | awk '{print $1}' | sort | uniq`;
do
    annot_module=`echo $line | awk '{print $1}'`
    echo $annot_cell $annot_module
    if [ ! -d $annot_cell/$annot_module ]
    then
        echo "Error: annotation module directory not found" > ldsc_logfile.log
        exit 100
    fi
    if [ ! -d $output_cell/$annot_module ]
    then
        mkdir $output_cell/$annot_module
    fi
    for ll in `ls $annot_cell/$annot_module | awk '{print $1}' | sort | uniq`;
    do
	annot_dir=`echo $ll | awk '{print $1}'`
	echo $annot_dir
	if [ ! -d $annot_cell/$annot_module/$annot_dir ]
	then
            echo "Error: annotation module directory not found" > ldsc_logfile.log
            exit 101
	fi
	if [ ! -d $output_cell/$annot_module/$annot_dir ]
	then
            mkdir $output_cell/$annot_module/$annot_dir
	fi
	for step in `cat $sumstats_taskfile | awk '{print $1}' | sort | uniq`;
	do
            sumstats_file=`echo $step | awk '{print $1}'`
            echo $sumstats_cell $sumstats_file
            if [ ! -f $sumstats_cell/$sumstats_file ]
            then
		echo "Error: sumstats file not found" > ldsc_logfile.log
		exit 102
            fi
            if [ ! -f $output_cell/$annot_module/$annot_dir/$sumstats_file.results ]
            then
            cmd="~/.conda/envs/ldsc/bin/python $ldsc_path/ldsc.py  --h2 $sumstats_cell/$sumstats_file --ref-ld-chr $annot_cell/$annot_module/$annot_dir/$annot_dir.,$baseline_cell/$baseline_version/baselineLD.  --frqfile-chr $freq_path/1000G.EUR.hg38. --w-ld-chr $weights_path/weights.hm3_noMHC. --overlap-annot --print-coefficients --print-delete-vals --out $output_cell/$annot_module/$annot_dir/$sumstats_file"
            bsub -W 300 -R "rusage[mem=20]" -e reg_max.err -o reg_max.out -n 1 "$cmd"
            fi
	done
     done
done



