annot_cell=/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Joint_Model_eQTL
baseline_cell=/data/deyk/kushal/ENCODE_Flagship2023/ANNOTATIONS/Baselines/baseline_Epi_hg38
#baseline_cell=/n/groups/price/kushal/LDSC/baselineLD_v2.1
ldsc_path=/data/deyk/kushal/ldsc/
#weights_path=/n/groups/price/kushal/LDSC/weights_hm3_no_hla
weights_path=/data/deyk/kushal/extras/1000G_EUR_Phase3_hg38/weights
freq_path=/data/deyk/kushal/extras/1000G_EUR_Phase3_hg38/plink_files
sumstats_cell=/data/deyk/kushal/ENCODE_Flagship2023/sumstats
output_cell=/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/Joint_Model_eQTL
IFS="
"

sumstats_taskfile=/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt
annot_taskfile=/data/deyk/kushal/Miao_DO/data/miao_combo.txt

module load conda2
source activate ldsc

for line in `cat $annot_taskfile | awk '{print $1}' | sort | uniq`;
do
    annotdir=`echo $line | awk '{print $1}'`
    echo $annot_cell $annotdir
    if [ ! -d $annot_cell/$annotdir ]
    then 
	echo "Error: annotation directory not found" > ldsc_logfile.log
	exit 100
    fi

    if [ ! -d $output_cell/$annotdir ]
    then
	mkdir $output_cell/$annotdir
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
#	if [ ! -f $output_cell/$annotdir/$sumstats_file.results ]
#	then
	    cmd="python $ldsc_path/ldsc.py  --h2 $sumstats_cell/$sumstats_file --ref-ld-chr $annot_cell/$annotdir/$annotdir.,$baseline_cell/baselineLD.  --frqfile-chr $freq_path/1000G.EUR.hg38. --w-ld-chr $weights_path/weights.hm3_noMHC. --overlap-annot --print-coefficients --print-delete-vals --out $output_cell/$annotdir/$sumstats_file"
	    #sbatch --time=120:00 --mem=10000 --output=reg_max.out --error=reg_max.err -p short -c 1 --wrap="$cmd"
            bsub -W 300 -R "rusage[mem=20]" -e reg_max.err -o reg_max.out -n 1 "$cmd"
#	fi
    done
done






#python ldsc.py  --h2 /n/groups/price/ldsc/sumstats_formatted/PASS_BMI1.sumstats --ref-ld-chr /n/groups/price/kushal/1000G/ANNOTATIONS_LD/H3K4me1_ABS_SCORE/H3K4me1_ABS_SCORE. --frqfile-chr /n/groups/price/kushal/LDSC/1000G_Phase3_frq/1000G.EUR.QC. --w-ld-chr /n/groups/price/kushal/LDSC/weights_hm3_no_hla/weights. --overlap-annot --print-coefficients --print-delete-vals --out /n/groups/price/kushal/OUTPUT
