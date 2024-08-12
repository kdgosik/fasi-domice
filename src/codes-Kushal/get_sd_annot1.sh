cellname=/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Miao_DO_propQTL
index=1
module load gcc/10.2.0
module load R

IFS="
"

cmd="Rscript get_sd_annot.R  $cellname $index"
bsub -W 450 -R "rusage[mem=20]" -e getsd.err -o getsd.out -n 2 "$cmd"




