annot_cell=/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Orr
ldsc_path=/data/deyk/kushal/ldsc/
bfile_path=/data/deyk/kushal/LDSC/1000G_EUR_Phase3_plink_hg38
hapmap_path=/data/deyk/kushal/LDSC/hapmap3_snps

IFS="
"

TASKFILE=/data/deyk/kushal/Miao_DO/data/orr.txt

module load conda2
source activate ldsc

for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
   annot_module=`echo $line | awk '{print $1}'`
   echo $annot_cell $annot_module
   for ll in `ls $annot_cell/$annot_module | awk '{print $1}' | sort | uniq`;
   do
       annot_dir=`echo $ll | awk '{print $1}'`
       echo $annot_dir
       if [ ! -d $annot_cell/$annot_module/$annot_dir ]
       then
	   mkdir $annot_cell/$annot_module/$annot_dir
       fi
       for chrom in {1..22}
       do
       if [ ! -f $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom.l2.ldscore.gz ]
       then
           cmd="~/.conda/envs/ldsc/bin/python $ldsc_path/ldsc.py --bfile $bfile_path/1000G.EUR.QC.$chrom --l2 --ld-wind-cm 1 --yes-really --annot $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom.annot.gz --print-snps $hapmap_path/hm.$chrom.snp --out $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom"
           bsub -W 300 -R "rusage[mem=20]" -e mega.err -o mega.out -n 1 "$cmd"
       fi
    done
  done
done

