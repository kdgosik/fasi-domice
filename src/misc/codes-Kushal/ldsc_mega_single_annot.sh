annot_cell=/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Joint_Model_eQTL
ldsc_path=/data/deyk/kushal/ldsc/
bfile_path=/data/deyk/kushal/LDSC/1000G_EUR_Phase3_plink_hg38
hapmap_path=/data/deyk/kushal/LDSC/hapmap3_snps

IFS="
"

TASKFILE=/data/deyk/kushal/Miao_DO/data/miao_combo.txt

module load conda2                                                                                    
source activate ldsc

for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
   annotdir=`echo $line | awk '{print $1}'`
   echo $annot_cell $annotdir
   if [ ! -d $annot_cell/$annotdir ]
   then
       mkdir $annot_cell/$annotdir
   fi
   for chrom in {1..22}
   do
       if [ ! -f $annot_cell/$annotdir/$annotdir.$chrom.l2.ldscore.gz ]
       then
           cmd="~/.conda/envs/ldsc/bin/python $ldsc_path/ldsc.py --bfile $bfile_path/1000G.EUR.QC.$chrom --l2 --ld-wind-cm 1 --yes-really --annot $annot_cell/$annotdir/$annotdir.$chrom.annot.gz --print-snps $hapmap_path/hm.$chrom.snp --out $annot_cell/$annotdir/$annotdir.$chrom"
           #sbatch --time=90:00 --mem=10000 --output=mega.out --error=mega.err -p short -c 1 --wrap="$cmd"
           bsub -W 300 -R "rusage[mem=20]" -e mega.err -o mega.out -n 1 "$cmd"
       fi
    done
done
