MAGMA_CELL=/data/deyk/kushal/Placenta/MAGMA
RAW_CELL=/data/deyk/kushal/Miao_DO/data/MAGMA/GENE_CELL_0
SET_CELL=/data/deyk/kushal/Miao_DO/data/MAGMA
GENESET_CELL=/data/deyk/kushal/Miao_DO/data/MAGMA/SETS_MIAO_NONSPECIFIC
sumstats_taskfile=/data/deyk/kushal/Miao_DO/data/MAGMA/traits_miao.txt

for step in `cat $sumstats_taskfile | awk '{print $1}' | sort | uniq`;
do
  sumstats_file=`echo $step | awk '{print $1}'`
  echo $sumstats_file
  traitname=$sumstats_file
  if [ ! -f $GENESET_CELL/${traitname}.gsa.out ];
  then
    $MAGMA_CELL/magma --gene-results $RAW_CELL/$traitname.gene_trait.genes.raw --set-annot $SET_CELL/miao_Nonspecific_celltypes.set --out $GENESET_CELL/$traitname
  fi
done


