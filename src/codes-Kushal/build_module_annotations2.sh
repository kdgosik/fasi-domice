genescore_cell=/data/deyk/kushal/Miao_DO/data/Gene_Scores/
bed_cell=/data/deyk/kushal/Miao_DO/data/BEDFILES/
enhancer_tissue=BLD

module load gcc/10.2.0
module load R

IFS="
"

TASKFILE=/data/deyk/kushal/Miao_DO/data/Tnames.txt

for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
   temp=`echo $line | awk '{print $1}'`
   bed_dir=$bed_cell/$temp
   if [ ! -d $bed_dir ]
   then
       mkdir $bed_dir
   fi
   genescore_dir=$genescore_cell/$temp
   if [ ! -d $genescore_dir ]
   then
       mkdir $genescore_dir
   fi
   for ll in `ls -1 $genescore_dir | sed 's/\.txt//g' | awk '{print $1}' | sort | uniq`;
   do
      annot_name=`echo $ll | awk '{print $1}'`
      echo $temp $annot_name
#     if [ ! -f $bed_dir/$annot_name/100kb.bed ]
#      then
	  cmd="Rscript build_module_annotations.R  $genescore_dir $bed_dir $annot_name $enhancer_tissue"
          bsub -W 90 -R "rusage[mem=20]" -e geneS2G.err -o geneS2G.out -n 1 "$cmd"
#      fi
   done
done

