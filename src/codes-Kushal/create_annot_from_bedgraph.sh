module load gcc/10.2.0
module load conda
source activate ldsc
bedfile_path=/data/deyk/kushal/Miao_DO/data/BEDFILES/Miao_DO_propQTL
bimfile_path=/data/deyk/kushal/extras/BIMS_hg38
annot_path=/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Miao_DO_propQTL
#ldsc_path=/n/groups/price/kushal/LDSC/ldsc

IFS="
"

#TASKFILE=/n/groups/price/kushal/LDSC-Average/data/Enames.txt
TASKFILE=/data/deyk/kushal/Miao_DO/data/miao2.txt
for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
    name=`echo $line | awk '{print $1}'`
    if [ ! -d $annot_path/$name ]
    then
	mkdir $annot_path/$name
    fi
    for bedline in `ls $bedfile_path/$name/ | cat | sort | uniq | awk -F .'bed' '{print $1}'`;
    do
	bedname=`echo $bedline | awk '{print $1}'`
	if [ ! -d $annot_path/$name/$bedname ]
	then
	    mkdir $annot_path/$name/$bedname
	fi
	if [ ! -f $annot_path/$name/$bedname/$bedname.22.annot.gz ]
	then
	    cmd="~/.conda/envs/ldsc/bin/python  make_annot_combine_from_bedgraph.py --bedname $bedname --bedfile_path $bedfile_path/$name --bimfile_path $bimfile_path --annot_path $annot_path/$name/$bedname"
            bsub -W 300 -R "rusage[mem=20]" -e annot.err -o annot.out -n 1 "$cmd"
	fi
    done
done
