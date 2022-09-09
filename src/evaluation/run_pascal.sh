 #!/bin/bash

#pvaldir=${1?Error: No path to pval given}
#genescoredir=${2?Error: No path to genescore given}
#modulefile=${3?Error: No module file given}
moduledir=${1?Error: No module path given}
outpath=${2?Error: No output path given}
ppi=${3?Error: No PPI specified (enter 1 for string, 2 for inweb)}

if [[ $ppi -eq 1 ]]; then
	setting=dream11_settings/settings_leaderboard_v3--1_ppi_anonym_v2.txt
elif [[ $ppi -eq 2 ]]; then
	setting=dream11_settings/settings_leaderboard_v3--2_ppi_anonym_v2.txt
else
   echo "invalid ppi"
fi


#./Pascal --set=dream11_settings/settings_leaderboard_v3--2_ppi_anonym_v2.txt --runpathway=on --genescoring=sum\
# --pval=$pvaldir/EUR.GLGC.jointGwasMc_TC.txt.gz\
# --genescorefile=$genescoredir/EUR.GLGC.jointGwasMc_TC.sum.genescores.txt\
# --genesetfile $modulefile --outdir $outpath
cd ../../data/DREAM/5_scoring_tools/PASCAL

genescoredir=../../4_gwas_datasets/gene_scores
pvaldir=../../4_gwas_datasets/snp_pvalues

#the path is provided relative to where run_pascal is. path needed should be relative to PASCAL
moduledir="../../../../src/evaluation/${moduledir}"
outpath="../../../../src/evaluation/${outpath}"

for k in $moduledir/module_*.txt
do
	modulefile=$k
	for i in $pvaldir/*.txt.gz
	do
		len=${#i}
		echo $i
		j=${i::len-7}.sum.genescores.txt
		#echo $j
		j=${j:34}
		#echo $j
		j=$genescoredir/$j
		echo $j

		./Pascal --set=$setting --runpathway=on --genescoring=sum\
		 --pval=$i\
		 --genescorefile=$j\
		 --genesetfile $modulefile --outdir $outpath
	done
done