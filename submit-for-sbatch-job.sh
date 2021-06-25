#!/bin/bash
#SBATCH -n 4
#SBATCH -c 8
#SBATCH -t 3-00:00:00
#SBATCH -p main
#SBATCH --mem=124000
#SBATCH --constraint=oarc 
#SBATCH --export=ALL
#SBATCH -D /projectsp/foran/Team_Chan/STAR_FUSION_results/FFPE-Sample-in-ORIEN-Analysis-Workflow-August-26-2020/analysis---leafcutter-counts---November-6-2020/code/only-extract-gene--python-get-kinase-genes/log
#SBATCH -o annotate-gene.out
#SBATCH -e annotate-gene.err
# this may work, check

# 
cd /projectsp/foran/Team_Chan/STAR_FUSION_results/FFPE-Sample-in-ORIEN-Analysis-Workflow-August-26-2020/analysis---leafcutter-counts---November-6-2020/code/only-extract-gene--python-get-kinase-genes
 
python 80-samples-only-extract-gene--python-annotate-gene-name.py