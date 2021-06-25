
cd /projectsp/foran/Team_Chan/STAR_FUSION_results/FFPE-new-data-10-2019/code-R-pick-up-OUTLIERS-2019-12-30

#
sbatch -p perceval -c 1 -n 1 -t 3-00:00:00 --mem=124000 --export=ALL --job-name="FGFR-gene" --error=log/"FGFR-gene".err --output=log/"FGFR-gene".out Rscript 1-all-pick-up-OUTLIERS-2019-12-30-gene.R

#
sbatch -p perceval -c 1 -n 1 -t 3-00:00:00 --mem=124000 --export=ALL --job-name=f-"FGFR-isoform" --error=log/"FGFR-isoform".err --output=log/"FGFR-isoform".out Rscript 2-all-pick-up-OUTLIERS-2019-12-30--isoforms.R


