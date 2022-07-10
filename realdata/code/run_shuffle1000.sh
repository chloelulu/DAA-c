#!/bin/bash
cd /research/bsi/projects/staff_analysis/m216453/correlated/code/

for i in {1..1000};
do
for method in LinDA MaAsLin2 GLMMPQL glmernb NBMM ZIGMM ZINBMM glmmadaptive glmmTMB LDM ZIBR
do
      qsub -N "${method}-${i}" -j y -cwd -q 1-day -m abe -l h_vmem=8G -V -o "${method}-${i}".out.smoker -e "${method}-${i}".err.smoker -b y R CMD BATCH \"--args ${method} ${i}\" shuffle_smoker.R "${method}-${i}".Rout.smoker
      qsub -N "${method}-${i}" -j y -cwd -q 1-day -m abe -l h_vmem=8G -V -o "${method}-${i}".out.Nicholas -e "${method}-${i}".err.Nicholas -b y R CMD BATCH \"--args ${method} ${i}\" shuffle_Nicholas2013.R "${method}-${i}".Rout.Nicholas
      qsub -N "${method}-${i}" -j y -cwd -q 1-day -m abe -l h_vmem=8G -V -o "${method}-${i}".out.IBD -e "${method}-${i}".err.IBD -b y R CMD BATCH \"--args ${method} ${i}\" shuffle_IBD2017.R "${method}-${i}".Rout.IBD

done
done
