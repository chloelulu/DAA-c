#!/bin/bash
cd /research/bsi/projects/staff_analysis/m216453/correlated/code/

for i in {1..100};
do
for method in MaAsLin2 GLMMPQL glmernb NBMM ZIGMM ZINBMM glmmadaptive glmmTMB LDM ZIBR LinDA
do
      qsub -N "${method}-${i}" -j y -cwd -q 1-day -m abe -l h_vmem=8G -V -o tmp/"${method}-${i}".out2 -e tmp/"${method}-${i}".err2 -b y R CMD BATCH \"--args ${method} ${i}\" subset_smoker.R tmp/"${method}-${i}".Rout2
      qsub -N "${method}-${i}" -j y -cwd -q 1-day -m abe -l h_vmem=8G -V -o tmp/"${method}-${i}".out.Nicholas -e tmp/"${method}-${i}".err.Nicholas -b y R CMD BATCH \"--args ${method} ${i}\" subset_Nicholas2013.R tmp/"${method}-${i}".Rout.Nicholas
      qsub -N "${method}-${i}" -j y -cwd -q 1-day -m abe -l h_vmem=8G -V -o tmp/"${method}-${i}".out.IBD -e tmp/"${method}-${i}".err.IBD -b y R CMD BATCH \"--args ${method} ${i}\" subset_IBD2017.R tmp/"${method}-${i}".Rout.IBD
done
done
