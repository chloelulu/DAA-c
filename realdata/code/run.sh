#!/bin/bash
cd /research/bsi/projects/staff_analysis/m216453/correlated/code/

for method in LinDA MaAsLin2 GLMMPQL glmernb NBMM ZIGMM ZINBMM glmmadaptive glmmTMB LDM ZIBR
do
   qsub -N "${method}" -j y -cwd -q 1-hour -m abe -l h_vmem=8G -V -o tmp/"${method}".out1 -e tmp/"${method}".err1 -M yang.lu@mayo.edu -b y R CMD BATCH \"--args ${method}\" smoker.R tmp/"${method}".Rout1
   qsub -N "${method}" -j y -cwd -q 1-hour -m abe -l h_vmem=8G -V -o tmp/"${method}".out.Nicholas -e tmp/"${method}".err.Nicholas -M yang.lu@mayo.edu -b y R CMD BATCH \"--args ${method}\" Nicholas_2013.R tmp/"${method}".Rout.Nicholas
   qsub -N "${method}" -j y -cwd -q 1-day -m abe -l h_vmem=8G -V -o tmp/"${method}".out.IBD -e tmp/"${method}".err.IBD  -M yang.lu@mayo.edu -b y R CMD BATCH \"--args ${method}\" IBD2017.R tmp/"${method}".Rout.IBD
done
