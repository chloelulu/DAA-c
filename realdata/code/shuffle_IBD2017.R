args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  method = args[1]
  method = as.character(method)
  i = as.integer(args[2])}



setwd('/research/bsi/projects/staff_analysis/m216453/correlated/')
source('code/CompetingMethods1.R') ## can adjust Z
load('Data/IBD_2017_dat.Rdata')

load('Data/dat.IBD2017.shuffle1000.meta1.Rdata')
dat$meta.dat <- meta.list[[i]]
dat$counts <- dat$counts[,rownames(dat$meta.dat), drop =F]
dat$prop <- dat$prop[,rownames(dat$meta.dat), drop =F]
dat$gmpr.size.factor <- dat$gmpr.size.factor[rownames(dat$meta.dat)]


## Initially I tried comparing  MZ vs DZ twin pairs, age is confounding variable, random effect: SubjectID 
methods_funs <- list('glmernb'='glmernb.func','lme4'='lme4.func','GLMMPQL'='GLMMPQL.func','ZIBR'='ZIBR.func','LDM'='LDM.func','NBMM'='NBMM.func',
                     'MaAsLin2'='MaAsLin2.func','IFAA' = 'IFAA.func',
                     'ZIGMM'='ZIGMM.func','ZINBMM'='ZINBMM.func','metamicrobiomeR'='metamicrobiomeR.func','glmmadaptive'='glmmadaptive.func',
                     'glmmTMB'='glmmTMB.func','LinDA'='LinDA.func','lme' = 'lme.func','lmer'='lmer.func')

wrapper <- match.fun(methods_funs[[method]])
if(method=='MaAsLin2'){
  wrapper.obj <- wrapper(dat,cutoff=0.05,output = '/research/bsi/projects/staff_analysis/m216453/')
}else{
  wrapper.obj <- wrapper(dat,cutoff=0.05)
}

save(wrapper.obj, file = paste0('/research/bsi/projects/staff_analysis/m216453/correlated/Data/realshuffle/', method,'_',i, '_IBD2017.res.Rdata'))


