args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  method = args[1]
  method = as.character(method)
  i = as.integer(args[2])}


setwd('/research/bsi/projects/staff_analysis/m216453/correlated/')
source('code/CompetingMethods1.R') ## can adjust Z
load('Data/dat.smoker_2010.Rdata')
load('Data/dat.smoker_2010.shuffle1000.meta.Rdata')
dat$meta.dat <- meta.list[[i]]
dat$meta.dat$X <- as.numeric(dat$meta.dat$X)-1
dat$meta.dat$time <- as.numeric(dat$meta.dat$time)-1



## Initially I tried comparing  MZ vs DZ twin pairs, age is confounding variable, random effect: SubjectID 
methods_funs <- list('glmernb'='glmernb.func','lme4'='lme4.func','GLMMPQL'='GLMMPQL.func','ZIBR'='ZIBR.func','LDM'='LDM.func','NBMM'='NBMM.func',
                     'MaAsLin2'='MaAsLin2.func','ZIGMM'='ZIGMM.func','ZINBMM'='ZINBMM.func','glmmadaptive'='glmmadaptive.func','IFAA' = 'IFAA.func',
                     'glmmTMB'='glmmTMB.func','LinDA'='LinDA.func','lme' = 'lme.func','lmer'='lmer.func')
# methods <- c('GLMMPQL','glmernb','ZIBR','LDM','LinDA','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMB','MaAsLin2')

wrapper <- match.fun(methods_funs[[method]])
if(method=='MaAsLin2'){
  wrapper.obj <- wrapper(dat,cutoff=0.05,output = '/research/bsi/projects/staff_analysis/m216453/', data.type ='repeated.measure')
}else{
  wrapper.obj <- wrapper(dat,cutoff=0.05, data.type ='repeated.measure')
}
save(wrapper.obj, file = paste0('/research/bsi/projects/staff_analysis/m216453/correlated/Data/realshuffle/', method,'_',i, '_smoker.res.Rdata'))
