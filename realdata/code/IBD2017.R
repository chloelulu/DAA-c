args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  method = args[1]
  method = as.character(method)}


setwd('/research/bsi/projects/staff_analysis/m216453/correlated/')
source('code/CompetingMethods1.R') ## can adjust Z
load('Data/IBD_2017_dat.Rdata')



## Shuffle disease status
# meta <- dat$meta.dat
# meta.list <- list()
# for(i in 1:1000){
#   cat(i,'\n')
#   length(unique(meta1[meta1$X==1,'SubjectID']));length(unique(meta1[meta1$X==0,'SubjectID']))
#   ## Shuffle time
#   meta1 <- as.data.frame(meta %>% rownames_to_column('sam.ids') %>% group_by(X,SubjectID) %>% dplyr::mutate(time=sample(time, replace = F)))
#   ## sample the disease status of each subject
#   uni <- unique(meta1[,c(2:3)])
#   uni$X <- sample(uni$X)
#   meta2 <- uni %>% right_join(meta1%>% dplyr::select(-X)) %>% column_to_rownames('sam.ids')
#   meta.list[[i]] <- meta2
# }
# save(meta.list, file = 'Data/dat.IBD2017.shuffle1000.meta3.Rdata')



methods_funs <- list('glmernb'='glmernb.func','lme4'='lme4.func','GLMMPQL'='GLMMPQL.func','ZIBR'='ZIBR.func','LDM'='LDM.func',
                      'MaAsLin2'='MaAsLin2.func','NBMM'='NBMM.func', 'ZIGMM'='ZIGMM.func','ZINBMM'='ZINBMM.func',
                     'glmmadaptive'='glmmadaptive.func', 'glmmTMB'='glmmTMB.func','IFAA'='IFAA.func',
                    'LinDA'='LinDA.func','lme' = 'lme.func','lmer'='lmer.func')
# methods <- c('GLMMPQL','glmernb','ZIBR','LDM','LinDA','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMB','MaAsLin2')
wrapper <- match.fun(methods_funs[[method]])
if(method=='MaAsLin2'){
  wrapper.obj <- wrapper(dat,cutoff=0.05,output = '/research/bsi/projects/staff_analysis/m216453/')
}else{
  wrapper.obj <- wrapper(dat,cutoff=0.05)
}
save(wrapper.obj, file = paste0('/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/', method, '_IBD2017.res.Rdata'))


