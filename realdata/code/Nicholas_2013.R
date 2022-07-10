args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  method = args[1]
  method = as.character(method)}


load('Data/dat.Nicholas_2013_genus.Rdata')
dat$counts <- as.matrix(dat$counts)

# ## Subset the data
# meta.list <- list()
# for(i in 1:100){
#   cat(i,'\n')
#   meta <- dat$meta.dat %>% base::split(dat$meta.dat$SubjectID) %>% sample(round(length(unique(dat$meta.dat$SubjectID)) * 0.8),replace = F) %>% bind_rows #%>% column_to_rownames('sam.ids')
#   meta.list[[i]] <- meta
# }
# save(meta.list , file = 'Data/dat.Nicholas2013.subset100.meta.Rdata')
# 
# ## Shuffle the meta data
# meta <- dat$meta.dat %>% dplyr::select(-X)
# meta.list <- list()
# for(i in 1:1000){
#   cat(i,'\n')
#   meta1 <- meta %>% rownames_to_column('sam.ids') %>%
#     group_by(SubjectID) %>% mutate(time=sample(time)) %>% as.data.frame
#   meta1 <- meta1 %>% base::split(meta1$SubjectID) %>% sample %>% bind_rows %>% column_to_rownames('sam.ids')
#   meta.list[[i]] <- meta1
# }
# save(meta.list, file = 'Data/dat.Nicholas2013.shuffle1000.meta.Rdata')



setwd('/research/bsi/projects/staff_analysis/m216453/correlated/')
source('code/CompetingMethods1.R') ## can adjust Z
load('Data/dat.Nicholas_2013_genus.Rdata')
dat$counts <- as.matrix(dat$counts)

methods_funs <- list('glmernb'='glmernb.func','lme4'='lme4.func','GLMMPQL'='GLMMPQL.func','ZIBR'='ZIBR.func','LDM'='LDM.func','NBMM'='NBMM.func',
                     'MaAsLin2'='MaAsLin2.func','IFAA'='IFAA.func',
                     'ZIGMM'='ZIGMM.func','ZINBMM'='ZINBMM.func','metamicrobiomeR'='metamicrobiomeR.func','glmmadaptive'='glmmadaptive.func',
                     'glmmTMB'='glmmTMB.func','LinDA'='LinDA.func','lme' = 'lme.func','lmer'='lmer.func')
# methods <- c('GLMMPQL','glmernb','ZIBR','LDM','LinDA','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMB','MaAsLin2')
wrapper <- match.fun(methods_funs[[method]])
if(method=='MaAsLin2'){
  wrapper.obj <- wrapper(dat,cutoff=0.05,output = '/research/bsi/projects/staff_analysis/m216453/', data.type ='matched.pair')
}else{
  wrapper.obj <- wrapper(dat,cutoff=0.05, data.type ='matched.pair')
}
save(wrapper.obj, file = paste0('/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/', method, '_Nicholas2013.res.Rdata'))







