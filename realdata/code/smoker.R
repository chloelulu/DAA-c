# setwd('/research/bsi/projects/staff_analysis/m216453/correlated/')
# load('Data/smoker_qiita_full.RData')
# ind <- smokers$meta$AIRWAYSITE == 'Throat' #& smokers$meta$SIDEOFBODY == 'Left'
# otu.tab <- as.matrix(smokers$otu[, ind])
# otu.tab <- otu.tab[rowSums(otu.tab)>0,colSums(otu.tab) > 1000]
# meta <- as.data.frame(as.matrix(smokers$meta))[,c('SMOKER','SEX','SIDEOFBODY','HOST_SUBJECT_ID')]
# colnames(meta) <- c('X','Z','time','SubjectID')
# meta <- meta[colnames(otu.tab),]
# rownames(otu.tab) <- paste0('otu',rownames(otu.tab))
# meta$X <- as.factor(meta$X)
# meta$Z <- as.factor(meta$Z)
# meta$time <- as.factor(meta$time)
# meta$SubjectID <- as.factor(meta$SubjectID)
# meta <- meta[meta$SubjectID %in% (names(table(meta$SubjectID)[table(meta$SubjectID)>1])),]
# otu.tab <- otu.tab[,rownames(meta)]
# 
# ## Basic pre-filtering
# diff.seq.p <- cor.test(colSums(otu.tab),as.numeric(as.factor(meta$X)), method = 'spearman')$p.value
# prev = 0.1; minp = 0.002
# size.factor <- GMPR::GMPR(t(otu.tab))
# size.factor[is.na(size.factor)] <- 1
# prop <- t(t(otu.tab) / colSums(otu.tab))
# colSums(prop)
# prop <- prop[rowSums(prop!=0) > prev * ncol(prop), , drop=FALSE]
# otu.tab <- otu.tab[rownames(prop), , drop=FALSE]
# prop <- prop[matrixStats::rowMaxs(prop) > minp, , drop=FALSE]
# otu.tab <- otu.tab[rownames(prop), , drop=FALSE]
# otu.tab <- otu.tab[,colSums(otu.tab)>0]# used for analysis
# prop <- prop[,colnames(otu.tab)]
# meta <- meta[colnames(otu.tab),,drop =T] %>% droplevels()  %>% dplyr::select(X,Z,time, SubjectID)
# size.factor <- size.factor[colnames(otu.tab)]
# dat <- list(counts = otu.tab, prop = prop, meta.dat = meta, gmpr.size.factor = size.factor)
# save(dat,file = 'Data/dat.smoker_2010.Rdata')


setwd('/research/bsi/projects/staff_analysis/m216453/correlated/')
source('code/CompetingMethods1.R') ## can adjust Z
load('Data/dat.smoker_2010.Rdata')
meta <- dat$meta.dat

### Shuffle the meta data
# meta.list <- list()
# for(i in 1:1000){
#   cat(i,'\n')
#   meta1 <- meta %>% #rownames_to_column('sam.ids') %>%
#     group_by(SubjectID) %>% mutate(X=sample(X)) %>% as.data.frame
#   meta1 <- meta1 %>% base::split(meta1$SubjectID) %>% sample %>% bind_rows #%>% column_to_rownames('sam.ids')
#   rownames(meta1) <- rownames(meta)
#   # meta1 <- meta1 %>% base::split(meta1$SubjectID) %>% sample %>% base::unsplit(meta1$SubjectID)
#   # meta1$X.old <- meta$X
#   meta.list[[i]] <- meta1
# }
# save(meta.list, file = 'Data/dat.smoker_2010.shuffle1000.meta.Rdata')

# load('Data/dat.smoker_2010.Rdata')
# # prepare subset meta.data
# meta.list <- list()
# for(i in 1:100){
#   cat(i,'\n')
#   meta <- dat$meta.dat %>% base::split(dat$meta.dat$SubjectID) %>% sample(round(length(unique(dat$meta.dat$SubjectID)) * 0.8),replace = F) %>% bind_rows #%>% column_to_rownames('sam.ids')
#   meta.list[[i]] <- meta
# }
# save(meta.list , file = 'Data/dat.smoker.subset100.meta.Rdata')



## Figure 1: Detection number vs FDR
## Overlap of findings among datasets
## Null scenario
## Random subset of samples, which methods shows higher stability (if different subset show totally different result, signifying the non-stability)

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  method = args[1]
  method = as.character(method)}


setwd('/research/bsi/projects/staff_analysis/m216453/correlated/')
source('code/CompetingMethods1.R') ## can adjust Z
load('Data/dat.smoker_2010.Rdata')
dat$meta.dat$X <- as.numeric(dat$meta.dat$X)-1
dat$meta.dat$time <- as.numeric(dat$meta.dat$time)-1


y <- (dat$meta.dat[dat$meta.dat$X=='y',]);length(unique(y$SubjectID))
n <- (dat$meta.dat[dat$meta.dat$X=='n',]);length(unique(n$SubjectID))

methods_funs <- list('glmernb'='glmernb.func','lme4'='lme4.func','GLMMPQL'='GLMMPQL.func','ZIBR'='ZIBR.func','LDM'='LDM.func','NBMM'='NBMM.func',
                     'MaAsLin2'='MaAsLin2.func','IFAA'='IFAA.func',
                     'ZIGMM'='ZIGMM.func','ZINBMM'='ZINBMM.func','metamicrobiomeR'='metamicrobiomeR.func','glmmadaptive'='glmmadaptive.func',
                     'glmmTMB'='glmmTMB.func','LinDA'='LinDA.func','lme' = 'lme.func','lmer'='lmer.func')
# methods <- c('GLMMPQL','glmernb','ZIBR','LDM','LinDA','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMB','MaAsLin2')

wrapper <- match.fun(methods_funs[[method]])
if(method=='MaAsLin2'){
  wrapper.obj <- wrapper(dat,cutoff=0.05,output = '/research/bsi/projects/staff_analysis/m216453/', data.type ='repeated.measure')
}else{
  wrapper.obj <- wrapper(dat,cutoff=0.05, data.type ='repeated.measure')
}
save(wrapper.obj, file = paste0('/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/', method, '_smoker.res.Rdata'))



