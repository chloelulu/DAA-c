prefix <- 'er4'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/correlated/smallsampleU/", prefix))
source("/research/bsi/projects/staff_analysis/m216453/correlated/code/Cluster_mayo.R")
temp <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35.RData',envir=.GlobalEnv)
temp0 <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35_dirmult.RData',envir=.GlobalEnv)
source('/research/bsi/projects/staff_analysis/m216453/correlated/code/pre_DA.R')


paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix
paras$methods <- c('MaAsLin2','glmernb','GLMMPQL','ZIBR','LDM','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMBP','LinDA')
#paras$methods <-c('MaAsLin2')

paras$nSubjects <- 20
paras$nTimes <- 2
paras$matched.pair <- F
paras$error.sds <- 4

paras$balanced.Xs <- F
paras$balanced.Ts <- F
paras$balanced.XTs <- F

#paras$MgXs <-  c(8,10)
paras$MgXs <-  c(2.5,3.5)
paras$SgXs <- 0
paras$X.diff.otu.pcts <- 0.1

paras$MgTs <- 0
paras$SgTs <- 0
paras$SbTs <- 0
paras$time.diff.otu.pcts <- 0

paras$MgXTs <- 0
paras$SgXTs <- 0
paras$interaction.diff.otu.pcts <- 0

setwd(resdir)
res <- clsapply(1:100, func, paras, queque='1-day', tempdir=file.path(resdir, 'tmpC'))
setwd(resdir)
save(res, file=file.path(resdir, paste(prefix, "_res.Rdata", sep="")))
