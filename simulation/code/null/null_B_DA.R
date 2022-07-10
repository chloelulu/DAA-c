prefix <- 'B'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/correlated/null/", prefix))
source("/research/bsi/projects/staff_analysis/m216453/correlated/null/U/Cluster_mayo.R")
temp <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35.RData',envir=.GlobalEnv)
temp0 <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35_dirmult.RData',envir=.GlobalEnv)
source('/research/bsi/projects/staff_analysis/m216453/correlated/code/pre_DA.R')


paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix
paras$methods <- c('glmernb','GLMMPQL','ZIBR','LDM','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMB','LinDA')
#paras$methods <- 'MaAsLin2'
paras$nSubjects <- 100
paras$nTimes <- 2
paras$matched.pair <- F
paras$error.sds <- c(1,4)

paras$balanced.Xs <- T
paras$balanced.Ts <- T
paras$balanced.XTs <- T

paras$MgXs <- 0
paras$SgXs <- 0
paras$X.diff.otu.pcts <- 0

paras$MgTs <- 0
paras$SgTs <- 0
paras$SbTs <- 0
paras$time.diff.otu.pcts <- 0

paras$MgXTs <- 0
paras$SgXTs <- 0
paras$interaction.diff.otu.pcts <- 0

setwd(resdir)
res <- clsapply(1:1000, func, paras, queque='1-day', tempdir=file.path(resdir, 'tmpC'))
setwd(resdir)
save(res, file=file.path(resdir, paste(prefix, "_res.Rdata", sep="")))
