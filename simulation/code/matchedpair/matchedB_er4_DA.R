prefix <- 'er4'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/correlated/matchedB/", prefix))
source("/research/bsi/projects/staff_analysis/m216453/correlated/code/Cluster_mayo.R")
temp <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35.RData',envir=.GlobalEnv)
temp0 <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35_dirmult.RData',envir=.GlobalEnv)
source('/research/bsi/projects/staff_analysis/m216453/correlated/code/pre_DA_T.R')

paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix
paras$methods <- c('MaAsLin2','lme','glmernb','GLMMPQL','ZIBR','LDM','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMB','LinDA','MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR')

paras$nSubjects <- 100
paras$nTimes <- 2
paras$matched.pair <- T
paras$error.sds <- 4

paras$balanced.Xs <- T
paras$balanced.Ts <- T
paras$balanced.XTs <- T

paras$MgXs <- 0
paras$SgXs <- 0
paras$X.diff.otu.pcts <- 0

paras$MgTs <- c(2.5,3.5)
paras$SgTs <- 0
paras$SbTs <- 0
paras$time.diff.otu.pcts <- 0.1

paras$MgXTs <- 0
paras$SgXTs <- 0
paras$interaction.diff.otu.pcts <- 0


setwd(resdir)
res <- clsapply(1:100, func, paras, queque='1-day', tempdir=file.path(resdir, 'tmpC'))
setwd(resdir)
save(res, file=file.path(resdir, paste(prefix, "_res.Rdata", sep="")))
