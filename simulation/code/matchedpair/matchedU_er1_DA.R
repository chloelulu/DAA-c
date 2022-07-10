prefix <- 'er1'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/correlated/matchedU/", prefix))
source("/research/bsi/projects/staff_analysis/m216453/correlated/code/Cluster_mayo.R")
temp <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35.RData',envir=.GlobalEnv)
temp0 <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35_dirmult.RData',envir=.GlobalEnv)
source('/research/bsi/projects/staff_analysis/m216453/correlated/code/pre_DA_T.R')

paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix
#paras$methods <- c('MaAsLin2','lme','glmernb','GLMMPQL','ZIBR','LDM','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMB','LinDA')
paras$methods <- 'IFAA'#c('MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR')

paras$nSubjects <- 100
paras$nTimes <- 2
paras$matched.pair <- T
paras$error.sds <- 1

paras$balanced.Xs <- F
paras$balanced.Ts <- F
paras$balanced.XTs <- F

paras$MgXs <- 0
paras$SgXs <- 0
paras$X.diff.otu.pcts <- 0

paras$MgTs <- c(0.55,0.8)
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
