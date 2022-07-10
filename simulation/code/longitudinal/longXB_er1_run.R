prefix <- 'longXB/er1/'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/correlated/", prefix))
source("/research/bsi/projects/staff_analysis/m216453/correlated/code/Cluster_mayo.R")
temp <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35.RData',envir=.GlobalEnv)
temp0 <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35_dirmult.RData',envir=.GlobalEnv)
source('/research/bsi/projects/staff_analysis/m216453/correlated/code/pre_run.R')

paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix

paras$nSubjects <- 40
paras$nTimes <- 5
paras$matched.pair <- FALSE
paras$error.sds <- 1

paras$balanced.Xs <- T
paras$balanced.Ts <- T
paras$balanced.XTs <- T

paras$MgXs <- c(2,3)
paras$SgXs <- 0
paras$X.diff.otu.pcts <- 0.1

paras$MgTs <- 0.5
paras$SgTs <- 0
paras$SbTs <- 0.5
paras$time.diff.otu.pcts <- 0.1

paras$MgXTs <- 0
paras$SgXTs <- 0
paras$interaction.diff.otu.pcts <- 0


setwd(resdir)
clsapply(1:100, func, paras, queque='1-hour', tempdir=file.path(resdir, 'tmpC'))
