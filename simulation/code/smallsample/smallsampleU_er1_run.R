prefix <- 'smallsampleU/er1/'
resdir <-  file.path(paste0("/research/bsi/projects/staff_analysis/m216453/correlated/", prefix))
temp <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35.RData',envir=.GlobalEnv)
temp0 <- load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35_dirmult.RData',envir=.GlobalEnv)
source('/research/bsi/projects/staff_analysis/m216453/correlated/code/pre_run.R')

paras <- list()
paras$model.paras <- NULL
paras$resdir <- resdir
paras$prefix <- prefix

paras$nSubjects <- 20
paras$nTimes <- 2
paras$matched.pair <- F
paras$error.sds <- 1

paras$balanced.Xs <- F
paras$balanced.Ts <- F
paras$balanced.XTs <- F

#paras$MgXs <- c(4,6)
paras$MgXs <- c(1,2)
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
clsapply(1:100, func, paras, queque='1-hour', tempdir=file.path(resdir, 'tmpC'))
