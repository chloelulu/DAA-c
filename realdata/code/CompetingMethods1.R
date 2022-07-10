## This script is different from CompetingMethods.R in following aspects:
## 1. no FDR calculation
## 2. GLMMPQL, glmadaptive, glmTMB.. needs to extract each OTU, they needd to be t() in CompetingMethods0.R for real data

## This the customized script for real datasets. matched pair data (Nicholas) does not have Z(confounder), reapted and longitudinal have Z
cal.FDR.TPR <- function(var.cols, data, alpha=0.05){
  res.obj <- list()
  for(var.col in var.cols){
    tp <- sum(data[,paste0(var.col,'.padj')] <= alpha & data[,paste0(var.col,'.truth')] ==TRUE, na.rm = T)
    tn <- sum(data[,paste0(var.col,'.padj')] > alpha & data[,paste0(var.col,'.truth')] ==FALSE, na.rm = T)
    fn <- sum(data[,paste0(var.col,'.padj')] > alpha &  data[,paste0(var.col,'.truth')] ==TRUE, na.rm = T)
    fp <- sum(data[,paste0(var.col,'.padj')] <= alpha &  data[,paste0(var.col,'.truth')] ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    obj <- c(tpr, fdr, fpr, fp, F1, MCC)
    names(obj) <- c('tpr', 'fdr', 'fpr', 'fp', 'F1', 'MCC')
    res.obj[[var.col]] <- obj
  }
  return(res.obj)
}

# source('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/code/SimulationEvaluation/ancom_v2.1.R')
pkg <- c('nlme','lme4','MASS','LinDA','ZIBR','LDM','aod','compositions','tidyverse','GMPR','glmmTMB','NBZIMM','GLMMadaptive','tibble','Maaslin2','IFAA')
suppressPackageStartupMessages(sapply(pkg, require, character = T))
# source('~/Documents/Mayo_Research/2021_03_29_CorrelatedData/Code/ancom2.R')

## For current test realdata, we dont have Z for matched pair

## lme4
## http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
glmernb.func <- function(dat, cutoff = 0.05, data.type='longitudinal'){
  otu.tab=dat[['counts']]
  meta=dat[['meta.dat']]
  t_start = Sys.time()
  pvs <- name <- NULL
  for (y in rownames(otu.tab)){
    data <- merge(merge(t(otu.tab), meta, by = 0) %>% column_to_rownames('Row.names'), cbind.data.frame(log.size.factor = log(dat$gmpr.size.factor)), by = 0)
    model <- NULL
    if(data.type =='repeated.measure'){
      colnames(data)[grep(y,colnames(data))] <- 'taxon'
      try(model <- glmer.nb(taxon ~ X + Z +(1|SubjectID)+ offset(log.size.factor), control=glmerControl(nAGQ0initStep=FALSE),data = data),silent = TRUE)
      if(class(model) %in% 'NULL'){pv <- NA}else{pv <- coef(summary(model))[,4][2]}
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) <- 'X'
    }
    
    if(data.type =='matched.pair'){
      try(model <- glmer.nb(formula(paste(y, '~ time+(1|SubjectID)+ offset(log.size.factor)')), control=glmerControl(nAGQ0initStep=FALSE),data = data),silent = TRUE)
      if(class(model) %in% 'NULL'){pv <- NA}else{pv <- coef(summary(model))[,4][2]}
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) <- 'time'
    }
    
    if(data.type =='longitudinal'){ ## needs to be checked
      try(model <- lme4::glmer.nb(formula(paste(y, "~ X*time+Z+ (time|SubjectID) + offset(log.size.factor)")), 
                                  control=glmerControl(nAGQ0initStep=FALSE),data = data),silent = TRUE)
      if(class(model) %in% 'NULL'){pv <- c(NA,NA,NA)}else{pv <- coef(summary(model))[,4][c(2:3,5)]}
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs)[c(1:3)] <- c('X','time','interaction')
    }
  }
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- as.data.frame(apply(pvs, 2, function(x) p.adjust(x,method ='fdr')))
  rownames(res) <- name
  colnames(res) = paste0(colnames(pvs),'.padj')
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  
  return(list(res=res, na = na, time = t_run, sig = sig))
}


## GLMMPQL
GLMMPQL.func <- function(dat, cutoff = 0.05,data.type='longitudinal'){
  otu.tab <- dat$counts
  meta <- dat$meta.dat
  t_start = Sys.time()
  pvs <- name <- NULL
  for (y in rownames(otu.tab)){
    df <- merge(merge(cbind(taxon = otu.tab[y,]),meta,by = 0) %>% column_to_rownames('Row.names'), cbind(log.size.factor = log(dat$gmpr.size.factor)),by = 0) %>% column_to_rownames('Row.names')
    control = lmeControl(maxIter = 1e8, msMaxIter = 1e8)
    model <- NULL
    if(data.type =='repeated.measure'){
      try(model <- MASS::glmmPQL(taxon~ X + Z + offset(log.size.factor), random = ~ 1 | SubjectID, data = df, verbose=F, family=quasipoisson(), control = control),silent = T)
      if(class(model) %in% 'NULL'){pv <- NA}else{pv <- coef(summary(model))[,5][2]}
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'X'
    }
    
    if(data.type =='matched.pair'){
      try(model <-  MASS::glmmPQL(taxon~ time +offset(log.size.factor), random = ~ 1 | SubjectID, data = df, verbose=F, family=quasipoisson(), control = control),silent = T)
      if(class(model) %in% 'NULL'){pv <- NA}else{pv <- coef(summary(model))[,5][2]}
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'time'
    }
    
    if(data.type =='longitudinal'){
      try(model <-  MASS::glmmPQL(taxon~ X*time+Z+offset(log.size.factor), random = ~ time | SubjectID, data = df, verbose=F, family=quasipoisson(), control = control),silent = T)
      if(class(model) %in% 'NULL'){pv <- c(NA,NA,NA)}else{pv <- coef(summary(model))[,5][c(2:3,5)]}
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs)[c(1:3)] = c('X','time','interaction')
    }
  }
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- apply(pvs, 2, function(x) p.adjust(x,method ='fdr'))
  cat(nrow(res),'\n')
  res <- as.data.frame(res)
  rownames(res) <- name
  colnames(res) <- paste0(colnames(res),'.padj')
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  res <- as.data.frame(res)
  return(list(res=res, na = na, time = t_run, sig = sig))
}



## ZIBR https://github.com/chvlyl/ZIBR
##  each subject must have the same number of time points. 
ZIBR.func <- function(dat, cutoff = 0.05, data.type ='longitudinal'){
  t_start = Sys.time()
  prop = dat[['prop']]
  meta = dat[['meta.dat']]
  time.ind <- meta$time
  subject.ind <- as.character(meta[,'SubjectID'])
  pvs <- name <- NULL
  for(i in 1:nrow(prop)){
    Y <- prop[i,]
    est <- NULL
    if(data.type=='repeated.measure'){
      X <- meta[,c('X'),drop=F]
      try(est <- zibr(logistic.cov=X,beta.cov=X,Y=Y,subject.ind=subject.ind,time.ind=time.ind),silent = T)
      if(class(est)=='NULL'){pv = NA}else{pv <- est$joint.p['X']}
      name <- c(name, rownames(prop)[i])
      pvs <- rbind(pvs, pv) 
    }
    
    if(data.type=='matched.pair'){
      X <- data.frame(meta[,c('time'),drop=F])
      try(est <- zibr(logistic.cov=X,beta.cov=X,Y=Y,subject.ind=subject.ind,time.ind=time.ind),silent = T)
      if(class(est)=='NULL'){pv = NA}else{pv <- est$joint.p['time']}
      name <- c(name, rownames(prop)[i])
      pvs <- rbind(pvs, pv) 
    }
    
    if(data.type =='longitudinal'){
      X <- meta[,c('X','time'),drop=F]
      X$interaction <- X$X * X$time
      try(est <- zibr(logistic.cov=X,beta.cov=X,Y=Y,subject.ind=subject.ind,time.ind=time.ind),silent = T)
      if(class(est)=='NULL'){pv = c(NA,NA,NA)}else{pv <- est$joint.p[c('X','time','interaction')]}
      name <- c(name, rownames(prop)[i])
      pvs <- rbind(pvs, pv) 
    }
  }
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- as.data.frame(apply(pvs, 2, function(x) p.adjust(x,method ='fdr')))
  rownames(res) <- name
  colnames(res) = paste0(colnames(res),'.padj')
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  res <- as.data.frame(res)
  return(list(res=res, na = na, time = t_run, sig = sig))
}

## LDM 
LDM.func <- function(dat, cutoff = 0.05, data.type ='longitudinal'){
  otu.tab <<- dat[['counts']]
  meta <- dat[['meta.dat']]
  # meta$SubjectID <- as.numeric(meta$SubjectID)#  LDM tutorial demo data is numeric for ID, while it will change into factor in matched case! It should be numeric here, otherwise will give error
  t_start = Sys.time()
  if(data.type == 'repeated.measure'){
    res.ldm <- ldm(formula=otu.tab|Z ~ X, data=meta, seed=34794, cluster.id = SubjectID, perm.within.type="none", perm.between.type="free", fdr.nominal=cutoff)
    res <- as.data.frame(apply(t(res.ldm$p.otu.omni), 2, function(x) p.adjust(x, method = 'fdr')))
    colnames(res) <- c('X.padj')
    t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
    na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
    sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  }
  
  if(data.type == 'matched.pair'){
    res.ldm <- ldm(formula=otu.tab |factor(SubjectID) ~ time, data=meta, cluster.id=SubjectID, perm.within.type="free", perm.between.type="none", fdr.nominal=cutoff, n.perm.max = 5000)
    res <- as.data.frame(apply(t(res.ldm$p.otu.omni), 2, function(x) p.adjust(x, method = 'fdr')))
    colnames(res) <- c('time.padj')
    t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
    na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
    sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  }
  
  if(data.type == 'longitudinal'){
    ## interaction (SubjectID should be equal size in each subject(Note that if perm.between.type has a value other than “none” then all clusters must have the same size.))
    res.ldm1 <- ldm(formula=otu.tab|Z ~ time + (X:time), data=meta, 
                    cluster.id=SubjectID, perm.within.type="free", perm.between.type="none", fdr.nominal=cutoff)
    
    res.ldm2 <- ldm(formula=otu.tab|Z ~ X, data=meta,
                    cluster.id=SubjectID, perm.within.type='none', perm.between.type="free", fdr.nominal=cutoff)
    
    res1 <- as.data.frame(apply(t(res.ldm1$p.otu.omni), 2, function(x) p.adjust(x, method = 'fdr')))
    colnames(res1) <- c('time.padj','interaction.padj')
    res2 <- as.data.frame(apply(t(res.ldm2$p.otu.omni), 2, function(x) p.adjust(x, method = 'fdr')))
    colnames(res2) <- c('X.padj')
    
    res <- merge(res1, res2, by = 0) %>% column_to_rownames('Row.names')
    res <- res[,c('X.padj','time.padj','interaction.padj'),drop =F]
    t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
    na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
    sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  }
  
  return(list(res=res, na = na, time = t_run, sig = sig))
}


## https://drizopoulos.github.io/GLMMadaptive/articles/ZeroInflated_and_TwoPart_Models.html#zero-inflated-poisson-mixed-effects-model-1
## encounter potential error, see https://stats.stackexchange.com/questions/397955/mixed-effect-zero-inflated-negative-binomial-model-the-leading-minor-of-order
glmmadaptive.func <- function(dat, cutoff = 0.05,data.type ='longitudinal'){
  t_start = Sys.time()
  otu.tab = (dat$counts)
  meta = dat$meta.dat
  pvs <- name <- NULL
  for (i in 1:nrow(otu.tab)){
    df <- merge(merge(cbind(taxon = otu.tab[i,]),meta,by = 0) %>% column_to_rownames('Row.names'), cbind(log.size.factor = log(dat$gmpr.size.factor)),by = 0) %>% column_to_rownames('Row.names')
    gm1 <- NULL
    if(data.type =='repeated.measure'){
      try(gm1 <- mixed_model(fixed = taxon ~ X +Z+ offset(log.size.factor), random = ~ 1|SubjectID, data = df, family=zi.negative.binomial(), zi_fixed = ~ 1,iter_EM = 0),silent = T)# default is 30, iter_EM = 0, When you set it to zero it means that you're only using the quasi-Newton algorithm
      if(class(gm1) =='NULL'){pv <- NA}else{pv <- summary(gm1)$coef_table[2,4]}  
      name <- c(name,rownames(otu.tab)[i])
      pvs <- rbind(pvs,pv)
      colnames(pvs) <- 'X'
    }
    
    if(data.type =='matched.pair') {
      try(gm1 <- mixed_model(fixed = taxon ~ time +offset(log.size.factor), random = ~ 1|SubjectID, data = df, family=zi.negative.binomial(), zi_fixed = ~ 1,iter_EM = 0),silent = T)
      if(class(gm1) =='NULL'){pv <- NA}else{pv <- summary(gm1)$coef_table[2,4]}  
      name <- c(name,rownames(otu.tab)[i])
      pvs <- rbind(pvs,pv)
      colnames(pvs) <- 'time'
    }
    
    if(data.type =='longitudinal') {
      try(gm1 <- mixed_model(fixed = taxon ~ X * time +Z+ offset(log.size.factor), random = ~ time|SubjectID, data = df, family=zi.negative.binomial(), zi_fixed = ~ 1,iter_EM = 0),silent = T)
      if(class(gm1) =='NULL'){pv <-c(NA,NA,NA)}else{pv <- summary(gm1)$coef_table[,4][c(2,3,5)]}
      name <- c(name,rownames(otu.tab)[i])
      pvs <- rbind(pvs,pv)
      colnames(pvs) <- c('X','time','interaction')
    }
  }
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- as.data.frame(apply(pvs, 2, function(x) p.adjust(x,method ='fdr')))
  colnames(res) <- paste0(colnames(res),'.padj')
  rownames(res) = name
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  return(list(res=res, na = na, time = t_run, sig = sig))
}


## https://fromthebottomoftheheap.net/2017/05/04/compare-mgcv-with-glmmtmb/
## the conditional output represents the zero portion (or a logistic regression) - the zero inflated output represents a "mixture" model of the two distributions - one for the subgroup who reports zero or close to zero and one for the subgroup who doesn't report zero.
glmmTMB.func <- function(dat, cutoff = 0.05,data.type ='longitudinal'){
  t_start = Sys.time()
  otu.tab = dat$counts
  meta = dat$meta.dat
  pvs <- name <- NULL
  for (i in 1:nrow(otu.tab)){
    df <- merge(merge(cbind(taxon = otu.tab[i,]),meta,by = 0) %>% column_to_rownames('Row.names'), cbind(log.size.factor = log(dat$gmpr.size.factor)),by = 0) %>% column_to_rownames('Row.names')
    gm1 <- NULL
    if(data.type =='repeated.measure'){
      try(gm1 <- glmmTMB::glmmTMB(formula = taxon ~ X + Z + (1| SubjectID)+offset(log.size.factor), zi=  ~ 1, family=nbinom2, data=df),silent = T)
      if(class(gm1) =='NULL'){pv <- NA}else{pv <- summary(gm1)$coefficients$cond[2,4]}
      name <- c(name,rownames(otu.tab)[i])
      pvs <- rbind(pvs,pv)
      colnames(pvs) <- 'X'}
    
    if(data.type =='matched.pair') {
      try(gm1 <- glmmTMB::glmmTMB(formula = taxon ~ time +  (1| SubjectID )+offset(log.size.factor), zi=  ~ 1,family=nbinom2, data=df),silent = T)
      if(class(gm1) =='NULL'){pv <- NA}else{pv <- summary(gm1)$coefficients$cond[2,4]}
      name <- c(name,rownames(otu.tab)[i])
      pvs <- rbind(pvs,pv)
      colnames(pvs) <- 'time'
    }
    
    if(data.type =='longitudinal') {
      try(gm1 <- glmmTMB::glmmTMB(formula = taxon  ~  X*time + Z+ (time| SubjectID)+offset(log.size.factor), zi=  ~ 1,family=nbinom2, data=df),silent = T)
      if(class(gm1) =='NULL'){pv <- c(NA,NA,NA)}else{pv <- summary(gm1)$coefficients$cond[,4][c(2,3,5)]}
      name <- c(name,rownames(otu.tab)[i])
      pvs <- rbind(pvs,pv)
      colnames(pvs) <- c('X','time','interaction')
    }
  }
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- as.data.frame(apply(pvs, 2, function(x) p.adjust(x,method ='fdr')))
  colnames(res) <- paste0(colnames(res),'.padj')
  rownames(res) = name
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  res <- as.data.frame(res)
  return(list(res=res, na = na, time = t_run, sig = sig))
}


## LinDA https://github.com/zhouhj1994/LinDA
LinDA.func <- function(dat, cutoff = 0.05, data.type ='longitudinal'){
  t_start = Sys.time()
  otu.tab = dat$counts
  meta = dat$meta.dat
  if(data.type =='repeated.measure'){
    linda.obj <- LinDA::linda(otu.tab, meta, formula = '~X+Z+(1|SubjectID)', alpha = cutoff, prev.cut = 0, lib.cut = 0, winsor.quan = NULL)
    res <- linda.obj$output[[1]][,c('padj'), drop = F]
    colnames(res) <- 'X.padj'
  }
  if(data.type =='matched.pair'){
    linda.obj <- LinDA::linda(otu.tab, meta, formula = '~time+(1|SubjectID)', alpha = cutoff, prev.cut = 0, lib.cut = 0, winsor.quan = NULL)
    res <- linda.obj$output[[1]][,c('padj'), drop = F]
    colnames(res) <- 'time.padj'
  } 
  if(data.type =='longitudinal'){
    linda.obj <- GUniFrac::linda(otu.tab, meta, formula = '~X*time+Z+(time|SubjectID)', 
                                 feature.dat.type = 'count')
    res <- merge(linda.obj$output$X[,c('padj'), drop = F] %>% dplyr::rename(X.padj = 'padj'), 
                 linda.obj$output$time[,c('padj'), drop = F]%>% dplyr::rename(time.padj = 'padj'),all.x =T, all.y = T, by = 0) %>% 
      column_to_rownames('Row.names') %>% merge(linda.obj$output$`X:time`[,c('padj'), drop = F]%>% dplyr::rename(interaction.padj = 'padj'),all.x =T, all.y = T, by = 0) %>% column_to_rownames('Row.names') 
  } 
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  res <- as.data.frame(res) 
  return(list(res=res, na = na, time = t_run, sig = sig))
}


lmer.func <- function(dat, cutoff = 0.05,data.type='longitudinal'){
  prop = asin(sqrt(dat[['prop']]))
  meta = dat[['meta.dat']]
  t_start = Sys.time()
  pvs <- name <- NULL
  for (y in rownames(prop)){
    df <- merge(merge(cbind(taxon = otu.tab[y,]),meta,by = 0) %>% column_to_rownames('Row.names'), cbind(log.size.factor = log(dat$gmpr.size.factor)),by = 0) %>% column_to_rownames('Row.names')
    model <- NULL
    if(data.type =='repeated.measure'){
      try(model <- lme4::lmer(taxon ~ X +Z+  (1|SubjectID) + offset(log.size.factor), data=df),silent = T)
      if(class(model) %in% 'NULL'){pv <- NA}else{
        pv <- wald.test(b = fixed.effects(model), Sigma = vcov( model), Terms = grep('^X$',names(fixed.effects( model))))$result$chi2['P']
      }
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'X'
    }
    
    if(data.type =='matched.pair'){
      try(model <- lme4::lmer(taxon ~ X +  (1|SubjectID) + offset(log.size.factor), data=df),silent = T)
      if(class(model) %in% 'NULL'){pv <- NA}else{
        pv <- wald.test(b = fixed.effects(model), Sigma = vcov( model), Terms = grep('^time$',names(fixed.effects( model))))$result$chi2['P']
      }
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'X'
    }
    
    if(data.type =='longitudinal'){
      try(model <- lme4::lmer(taxon ~ X*time +  (time|SubjectID) + offset(log.size.factor), data=df),silent = T)
      if(class(model) %in% 'NULL'){pv <- c(NA,NA,NA)}else{
        pv.X <- wald.test(b = fixed.effects( model), Sigma = vcov( model), Terms = grep('^X$',names(fixed.effects(model))))$result$chi2['P']
        pv.T <- wald.test(b = fixed.effects( model), Sigma = vcov( model), Terms = grep('^time$',names(fixed.effects(model))))$result$chi2['P']
        pv.XT <- wald.test(b = fixed.effects( model), Sigma = vcov( model), Terms = grep('X:time',names(fixed.effects(model))))$result$chi2['P']
        pv <- c(pv.X, pv.T, pv.XT)
      }
      names(pv) = c('X','time','interaction')
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
    }
  }
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- as.data.frame(apply(pvs, 2, function(x) p.adjust(x,method ='fdr')))
  rownames(res) <- name
  colnames(res) = paste0(colnames(res),'.padj')
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  return(list(res=res, na = na, time = t_run, sig = sig))
}



lme.func <- function(dat, cutoff = 0.05,data.type='longitudinal'){
  prop = asin(sqrt(dat[['prop']]))
  meta = dat[['meta.dat']]
  t_start = Sys.time()
  pvs <- name <- NULL
  for (y in rownames(prop)){
    df <- merge(merge(cbind(taxon = otu.tab[y,]),meta,by = 0) %>% column_to_rownames('Row.names'), cbind(log.size.factor = log(dat$gmpr.size.factor)),by = 0) %>% column_to_rownames('Row.names')
    model <- NULL
    if(data.type =='repeated.measure'){
      try(model <- lme(taxon ~ X +  offset(log.size.factor), random = ~1|SubjectID, data=df),silent = T)
      if(class(model) %in% 'NULL'){pv = NA}else{pv <- summary(model)$tTable[2, 5]}
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'X'
    }
    
    if(data.type =='matched.pair'){
      try(model <- lme(taxon ~ time +  offset(log.size.factor), random = ~1|SubjectID, data=df),silent = T)
      if(class(model) %in% 'NULL'){pv = NA}else{pv <- summary(model)$tTable[2, 5]}
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'X'
    }
    
    if(data.type =='longitudinal'){
      try(model <- lme(taxon ~ X*time +  offset(log.size.factor), list(SubjectID = pdDiag(~time)), data=df),silent = T)
      if(class(model) %in% 'NULL'){pv = NA}else{pv <- summary(model)$tTable[c(2:3,5), 5]}
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = c('X','time','interaction')
    }
  }
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- as.data.frame(apply(pvs, 2, function(x) p.adjust(x,method ='fdr')))
  rownames(res) <- name
  colnames(res) = paste0(colnames(res),'.padj')
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  return(list(res=res, na = na, time = t_run, sig = sig))
}


NBMM.func <- function(dat, cutoff = 0.05, data.type ='longitudinal'){
  otu.tab = dat$counts
  meta = dat[['meta.dat']]
  t_start = Sys.time()
  pvs <- name <- NULL
  for (y in rownames(otu.tab)){
    df <- merge(merge(cbind(taxon = otu.tab[y,]),meta,by = 0) %>% column_to_rownames('Row.names'), cbind(log.size.factor = log(dat$gmpr.size.factor)),by = 0) %>% column_to_rownames('Row.names')
    model <- NULL
    if(data.type =='repeated.measure'){
      try(model <- glmm.nb(taxon ~ X +Z+  offset(log.size.factor), random = ~1|SubjectID, data=df),silent = T)
      if(class(model) %in% 'NULL'){pv <- NA}else{
        pv <- summary(model)$tTable[2, 5]
      }
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'X'
    }
    
    if(data.type =='matched.pair'){
      try(model <- glmm.nb(taxon ~ time + offset(log.size.factor), random = ~1|SubjectID, data=df),silent = T)
      if(class(model) %in% 'NULL'){
        pv <- NA
      }else{
        pv <- summary(model)$tTable[2, 5]
      }
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'time'
    }
    
    if(data.type =='longitudinal'){
      try(model <- glmm.nb(taxon ~ X*time+Z+offset(log.size.factor), random = list(SubjectID = pdDiag(~time)), data=df),silent = T)
      if(class(model) %in% 'NULL'){
        pv <- c(NA,NA,NA) 
      }else{ 
        pv <- summary(model)$tTable[c(2:3,5), 5]
      }
      names(pv) = c('X','time','interaction')
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
    }
  }
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- as.data.frame(apply(pvs, 2, function(x) p.adjust(x,method ='fdr')))
  rownames(res) <- name
  colnames(res) = paste0(colnames(res),'.padj')
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  return(list(res=res, na = na, time = t_run, sig = sig))
}


ZIGMM.func <- function(dat, cutoff = 0.05, data.type ='longitudinal'){
  otu.tab=log(dat$counts+1)
  meta=dat[['meta.dat']]
  t_start = Sys.time()
  pvs <- name <- NULL
  for (y in rownames(otu.tab)){
    df <- merge(merge(cbind(taxon = otu.tab[y,]),meta,by = 0) %>% column_to_rownames('Row.names'), cbind(log.size.factor = log(dat$gmpr.size.factor)),by = 0) %>% column_to_rownames('Row.names')
    model <- NULL
    if(data.type =='repeated.measure'){
      try(model <- lme.zig(taxon ~ X +Z+  offset(log.size.factor), random = ~1|SubjectID, data=df),silent = T)
      if(class(model) %in% 'NULL'){pv <- NA}else{
        pv <- summary(model)$tTable[2, 5]
      }
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'X'
    }
    
    if(data.type =='matched.pair'){
      try(model <- lme.zig(taxon ~ time +offset(log.size.factor), random = ~1|SubjectID, data=df),silent = T)
      if(class(model) %in% 'NULL'){pv <- NA}else{
        pv <- summary(model)$tTable[2, 5]
      }
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'time'
    }
    
    if(data.type =='longitudinal'){
      try(model <- lme.zig(taxon ~ X*time +  Z+offset(log.size.factor), random = list(SubjectID = pdDiag(~time)), data=df),silent = T)
      if(class(model) %in% 'NULL'){ pv <- c(NA,NA,NA) }else{ pv <- summary(model)$tTable[c(2:3,5), 5]}
      names(pv) = c('X','time','interaction')
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
    }
  }
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- as.data.frame(apply(pvs, 2, function(x) p.adjust(x,method ='fdr')))
  rownames(res) <- name
  colnames(res) = paste0(colnames(res),'.padj')
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  return(list(res=res, na = na, time = t_run, sig = sig))
}



ZINBMM.func <- function(dat, cutoff = 0.05, data.type ='longitudinal'){
  otu.tab = dat$counts
  meta = dat[['meta.dat']]
  t_start = Sys.time()
  pvs <- name <- NULL
  for (y in rownames(otu.tab)){
    df <- merge(merge(cbind(taxon = otu.tab[y,]),meta,by = 0) %>% column_to_rownames('Row.names'), cbind(log.size.factor = log(dat$gmpr.size.factor)),by = 0) %>% column_to_rownames('Row.names')
    model <- NULL
    if(data.type =='repeated.measure'){
      try(model <- glmm.zinb(taxon ~ X +Z+  offset(log.size.factor), random = ~1|SubjectID, data=df),silent = T)
      if(class(model) %in% 'NULL'){ pv <- NA }else{
        pv <- summary(model)$tTable[2, 5]
      }
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'X'
    }
    
    if(data.type =='matched.pair'){
      try(model <- glmm.zinb(taxon ~ time + offset(log.size.factor), random = ~1|SubjectID, data=df),silent = T)
      if(class(model) %in% 'NULL'){ pv <- NA }else{
        pv <- summary(model)$tTable[2, 5]
      }
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
      colnames(pvs) = 'time'
    }
    
    if(data.type =='longitudinal'){
      try(model <- glmm.zinb(taxon ~ X*time +  Z+offset(log.size.factor), random = list(SubjectID = pdDiag(~time)), data=df),silent = T)
      if(class(model) %in% 'NULL'){ pv <- c(NA,NA,NA) }else{ pv <- summary(model)$tTable[c(2:3,5), 5]}
      names(pv) = c('X','time','interaction')
      pvs <- rbind(pvs,pv)
      name <- c(name, y)
    }
  }
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  res <- as.data.frame(apply(pvs, 2, function(x) p.adjust(x,method ='fdr')))
  rownames(res) <- name
  colnames(res) = paste0(colnames(res),'.padj')
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  return(list(res=res, na = na, time = t_run, sig = sig))
}



MaAsLin2.func <- function(dat, cutoff = 0.05, data.type ='longitudinal', output='~/Desktop/'){
  t_start = Sys.time()
  otu <- t(dat$counts)
  meta <- dat$meta.dat
  if(length(unique(meta$X))==2){meta$X <- as.factor(meta$X)}
  
  if(data.type=='repeated.measure'){
    fit <- Maaslin2(input_data = otu, input_metadata = meta, output = output,min_prevalence = 0,
                   fixed_effects = c("X","Z"),random_effects ='SubjectID',plot_heatmap = F, plot_scatter = F)
    res <- fit$results[,c('feature','pval','metadata'), drop =F] %>% dplyr::filter(metadata =='X') %>% column_to_rownames('feature') %>% dplyr::select(-metadata)
    res$X.padj <- p.adjust(res$pval, method = 'fdr') 
    res <- res %>% dplyr::select(-pval)  
  }
  
  if(data.type=='matched.pair'){
    fit = Maaslin2(input_data = otu, input_metadata = meta, output = output,min_prevalence = 0,
                   fixed_effects = c("time"),random_effects ='SubjectID',plot_heatmap = F, plot_scatter = F)
    res <- fit$results[,c('feature','pval','metadata'), drop =F] %>% filter(metadata =='time') %>% column_to_rownames('feature') %>% dplyr::select(-metadata)
    res$time.padj <- p.adjust(res$pval, method = 'fdr') 
    res <- res %>% dplyr::select(-pval)
  }
  
  if(data.type=='longitudinal'){
    meta$interaction <- as.numeric(meta$X) * meta$time
    fit = Maaslin2(input_data = otu, input_metadata = meta, output = output, min_prevalence = 0,
                   fixed_effects = c("X","time","interaction","Z"),random_effects ='SubjectID',plot_heatmap = F, plot_scatter = F)
    X <- fit$results[,c('feature','pval','metadata'), drop =F] %>% dplyr::filter(metadata =='X') %>% dplyr::select(-metadata)
    X$X.padj <- p.adjust(X$pval, method = 'fdr')
    X <- X %>% dplyr::select(-pval)
    time <- fit$results[,c('feature','pval','metadata'), drop =F] %>% dplyr::filter(metadata =='time') %>% dplyr::select(-metadata) 
    time$time.padj <- p.adjust(time$pval, method = 'fdr') 
    time <- time %>% dplyr::select(-pval)
    interaction <- fit$results[,c('feature','pval','metadata'), drop =F] %>% dplyr::filter(metadata =='interaction') %>% dplyr::select(-metadata)
    interaction$interaction.padj <- p.adjust(interaction$pval, method = 'fdr')
    interaction <- interaction %>% dplyr::select(-pval)
    res <- full_join(X, time) %>% full_join(interaction) %>% column_to_rownames('feature')
  }
  
  ## does not support longitudinal: https://forum.biobakery.org/t/longitudinal-study-and-correcting-for-cage-effect/2063
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  return(list(res=res, na = na, time = t_run, sig = sig))
}




## Add 06/24/2022 
IFAA.func <- function(dat, cutoff = 0.05, data.type ='longitudinal'){
  t_start = Sys.time()
  otu <- as.data.frame(t(dat$counts))
  otu$id <- rownames(otu)
  meta <- dat$meta.dat
  meta$id <- rownames(meta)
  if(data.type=='repeated.measure'){
    fit <- IFAA(MicrobData = otu, CovData = meta, linkIDname = "id", fdrRate = cutoff, 
                adjust_method = "BH", testCov = c("X"), ctrlCov = c("Z"),
                refReadsThresh = 0.01, SDThresh = 0.01, balanceCut = 0.01)
    res <- as.data.frame(fit$full_results$X[,'adj p-value',drop =F])
    colnames(res) <- 'X.padj'
  }
  
  if(data.type=='matched.pair'){
    meta <- meta[,c("time","id"), drop =F]
    fit <- IFAA(MicrobData = otu, CovData = meta, linkIDname = "id", fdrRate = cutoff, 
                adjust_method = "BH", testCov = c("time"), ctrlCov =NULL,
                refReadsThresh = 0.01, SDThresh = 0.01, balanceCut = 0.01)
    res <- as.data.frame(fit$full_results$time[,'adj p-value',drop =F])
    colnames(res) <- 'time.padj'
  }
  
  if(data.type=='longitudinal'){
    meta$interaction <- meta$X * meta$time
    fit <- IFAA(MicrobData = otu, CovData = meta, linkIDname = "id", fdrRate = cutoff, 
                adjust_method = "BH", testCov = c("X","time","interaction"), ctrlCov = c("Z"),
                refReadsThresh = 0.01, SDThresh = 0.01, balanceCut = 0.01)
    res <- as.data.frame(fit$full_results$X[,'adj p-value',drop =F]) %>% dplyr::rename(X.padj = 'adj p-value') %>% rownames_to_column('otu') %>%
      full_join(as.data.frame(fit$full_results$time[,'adj p-value',drop =F]) %>% dplyr::rename(time.padj = 'adj p-value') %>% rownames_to_column('otu') ) %>%
      full_join(as.data.frame(fit$full_results$interaction[,'adj p-value',drop =F])%>% dplyr::rename(interaction.padj = 'adj p-value') %>% rownames_to_column('otu') ) %>% column_to_rownames('otu')
  }
  
  t_run = as.numeric(difftime(Sys.time(), t_start,units = 'secs'))
  na <- apply(res, 2, function(x)sum(is.na(x))/length(x))
  sig <- apply(res, 2, function(x) sum(x[!is.na(x)] <= cutoff)) 
  return(list(res=res, na = na, time = t_run, sig = sig))
}






