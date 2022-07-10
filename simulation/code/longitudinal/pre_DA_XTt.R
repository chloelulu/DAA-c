func <- function(part, method, paras) {
  
  set.seed(part)
  methods <- paras$methods
  model.paras <- paras$model.paras
  
  nSubjects <- paras$nSubjects
  nTimes <- paras$nTimes
  balanced.Xs <- paras$balanced.Xs
  balanced.Ts <- paras$balanced.Ts
  balanced.XTs <- paras$balanced.XTs
  matched.pair <- paras$matched.pair
  error.sds <- paras$error.sds
  
  MgXs <- paras$MgXs
  SgXs <- paras$SgXs
  X.diff.otu.pcts <- paras$X.diff.otu.pcts
  
  MgTs <- paras$MgTs
  SgTs <- paras$SgTs
  SbTs <- paras$SbTs
  time.diff.otu.pcts <- paras$time.diff.otu.pcts
  
  MgXTs <- paras$MgXTs
  SgXTs <- paras$SgXTs
  interaction.diff.otu.pcts <- paras$interaction.diff.otu.pcts
  
  measures <- c('TPR', 'FDR', 'FP','FPR', 'na','time','MCC','F1')
  
  resdir <- paras$resdir
  prefix <- paras$prefix
  resdir <- gsub('/$', '', resdir)
  
  # 23 methods totally
  methods_funs <- list('glmernb'='glmernb.func','lme4'='lme4.func','GLMMPQL'='GLMMPQL.func','ZIBR'='ZIBR.func','LDM'='LDM.func','NBMM'='NBMM.func','MaAsLin2'='MaAsLin2.func',
                       'ZIGMM'='ZIGMM.func','ZINBMM'='ZINBMM.func','metamicrobiomeR'='metamicrobiomeR.func','glmmadaptive'='glmmadaptive.func',
                       'glmmTMB'='glmmTMB.func','LinDA'='LinDA.func','lme' = 'lme.func','lmer'='lmer.func',
                       'MaAsLin2CSS'='MaAsLin2CSS.func','MaAsLin2TMM'='MaAsLin2TMM.func','MaAsLin2GMPR'='MaAsLin2GMPR.func','IFAA'='IFAA.func')
  res <- array(NA, c(length(balanced.Xs),length(balanced.Ts),length(balanced.XTs), length(nSubjects), length(nTimes), length(error.sds), length(MgXs), length(SgXs), length(X.diff.otu.pcts), length(MgTs),
                     length(SgTs), length(SbTs), length(time.diff.otu.pcts), length(MgXTs), length(SgXTs),length(interaction.diff.otu.pcts),
                     length(measures), length(method)),
               dimnames=list(balanced.Xs,balanced.Ts,balanced.XTs,nSubjects,nTimes,error.sds,MgXs,SgXs,X.diff.otu.pcts,MgTs,SgTs,SbTs,time.diff.otu.pcts,MgXTs,SgXTs,interaction.diff.otu.pcts,
                             measures, method))
  
  res <- array(NA, c(length(balanced.Xs),length(balanced.Ts),length(balanced.XTs), length(nSubjects), length(nTimes), length(error.sds), length(MgXs), length(SgXs), length(X.diff.otu.pcts), length(MgTs),
                     length(SgTs), length(SbTs), length(time.diff.otu.pcts), length(MgXTs), length(SgXTs),length(interaction.diff.otu.pcts),
                     length(measures), length(method)),
               dimnames=list(balanced.Xs,balanced.Ts,balanced.XTs, nSubjects,nTimes,error.sds,MgXs,SgXs,X.diff.otu.pcts,MgTs,SgTs,SbTs,time.diff.otu.pcts,MgXTs,SgXTs,interaction.diff.otu.pcts,
                             measures, method))
  
  sink(file.path(resdir, paste(prefix, "_",  part, '_',method,".log", sep="")))
  cat(date(), '\n')
  path0 = paste0(resdir,'/demo/',part)
if(!dir.exists(path0)){dir.create(path0)}

  # nSubject = nSubjects[1];nTime = nTimes[1];matched.pair=F;error.sd = error.sds[1];MgX=MgXs[1];SgX=SgXs[1];X.diff.otu.pct=X.diff.otu.pcts[1];
  # MgT=MgTs[1];SgT=SgTs[1];SbT=SbTs[1];time.diff.otu.pct=time.diff.otu.pcts[1];MgXT=MgXTs[1];SgXT=SgXTs[1];interaction.diff.otu.pct=interaction.diff.otu.pcts[1];
  
  res_seqs <- list()
  for(balanced.X in balanced.Xs){
    for(balanced.T in balanced.Ts){
      for(balanced.XT in balanced.XTs){
        for (nSubject in nSubjects){
          for (nTime in nTimes){
            for (matched.pair in matched.pair){
              for (error.sd in error.sds){
                for (MgX in MgXs){
                  for (SgX in SgXs){
                    for (X.diff.otu.pct in X.diff.otu.pcts){
                      for (MgT in MgTs){
                        for (SgT in SgTs){
                          for (SbT in SbTs){
                            for (time.diff.otu.pct in time.diff.otu.pcts){
                              for (MgXT in MgXTs){
                                for (SgXT in SgXTs){
                                    for (interaction.diff.otu.pct in interaction.diff.otu.pcts){
                                      cat(getwd())
                                      load(paste0('../iter',part,'/','data',balanced.X,'_',balanced.T,'_',balanced.XT,'_',nSubject,'_',nTime,'_',error.sd,'_',MgX,'_',SgX,'_',X.diff.otu.pct,'_',MgT,'_',SgT,'_',SbT,'_',time.diff.otu.pct,'_',MgXT,'_',SgXT,'_',interaction.diff.otu.pct,'.Rdata'))
output = paste0(resdir,'/demo/',part,'/', 'data',balanced.X,'_',balanced.T,'_',balanced.XT,'_',nSubject,'_',nTime,'_',error.sd,'_',MgX,'_',SgX,'_',X.diff.otu.pct,'_',MgT,'_',SgT,'_',SbT,'_',time.diff.otu.pct,'_',MgXT,'_',SgXT,'_',interaction.diff.otu.pct)

if(!dir.exists(output)){dir.create(output)}
                                      
                                      # apply each DA method
                                      MCC <- F1 <- fpr <- fdr <- tpr <- fp <- na <- time <- NA;res_seq <- NULL
                                      wrapper <- match.fun(methods_funs[[method]])
                                      tryCatch({
                                        cat(paste(method, '\n'))
if(method %in% c('MaAsLin2','MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR')){
    wrapper.obj <- wrapper(dat,cutoff=0.05,output = output, data.type ='longitudinal')
  }else{
    wrapper.obj <- wrapper(dat,cutoff=0.05, data.type ='longitudinal')
  }
                                        
fdr <- wrapper.obj$res.obj$time['fdr']
                                        tpr <- wrapper.obj$res.obj$time['tpr']
                                        fp <- wrapper.obj$res.obj$time['fp']
                                        fpr <- wrapper.obj$res.obj$time['fpr']
                                        MCC <- wrapper.obj$res.obj$time['MCC']
                                        F1 <- wrapper.obj$res.obj$time['F1']
                                        res_seq <- as.data.frame(wrapper.obj$res)
                                        na <- wrapper.obj$na
                                        time <- wrapper.obj$time
                                        res[as.character(balanced.X),as.character(balanced.T),as.character(balanced.XT),as.character(nSubject),as.character(nTime),as.character(error.sd),as.character(MgX),as.character(SgX),as.character(X.diff.otu.pct),as.character(MgT),
                                            as.character(SgT),as.character(SbT),as.character(time.diff.otu.pct),as.character(MgXT),as.character(SgXT),
                                            as.character(interaction.diff.otu.pct),'TPR', method] <- tpr
                                        res[as.character(balanced.X),as.character(balanced.T),as.character(balanced.XT),as.character(nSubject),as.character(nTime),as.character(error.sd),as.character(MgX),as.character(SgX),as.character(X.diff.otu.pct),as.character(MgT),
                                            as.character(SgT),as.character(SbT),as.character(time.diff.otu.pct),as.character(MgXT),as.character(SgXT),
                                            as.character(interaction.diff.otu.pct), 'FDR', method] <-  fdr
                                        
                                        res[as.character(balanced.X),as.character(balanced.T),as.character(balanced.XT),as.character(nSubject),as.character(nTime),as.character(error.sd),as.character(MgX),as.character(SgX),as.character(X.diff.otu.pct),as.character(MgT),
                                            as.character(SgT),as.character(SbT),as.character(time.diff.otu.pct),as.character(MgXT),as.character(SgXT),
                                            as.character(interaction.diff.otu.pct), 'FP', method] <-  fp
                                        
                                        res[as.character(balanced.X),as.character(balanced.T),as.character(balanced.XT),as.character(nSubject),as.character(nTime),as.character(error.sd),as.character(MgX),as.character(SgX),as.character(X.diff.otu.pct),as.character(MgT),
                                            as.character(SgT),as.character(SbT),as.character(time.diff.otu.pct),as.character(MgXT),as.character(SgXT),
                                            as.character(interaction.diff.otu.pct), 'FPR', method] <-  fpr
                                        
                                        res[as.character(balanced.X),as.character(balanced.T),as.character(balanced.XT),as.character(nSubject),as.character(nTime),as.character(error.sd),as.character(MgX),as.character(SgX),as.character(X.diff.otu.pct),as.character(MgT),
                                            as.character(SgT),as.character(SbT),as.character(time.diff.otu.pct),as.character(MgXT),as.character(SgXT),
                                            as.character(interaction.diff.otu.pct), 'na', method] <-  na
                                        
                                        res[as.character(balanced.X),as.character(balanced.T),as.character(balanced.XT),as.character(nSubject),as.character(nTime),as.character(error.sd),as.character(MgX),as.character(SgX),as.character(X.diff.otu.pct),as.character(MgT),
                                            as.character(SgT),as.character(SbT),as.character(time.diff.otu.pct),as.character(MgXT),as.character(SgXT),
                                            as.character(interaction.diff.otu.pct), 'time', method] <-  time
                                        
                                        res[as.character(balanced.X),as.character(balanced.T),as.character(balanced.XT),as.character(nSubject),as.character(nTime),as.character(error.sd),as.character(MgX),as.character(SgX),as.character(X.diff.otu.pct),as.character(MgT),
                                            as.character(SgT),as.character(SbT),as.character(time.diff.otu.pct),as.character(MgXT),as.character(SgXT),
                                            as.character(interaction.diff.otu.pct), 'MCC', method] <-  MCC
                                        
                                        res[as.character(balanced.X),as.character(balanced.T),as.character(balanced.XT),as.character(nSubject),as.character(nTime),as.character(error.sd),as.character(MgX),as.character(SgX),as.character(X.diff.otu.pct),as.character(MgT),
                                            as.character(SgT),as.character(SbT),as.character(time.diff.otu.pct),as.character(MgXT),as.character(SgXT),
                                            as.character(interaction.diff.otu.pct), 'F1', method] <-  F1
                                        
                                        res_seqs[[paste0(balanced.X,balanced.T, balanced.XT,nSubject,nTime,error.sd,MgX,SgX,X.diff.otu.pct,MgT,SgT,SbT,time.diff.otu.pct,MgXT,SgXT,interaction.diff.otu.pct, method)]] <- res_seq
                                      }, error =function(e){cat(paste0(paste0(balanced.Xs,balanced.Ts,balanced.XTs,nSubject,nTime,error.sd,MgX,SgX,X.diff.otu.pct,MgT,SgT,SbT,time.diff.otu.pct,MgXT,SgXT,interaction.diff.otu.pct), method, ' ERROR : '),conditionMessage(e), "\n")})
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
    }
    warnings()
    cat('\n', date(), '\n', method,'Finished!')
    sink()
    save(res_seqs,file=file.path(resdir, paste(prefix, "_summarymatrix",  part, '-', method, ".Rdata", sep="")))
    #save(res, file=file.path(resdir, paste(prefix, "_res",  part, '-', method, ".Rdata", sep="")))
    return(res)
  }
}



