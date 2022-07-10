func <- function(part, paras) {
  set.seed(part)
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
  
  resdir <- paras$resdir
  prefix <- paras$prefix
  resdir <- gsub('/$', '', resdir)
  
  sum.list <- list()
  cat(date(), '\n')
  dir.create(paste0('../iter',part))
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
                                    Sim.obj <- SimulateSeq(otu.tab = otu.tab, nOTU = 500, nSubject = nSubject, nTime =nTime,
                                                           error.sd =error.sd, error.mean = 0, # 2,0.2
                                                           matched.pair = matched.pair, 
                                                           balanced.X =balanced.X, balanced.T =balanced.T,balanced.XT =balanced.XT,
                                                           # Covariate effect
                                                           MgX = MgX, SgX = SgX, X.diff.otu.pct = X.diff.otu.pct,
                                                           # Time effect
                                                           MgT = MgT, SgT = SgT, SbT = SbT, time.diff.otu.pct = time.diff.otu.pct,
                                                           # Interaction effect
                                                           MgXT = MgXT, SgXT = SgXT, interaction.diff.otu.pct = interaction.diff.otu.pct)
                                    
                                    Sim.obj$meta$SubjectID <- as.factor(Sim.obj$meta$SubjectID)
                                    
                                    otu.tab.sim <- Sim.obj$otu.tab.sim;dim(otu.tab.sim)
                                    truth <- as.data.frame(cbind(Sim.obj$X.diff.otu.ind,Sim.obj$time.diff.otu.ind,Sim.obj$interaction.diff.otu.ind))
                                    rownames(truth) <- Sim.obj$otu.names
                                    colnames(truth) <- c('X.truth','time.truth','interaction.truth')
                                    meta.dat <- as.data.frame(Sim.obj$meta)
                                    ##-- Preprocessing
                                    gmpr.size.factor <- GMPR::GMPR(t(otu.tab.sim))
                                    gmpr.size.factor[is.na(gmpr.size.factor)] <- 1
                                    ## raw otu table prevelance filtration
                                    prop <- t(t(otu.tab.sim) / colSums(otu.tab.sim));dim(prop)
                                    prop <- prop[rowSums(prop!=0) > 0.1 * ncol(prop), , drop=FALSE];dim(prop)
                                    otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE];dim(prop)
                                    prop <- prop[rowMaxs(prop) > 0.002, , drop=FALSE]
                                    otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE] # used for analysis
                                    
                                    gmpr.size.factor <- gmpr.size.factor[colnames(otu.tab.sim)]
			            truth <- truth[rownames(otu.tab.sim),,drop =F]                                    
                                    dat <- list(counts = otu.tab.sim, prop = prop, truth = truth, meta.dat = meta.dat, gmpr.size.factor = gmpr.size.factor)
                                    save(dat,file=paste0('../iter',part,'/','data',balanced.X,'_',balanced.T,'_',balanced.XT,'_',nSubject,'_',nTime,'_',error.sd,'_',MgX,'_',SgX,'_',X.diff.otu.pct,'_',MgT,'_',SgT,'_',SbT,'_',time.diff.otu.pct,'_',MgXT,'_',SgXT,'_',interaction.diff.otu.pct,'.Rdata'))
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
  }
}




clsapply <- function(dat, func, df=NULL, tempdir="~/project/tempZ",
                     queque="7-day", timespan="540", mem="8G", req.pack='stats',
                     tidy=F) {
  # Apply to the row of dat by default (dim=2 for column)
  # Create temp dir
  if (!file.exists(tempdir)) {
    dir.create(tempdir)
  }
  setwd(tempdir)
  
  save(dat, func, df, otu.tab, dirmult.paras, file=("dat.RData"))
  # methods <- df$methods
  nJob <- length(dat)
  # Create R script
  txt <- paste('
               args=(commandArgs(TRUE))
               if(length(args)==0){
               print("No arguments supplied.")
               ##supply default values
               } else{
               part = args[1]
               part = as.integer(part)
                              }',
               paste(paste0('\nrequire(',req.pack, ')'), collapse="\n"),
               '\n date()
      load("dat.RData") 
      
                        if (is.list(dat)) {
                        item <- dat[[part]]
                        } 
                        if (is.vector(dat)) { 
                        item <- dat[part]
                        } 

                        source("/research/bsi/projects/staff_analysis/m216453/correlated/code/DM1.R")
                        pkg <- c("matrixStats","nlme","lme4","LinDA","ZIBR","miLineage","LDM","aod","compositions","tidyverse","GMPR","glmmTMB","NBZIMM","GLMMadaptive","metamicrobiomeR","tibble")
                        lapply(pkg, require, character.only = TRUE)

                        res0 <- func(item, df)
                        date()
                        ')
  writeLines(txt, "dat_script.R")
  # Submit job
  cat("Submit jobs ...")
  prefix <- paste(c("J", sample(LETTERS, 4, repl=TRUE)), collapse="")
  for (part in 1:nJob) {
    rfile <- "dat_script.R"
    rout <- paste0(part, ".Rout")
    sh <- paste(
      "qsub",
      paste0("-N ", prefix, part),
      "-j y",
      "-cwd",
      "-q", queque,
      "-m abe",
      paste0("-l h_vmem=", mem),
      "-V",
      paste("-o ",  part, ".out", sep=""),
      paste("-e ",  part, ".err", sep=""),
      "-b y",
      paste("\'R CMD BATCH --no-save --no-restore \"--args ", part, "\" ",
            rfile, " ", rout, "\'", sep="")
    )
    print(sh)
    system(sh)
  }
}

