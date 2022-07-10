require(MASS); require(rmutil)
rdirichlet.m <- function (alpha) {
  Gam <- matrix(rgamma(length(alpha), shape = alpha), nrow(alpha), ncol(alpha))
  t(t(Gam) / colSums(Gam))
}
EstPara0 <- function (otu.tab) {
  
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste0('OTU', 1 : nrow(otu.tab))
  } # otu * sample
  samplenames = colnames(otu.tab)
  taxnames = rownames(otu.tab)
  
  dirmult.paras <- dirmult::dirmult(t(otu.tab))
  
  gamma = dirmult.paras$gamma
  names(gamma) = names(dirmult.paras$pi)
  
  # Add pseduo count(each OTU add gamma estimated from dirmult)
  otu.tab = sapply(1:ncol(otu.tab), function (i) gamma + otu.tab[,i]) # C_ij otu * sample
  
  # back to dirchlet, calculate the true proportion
  otu.tab.p <- rdirichlet.m(otu.tab) # P_ij nOTU*nSam
  colnames(otu.tab.p) = samplenames
  rownames(otu.tab.p) = taxnames
  
  # order OTUs by mean OTU proportion, for later selection
  ord = order(rowMeans(otu.tab.p), decreasing = TRUE)
  otu.tab.p =  otu.tab.p[ord,]
  
  # apply size factor
  Si = exp(rnorm(ncol(otu.tab.p)))
  otu.tab0 = t(t(otu.tab.p)*Si)
  colnames(otu.tab0) = colnames(otu.tab.p)
  return(list(mu = otu.tab.p, otu.tab = otu.tab0))
}


EstPara <- function (otu.tab, dirmult.paras) {
  
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste0('OTU', 1 : nrow(otu.tab))
  } # otu * sample
  
  samplenames = colnames(otu.tab)
  taxnames = rownames(otu.tab)
  
  # For saving time, here will use the parameters already estimated
  gamma = dirmult.paras$gamma
  names(gamma) = names(dirmult.paras$pi)
  gamma = gamma[taxnames]
  
  # Add pseduo count(each OTU add gamma estimated from dirmult)
  otu.tab = sapply(1:ncol(otu.tab), function (i) gamma + otu.tab[,i]) # C_ij otu * sample
  
  # back to dirchlet, calculate the true proportion
  otu.tab.p <- rdirichlet.m(otu.tab) # P_ij nOTU*nSam
  colnames(otu.tab.p) = samplenames
  rownames(otu.tab.p) = taxnames
  
  # order OTUs by mean OTU proportion, for later selection
  ord = order(rowMeans(otu.tab.p), decreasing = TRUE)
  otu.tab.p =  otu.tab.p[ord,]
  
  # apply size factor
  Si = exp(rnorm(ncol(otu.tab.p)))
  otu.tab0 = t(t(otu.tab.p)*Si)
  colnames(otu.tab0) = colnames(otu.tab.p)
  
  return(list(mu = otu.tab.p, otu.tab = otu.tab0))
}

SimulateSeq <- function (
  # Input data
  otu.tab, # otu * sample
  # design matrix 
  # [change to nTime for consistency]
  nSubject = 100, nOTU = 500, nTime = 5, matched.pair=F,
  # random error
  error.sd = 1, error.mean = 0, 
  # Covariate effect alpha = 0.07 is enough to have around 10 fold change
  MgX = 0.5, SgX = 0, X.diff.otu.pct = 0, grp.ratio = 1, alpha = 0, balanced.X = T, 
  # Time effect
  MgT = 1, SgT = 0, SbT = 0, time.diff.otu.pct =0, balanced.T = T,
  # Interaction effect
  MgXT = 0, SgXT = 0, interaction.diff.otu.pct = 0, balanced.XT = T,
  # Confounding effect
  MgZ = 0.5, SgZ = 0,  Z.diff.otu.pct = 0.05, conf.cov.cor = 0.6, Z.nondiff.otu.pct = 0.1,
  # Sequence depth
  depth.mu = 10000, depth.theta = 5,  depth.conf.factor = 0
){
  

  ## Estimated parameters
  model.paras <- EstPara(otu.tab, dirmult.paras = dirmult.paras)
  
  #### select nSubject first then top nOTU
  sample.names <- colnames(model.paras$otu.tab)
  idx.sample <- sample(sample.names, nSubject)
  OTU <- model.paras$otu.tab[,idx.sample, drop=F]
  otu.names <- names(sort(rowMeans(OTU),decreasing=T))[1:nOTU]
  otu_tab <- OTU[otu.names,, drop=F]

  
  ## Replicate template sample in each Subject, generate baseline otu table
  otu.tab <- otu_tab[, rep(1:nSubject, each=nTime)]
  colnames(otu.tab) <- paste0(paste0('rep',1:nTime), '_',rep(paste0('Subject',1:nSubject), each = nTime))
 
  prop.max <- otu.tab[which.max(rowMeans(otu.tab)),]
  adj.para <- t(apply(otu.tab,1,function(x) (prop.max / x)^alpha))
  max(adj.para);min(adj.para);mean(adj.para);median(adj.para)
  
  ## Generate meta data 
  nSam <- ncol(otu.tab)
  SubjectID <- as.numeric(gsub('.*Subject','',colnames(otu.tab)))
  time <- as.numeric(as.factor(gsub('_.*','',colnames(otu.tab))))-1 
  if(matched.pair){
    X <- time
  }else{
    X <- as.vector(ifelse(SubjectID <= quantile(unique(SubjectID), grp.ratio / (1 + grp.ratio)), 0, 1))
  }
  rho <- sqrt(conf.cov.cor ^ 2 / (1 - conf.cov.cor ^ 2))
  Z <- as.vector(rho * scale(X) + rnorm(nSam))
  meta <- as.data.frame(cbind(X, Z, SubjectID, time))
  rownames(meta) <- colnames(otu.tab)
  
  
  ## Generate random error
  error <- replicate(nSam,rnorm(nOTU, error.mean, error.sd))
  
  ## Generate diff.otu
  otu.ord <- 1:nOTU
  X.diff.otu.ind <- time.diff.otu.ind <- interaction.diff.otu.ind <- NULL
  
  X.diff.otu.num <- round(X.diff.otu.pct * nOTU)
  X.diff.otu.ind <- c(X.diff.otu.ind, sample(otu.ord, X.diff.otu.num))
  
  time.diff.otu.num <- round(time.diff.otu.pct * nOTU)
  time.diff.otu.ind <- c(time.diff.otu.ind, sample(otu.ord, time.diff.otu.num))
  
  Z.diff.otu.num <- round(Z.diff.otu.pct * nOTU)
  Z.nondiff.otu.num <- round(Z.nondiff.otu.pct * nOTU)
  if(length(X.diff.otu.ind)==0){
    # for matched case
    Z.diff.otu.ind <- c(sample(time.diff.otu.ind, Z.diff.otu.num), sample(setdiff(otu.ord, X.diff.otu.ind), Z.nondiff.otu.num))
  }else{
    Z.diff.otu.ind <- c(sample(X.diff.otu.ind, Z.diff.otu.num), sample(setdiff(otu.ord, X.diff.otu.ind), Z.nondiff.otu.num))
  }
  
  
  interaction.diff.otu.num <- round(interaction.diff.otu.pct * nOTU)
  interaction.diff.otu.ind <- c(interaction.diff.otu.ind, sample(otu.ord, interaction.diff.otu.num))

  ## Generate coefficient for X
  if(balanced.X){
    coef.X <- sample(c(rnorm(floor(nOTU / 2), mean = -MgX, sd = SgX), rnorm(nOTU - floor(nOTU / 2), mean = MgX, sd = SgX)))
  }else{
    coef.X <- rnorm(nOTU, MgX, SgX)
  }
  
  coef.X[setdiff(otu.ord, X.diff.otu.ind)] <- 0
  eta.X <- coef.X  %*% t(scale(X)) * adj.para

  ## Generate coefficient for time
  # generate population mean
  if(balanced.T){
    coef.T <- sample(c(rnorm(floor(nOTU / 2), mean = -MgT, sd = SgT), rnorm(nOTU - floor(nOTU / 2), mean = MgT, sd = SgT)))
  }else{
    coef.T <- rnorm(nOTU, MgT, SgT)
  }
  coef.T[setdiff(otu.ord, time.diff.otu.ind)] <- 0 
  coef.T.Subject <- replicate(nSubject, rnorm(nOTU, coef.T, SbT)) 
  coef.T <- matrix(data = apply(coef.T.Subject, 2, function(x) rep(x, nTime)), ncol = ncol(coef.T.Subject)*nTime)
  eta.time <- t(t(coef.T) * as.vector(scale(time))) * adj.para


  ## Generate coefficient for Z, assume balanced change
  coef.Z <- sample(c(rnorm(floor(nOTU / 2), mean = -MgZ, sd = SgZ), rnorm(nOTU - floor(nOTU / 2), mean = MgZ, sd = SgZ)))
  coef.Z[setdiff(otu.ord, Z.diff.otu.ind)] <- 0 
  eta.Z <- coef.Z %*% t(scale(Z)) * adj.para
  
  ## Generate coefficient for interaction 
  if(balanced.XT){
    coef.interaction <- sample(c(rnorm(floor(nOTU / 2), mean = -MgXT, sd = SgXT), rnorm(nOTU - floor(nOTU / 2), mean = MgXT, sd = SgXT)))
  }else{
    coef.interaction  <- rnorm(nOTU, MgXT, SgXT)
  }
  coef.interaction[setdiff(otu.ord, interaction.diff.otu.ind)] <- 0 
  coef.interaction.Subject <- replicate(nSubject, rnorm(nOTU, coef.interaction, 0)) 
  coef.interaction <- matrix(data = apply(coef.interaction.Subject, 2, function(x) rep(x, nTime)), ncol = ncol(coef.interaction.Subject)*nTime)
  eta.interaction <-  t(t(coef.interaction) * as.vector(scale(time) *scale(X))) * adj.para 

  ## combine index matrix, coefficient matrix and design matrix
  eta.exp <- exp(t(eta.X + eta.Z + eta.time + eta.interaction + error))
  eta.exp <- eta.exp * t(otu.tab)
  not.trans.eta.exp <- eta.exp
  
  # Renormalize
  otu.tab.prop <- eta.exp / rowSums(eta.exp)
  otu.tab.prop <- t(otu.tab.prop) 
  
  # Generate the sequence depth
  nSeq <- rnegbin(ncol(otu.tab.prop), mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta)
  otu.tab.sim <- sapply(1:ncol(otu.tab.prop), function (i) rmultinom(1, nSeq[i], otu.tab.prop[,i]))
  
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- rownames(otu.tab)
  
  X.diff.otu.ind = (1 : nOTU) %in% X.diff.otu.ind
  time.diff.otu.ind = (1 : nOTU) %in% time.diff.otu.ind
  interaction.diff.otu.ind = (1 : nOTU) %in% interaction.diff.otu.ind
  Z.diff.otu.ind = (1 : nOTU) %in% Z.diff.otu.ind
  
  return(list(otu.tab.sim = otu.tab.sim, meta = meta, otu.names = otu.names, error = error,
              not.trans.eta.exp = not.trans.eta.exp, 
              otu.tab = otu.tab, # baseline otu table
              eta.X = eta.X, eta.time = eta.time, eta.interaction= eta.interaction, eta.Z = eta.Z,
              X.diff.otu.ind = X.diff.otu.ind, time.diff.otu.ind = time.diff.otu.ind, 
              Z.diff.otu.ind = Z.diff.otu.ind, interaction.diff.otu.ind = interaction.diff.otu.ind))
}

