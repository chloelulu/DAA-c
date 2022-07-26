rm(list = ls())
library(gtable)
library(gridExtra)
library(circlize)

source('/research/bsi/projects/staff_analysis/m216453/correlated/code/DM1.R')
tm = load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35.RData')
tm1 = load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35_dirmult.RData')

# Count generation modeL
nSubject =2; nOTU = 10; nTime = 3;matched.pair=F; balanced.X = T;balanced.T = T;balanced.XT = T;
# random error
error.sd = 4; error.mean = 0;
# Covariate effect
MgX = 1; SgX = 0; X.diff.otu.pct = 0.1;
# Interaction effect
MgXT = 0; SgXT = 0; interaction.diff.otu.pct = 0;
# Time effect
MgT = 1; SgT = 0; SbT = 1; time.diff.otu.pct = 0.1;
# Confounding effect
grp.ratio = 1; conf.cov.cor = 0.6; Z.diff.otu.pct = 0.05; MgZ = 1; SgZ = 0;Z.nondiff.otu.pct = 0.1;
# Sequence depth
depth.mu = 10000; depth.theta = 5;  depth.conf.factor = 0



#### ---- Step 1: The observed counts from a reference dataset Cki (1≤k≤M, 1≤i≤N) --- ####
## show 1st 2 and last 2 samples in figures
idx.sam1 <- colnames(otu.tab)[1:2]
idx.sam2 <- colnames(otu.tab)[294:295]

sub.count1 <- otu.tab[,idx.sam1]
sub.count1 <- sub.count1[order(-rowMeans(sub.count1)),]
col_fun = colorRamp2(c(min(sub.count1), mean(sub.count1),max(sub.count1)), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
colnames(sub.count1) = c('S1','S2')
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/observedcount1.pdf", width = 1, height = 10)
ComplexHeatmap::Heatmap(sub.count1, col = col_fun,
                        row_order = rownames(sub.count1),
                        column_order = colnames(sub.count1),
                        show_row_dend = F, show_column_dend = F, name = " ",
                        show_heatmap_legend = T, show_row_names = F, 
                        show_column_names =T)
dev.off()


sub.count2 <- otu.tab[,idx.sam2]
colnames(sub.count2) = c('S294','S295')
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/observedcount2.pdf", width = 1, height = 10)
ComplexHeatmap::Heatmap(sub.count2, col = col_fun,show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(sub.count1),
                        column_order = colnames(sub.count2),
                        name = " ",show_heatmap_legend = T, show_row_names = F, show_column_names = T)
dev.off()


#### ----- Step 2: Estimate the underlying true composition (Pki) using empirical Bayes ---- ####
set.seed(1) ### If no seed here, the result will be non-reproducible. rgamma() random generation 
model.paras <- EstPara(otu.tab,dirmult.paras = dirmult.paras)
absolute <- model.paras$otu.tab
max(absolute)
absolute <- absolute[rownames(otu.tab),colnames(otu.tab)]

sub.absolute1 <- absolute[rownames(sub.count1),idx.sam1]
col_fun = colorRamp2(c(min(sub.absolute1), mean(sub.absolute1),max(sub.absolute1)), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/absolute_abundance1.pdf", width = 1, height = 10)
ComplexHeatmap::Heatmap(sub.absolute1, col = col_fun, 
                        row_order = rownames(sub.absolute1),
                        column_order = colnames(sub.absolute1),
                        show_row_dend = F, show_column_dend = F, show_heatmap_legend = T, 
                        name = " ",show_row_names = F, show_column_names = F)
dev.off()

sub.absolute2 <- absolute[rownames(sub.count2),idx.sam2]
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/absolute_abundance2.pdf", width = 1, height = 10)
ComplexHeatmap::Heatmap(sub.absolute2, 
                        row_order = rownames(sub.absolute2),
                        column_order = colnames(sub.absolute2),
                        col = col_fun, show_row_dend = F, show_column_dend = F, show_heatmap_legend = T, 
                        name = " ",show_row_names = F, show_column_names = F)
dev.off()


otu_tab <- (sub.absolute1[order(-rowMeans(sub.absolute1)),])[c(1:10),]
colnames(otu_tab) <- paste0('S',1:ncol(otu_tab))
rownames(otu_tab) <- paste0('taxon',1:nrow(otu_tab))

col_fun = colorRamp2(c(min(otu_tab), mean(otu_tab),max(otu_tab)), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/baseline.pdf", width = 3, height = 3)
ComplexHeatmap::Heatmap(otu_tab, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(otu_tab),
                        column_order = colnames(otu_tab),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


#### ----- Step 3: Randomly sample I subjects with K most abundant taxa and replicate their abundance profiles J times (K=10, I=2, J=3 in this example) ---####
otu.tab <- otu_tab[, rep(1:nSubject, each=nTime)]
colnames(otu.tab) <- paste0( rep(paste0('S',1:nSubject), each = nTime),':',paste0('R',1:nTime))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/baseline_rep.pdf", width = 3, height = 3)
ComplexHeatmap::Heatmap(otu.tab, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(otu.tab),
                        column_order = colnames(otu.tab),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()



#### ------- Step 4: Generate covariate/confounder/time variables and their coefficients -----####
nSam <- ncol(otu.tab)
SubjectID <- as.numeric(as.factor(gsub(':.*','',colnames(otu.tab))))
time <- c(1,2,3,1,2,3)#as.numeric(as.factor(gsub(':.*','',colnames(otu.tab))))-1 
X <- as.vector(ifelse(SubjectID <= quantile(unique(SubjectID), grp.ratio / (1 + grp.ratio)), 0, 1))
rho <- sqrt(conf.cov.cor ^ 2 / (1 - conf.cov.cor ^ 2))
Z <- c(0.6, 0.6, 0.6, -0.6, -0.6, -0.6)#as.vector(rho * scale(X) + rnorm(nSam))
meta <- as.data.frame(cbind(X, Z, SubjectID, time))
rownames(meta) <- colnames(otu.tab)


## X variable
X.m <- t(replicate(nOTU, X))
colnames(X.m) <- colnames(otu.tab)
rownames(X.m) <- (rownames(otu.tab))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/X.pdf")
g <- tableGrob(X.m)
tg <- grid.draw(g)
dev.off()
## time variable
time.m <- t(replicate(nOTU, time))
colnames(time.m) <- colnames(otu.tab)
rownames(time.m) <- (rownames(otu.tab))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/time.pdf")
g <- tableGrob(time.m)
tg <- grid.draw(g)
dev.off()

## Z variable
Z.m <- t(replicate(nOTU, Z))
colnames(Z.m) <- colnames(otu.tab)
rownames(Z.m) <- (rownames(otu.tab))
Z.m <- round(Z.m, 2)
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/Z.pdf")
g <- tableGrob(Z.m)
tg <- grid.draw(g)
dev.off()


## Generate diff.otu
otu.ord <- 1:nOTU
X.diff.otu.ind <- time.diff.otu.ind <- interaction.diff.otu.ind <- NULL

X.diff.otu.num <- round(X.diff.otu.pct * nOTU)
X.diff.otu.ind <- c(X.diff.otu.ind, sample(otu.ord, X.diff.otu.num))

time.diff.otu.num <- round(time.diff.otu.pct * nOTU)
time.diff.otu.ind <- c(time.diff.otu.ind, sample(otu.ord, time.diff.otu.num))

Z.diff.otu.num <- round(Z.diff.otu.pct * nOTU)
Z.nondiff.otu.num <- round(Z.nondiff.otu.pct * nOTU)
Z.diff.otu.ind <- c(sample(X.diff.otu.ind, Z.diff.otu.num), sample(setdiff(otu.ord, X.diff.otu.ind), Z.nondiff.otu.num))


## X coefficient
coef.X <- sample(c(rnorm(floor(nOTU / 2), mean = -MgX, sd = SgX), rnorm(nOTU - floor(nOTU / 2), mean = MgX, sd = SgX)))
coef.X[setdiff(otu.ord, X.diff.otu.ind)] <- 0 

coef.X.m <- replicate(nSam, coef.X)
rownames(coef.X.m) <- rownames(otu.tab)
colnames(coef.X.m) <- rev(colnames(otu.tab))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/coef.X.pdf")
g <- tableGrob(coef.X.m)
tg <- grid.draw(g)
dev.off()


## aX
eta.X <- coef.X  %*% t(scale(X)) 
eta.X <- round(eta.X)
rownames(eta.X) <- rownames(otu.tab)
colnames(eta.X) <- colnames(otu.tab)
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/eta.X.pdf")
g <- tableGrob(eta.X)
tg <- grid.draw(g)
dev.off()


col_fun = colorRamp2(c(min(eta.X), 0,max(eta.X)), c(brewer.pal(9,'Blues')[7],'white', brewer.pal(9,'Reds')[8]))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/eta.X.pdf", width = 3.5, height = 3)
ComplexHeatmap::Heatmap(eta.X, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(eta.X),
                        column_order = colnames(eta.X),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


## time coefficient
set.seed(123)
coef.T <- sample(c(rnorm(floor(nOTU / 2), mean = -MgT, sd = SgT), rnorm(nOTU - floor(nOTU / 2), mean = MgT, sd = SgT)))
coef.T[setdiff(otu.ord, time.diff.otu.ind)] <- 0 
vc <- c(0.5,0.5,0.5,0.8,0.8,0.8)
vc0 <- rep(0,6)
coef.T <- t(replicate(nOTU,vc0)) # for illustration purpose
coef.T[2,] <- vc
rownames(coef.T) <- rownames(otu.tab)
colnames(coef.T) <- colnames(otu.tab)
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/coef.T.pdf")
g <- tableGrob(round(coef.T,2))
tg <- grid.draw(g)
dev.off()

## c*T
eta.time <- t(t(coef.T) * as.vector(time))
col_fun = colorRamp2(c(min(eta.time), median(eta.time),max(eta.time)), c('white',brewer.pal(9,'Blues')[5], brewer.pal(9,'Blues')[8]))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/eta.time.pdf", width = 3.5, height = 3)
ComplexHeatmap::Heatmap(eta.time, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(eta.time),
                        column_order = colnames(eta.time),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


## Generate coefficient for Z, assume balanced change
set.seed(123)
coef.Z <- sample(c(rnorm(floor(nOTU / 2), mean = -MgZ, sd = SgZ), rnorm(nOTU - floor(nOTU / 2), mean = MgZ, sd = SgZ)))
coef.Z[setdiff(otu.ord, Z.diff.otu.ind)] <- 0
coef.Z.m <- replicate(nSam, coef.Z)
rownames(coef.Z.m) <- rownames(otu.tab)
colnames(coef.Z.m) <- colnames(otu.tab)
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/coef.Z.pdf")
g <- tableGrob(round(coef.Z.m,2))
tg <- grid.draw(g)
dev.off()


## b*Z
eta.Z <- coef.Z %*% t(scale(Z)) 
rownames(eta.Z) <- rownames(otu.tab)
colnames(eta.Z) <- colnames(otu.tab)
col_fun = colorRamp2(c(min(eta.Z), 0,max(eta.Z)), c(brewer.pal(9,'Blues')[7],'white', brewer.pal(9,'Reds')[9]))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/eta.Z.pdf", width = 3.5, height = 3)
ComplexHeatmap::Heatmap(eta.Z, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(eta.Z),
                        column_order = colnames(eta.Z),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()

## Generate random error
error <- replicate(nSam,rnorm(nOTU, error.mean, error.sd))
rownames(error) <- rownames(otu.tab)
colnames(error) <- colnames(otu.tab)
col_fun = colorRamp2(c(min(error), 0,max(error)), c(brewer.pal(9,'Blues')[7],'white', brewer.pal(9,'Reds')[8]))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/error.pdf", width = 3.5, height = 3)
ComplexHeatmap::Heatmap(error, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(error),
                        column_order = colnames(error),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


#### ---- Step 5: Generate log fold change due to the covariate/confounder/time effects ---- ####
eta.exp <- (eta.X + eta.Z + eta.time + error)
col_fun = colorRamp2(c(min(eta.exp), 0,max(eta.exp)), c(brewer.pal(9,'Blues')[7],'white', brewer.pal(9,'Reds')[8]))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/eta.pdf", width = 3.5, height = 3)
ComplexHeatmap::Heatmap(eta.exp, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(eta.exp),
                        column_order = colnames(eta.exp),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


#### ---- Step 6: Obtain the composition after adding covariate/confounder/time effects and random errors ---- ####
eta.exp <- exp(t(eta.exp)) * t(otu.tab)
otu.tab.prop <- eta.exp / rowSums(eta.exp)
otu.tab.prop <- t(otu.tab.prop) 
col_fun = colorRamp2(c(min(otu.tab.prop), median(otu.tab.prop),max(otu.tab.prop)), c(brewer.pal(9,'Blues')[1],brewer.pal(9,'Reds')[3], brewer.pal(9,'Reds')[8]))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/eta_norm.pdf", width = 3.5, height = 3)
ComplexHeatmap::Heatmap(otu.tab.prop, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(otu.tab.prop),
                        column_order = colnames(otu.tab.prop),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


#### ---- Step 7. Given the sequencing depth Dij, generate the read counts C’kij for sample i based on a multinomial distribution with parameters (Dij, P’1ij,..., P’kij) ----####
nSeq <- rnegbin(ncol(otu.tab.prop), mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta)
otu.tab.sim <- sapply(1:ncol(otu.tab.prop), function (i) rmultinom(1, nSeq[i], otu.tab.prop[,i]))
colnames(otu.tab.sim) <- rownames(eta.exp)
rownames(otu.tab.sim) <- rownames(otu.tab)
col_fun = colorRamp2(c(min(otu.tab.sim), median(otu.tab.sim),max(otu.tab.sim)), c(brewer.pal(9,'Greens')[1],brewer.pal(9,'Greens')[5], brewer.pal(9,'Greens')[9]))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/otu.tab.sim.pdf", width = 3.5, height = 3)
ComplexHeatmap::Heatmap(otu.tab.sim, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(otu.tab.sim),
                        column_order = colnames(otu.tab.sim),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


