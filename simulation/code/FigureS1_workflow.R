load('/research/bsi/projects/staff_analysis/m216453/correlated/code/Stool_V35.RData')
# Count generation modeL
nSubject =2; nOTU = 10; nTime = 3;matched.pair=F; alpha = 0; balanced.X = T;balanced.T = T;balanced.XT = T;
# random error
error.sd = 4; error.mean = 0;
# Covariate effect
MgX = 1; SgX = 0; X.diff.otu.pct = 0.1;
# Interaction effect
MgXT = 1; SgXT = 0; interaction.diff.otu.pct = 0.1;
# Time effect
MgT = 3; SgT = 0; SbT = 1; time.diff.otu.pct = 0.1;
# Confounding effect
grp.ratio = 1; conf.cov.cor = 0.6; Z.diff.otu.pct = 0.05; MgZ = 1; SgZ = 0;Z.nondiff.otu.pct = 0.1;
# Sequence depth
depth.mu = 10000; depth.theta = 5;  depth.conf.factor = 0

prop <- apply(otu.tab, 2, function(x) x/sum(x))
idx <- names(sort(rowMeans(prop), decreasing = T))[1:500]
prop <- prop[idx,, drop=F]
library(circlize)

sub.count <- otu.tab[,c(1:2)]
col_fun = colorRamp2(c(0, mean(sub.count),max(sub.count)), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
colnames(sub.count) = c('S1','S2')
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/observedcount1.pdf", width = 1, height = 10)
ComplexHeatmap::Heatmap(sub.count, col = col_fun,
                        show_row_dend = F, show_column_dend = F, name = " ",
                        show_heatmap_legend = T, show_row_names = F, 
                        show_column_names =T)
dev.off()


sub.count1 <- otu.tab[,c((ncol(otu.tab)-1):ncol(otu.tab))]
colnames(sub.count1) = c('S294','S295')
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/observedcount2.pdf", width = 2, height = 10)
ComplexHeatmap::Heatmap(sub.count1, col = col_fun,show_row_dend = F, show_column_dend = F, 
                        name = " ",show_heatmap_legend = T, show_row_names = F, show_column_names = T)
dev.off()


col_fun = colorRamp2(c(0, mean(otu.tab),max(otu.tab)), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/observedcount.pdf", width = 10, height = 10)
ComplexHeatmap::Heatmap(otu.tab, col = col_fun,show_row_dend = F, show_column_dend = F, 
                        name = " ",
                        show_heatmap_legend = T, show_row_names = F, show_column_names =F)
dev.off()




col_fun = colorRamp2(c(0, mean(prop[,c(1:2)]),max(prop[,c(1:2)])), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/input_prop1.pdf", width = 1, height = 10)
ComplexHeatmap::Heatmap(prop[,c(1:2)], col = col_fun,show_row_dend = F, name = " ",show_column_dend = F, show_heatmap_legend = F, show_row_names = F, show_column_names = F)
dev.off()

col_fun = colorRamp2(c(0, mean(prop[,ncol(prop)]),max(prop[,ncol(prop)])), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/input_prop2.pdf", width = 1, height = 10)
ComplexHeatmap::Heatmap(prop[,c(ncol(prop)), drop =F],name = " ", col = col_fun,show_row_dend = F, show_column_dend = F, show_heatmap_legend = F, show_row_names = F, show_column_names = F)
dev.off()
### select otu.tab and sample firstly, then estimate parameters
# otu.tab.p <- apply(otu.tab, 2, function(x) x/sum(x))
# sample.names <- colnames(otu.tab.p)
# idx.sample <- sample(sample.names, nSubject)
# OTU <- otu.tab[,idx.sample, drop=F]
# otu.names <- names(sort(rowMeans(otu.tab.p), decreasing = T))[1:nOTU]
# otu_tab <- OTU[otu.names,, drop=F]
# model.paras <- EstPara0(otu_tab)
# otu_tab <- model.paras$otu.tab
# rownames(otu_tab) <- paste0('taxon',1:nrow(otu_tab))
# colnames(otu_tab) <- paste0('S',1:ncol(otu_tab))

## Estimated parameters
model.paras <- EstPara(otu.tab, dirmult.paras = dirmult.paras)

# col_fun = colorRamp2(c(min(model.paras$otu.tab), mean(prop),max(prop)), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
# pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_stepsabsolute_abundance.pdf", width = 10, height = 10)
# ComplexHeatmap::Heatmap(model.paras$otu.tab, col = col_fun, show_row_dend = F, show_column_dend = F, show_heatmap_legend = T, 
#                         name = " ",show_row_names = F, show_column_names = F)
# dev.off()

col_fun = colorRamp2(c(min(model.paras$otu.tab[,c(1:2)]), mean(model.paras$otu.tab[,c(1:2)]),max(model.paras$otu.tab[,c(1:2)])), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/absolute_abundance1.pdf", width = 1, height = 10)
ComplexHeatmap::Heatmap(model.paras$otu.tab[,c(1:2)], col = col_fun, show_row_dend = F, show_column_dend = F, show_heatmap_legend = T, 
                        name = " ",show_row_names = F, show_column_names = F)
dev.off()

pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/absolute_abundance2.pdf", width = 1, height = 10)
ComplexHeatmap::Heatmap(model.paras$otu.tab[,c(264:265)], col = col_fun, show_row_dend = F, show_column_dend = F, show_heatmap_legend = T, 
                        name = " ",show_row_names = F, show_column_names = F)
dev.off()



OTU <- model.paras$otu.tab[(1 : (nOTU)),] # otu * sample
sample.names <- colnames(OTU)
otu.names <- rownames(OTU)
set.seed(123)
idx.sample <- sample(sample.names, nSubject)
otu_tab <- OTU[,idx.sample, drop=F] # otu * sample, template otu table: each column represents template sample for each Subject

#### select nSubject first then top nOTU
# otu0 <- model.paras$otu.tab
# rownames(otu0) <- paste0('taxon',1:nrow(otu0))
# colnames(otu0) <- paste0('S',1:ncol(otu0))
# sample.names <- colnames(otu0)
# idx.sample <- sample(sample.names, nSubject)
# OTU <- otu0[,idx.sample, drop=F]
# otu.names <- names(sort(rowMeans(OTU),decreasing=T))[1:nOTU]
# otu_tab <- OTU[otu.names,, drop=F]
# otu_tab.prop <- apply(otu_tab, 2, function(x) x/sum(x))

# m0 <- reshape2::melt(otu0)
# colnames(m0) <- c('otu','sample','value')
# p0 = ggplot(m0, aes(sample, otu, fill = value)) + 
#   geom_tile(color = 'black') + 
#   scale_fill_gradient(low = "white", high = "forestgreen") + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   labs(x = '', y = '')
# p0


colnames(otu_tab) <- paste0('S',1:ncol(otu_tab))
rownames(otu_tab) <- paste0('taxon',1:nrow(otu_tab))

col_fun = colorRamp2(c(min(otu_tab), mean(otu_tab),max(otu_tab)), c(brewer.pal(9,'Greens')[1],c(brewer.pal(9,'Greens')[7], brewer.pal(9,'Greens')[9])))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/baseline.pdf", width = 3, height = 3)
ComplexHeatmap::Heatmap(otu_tab, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(otu_tab),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


## Replicate template sample in each Subject, generate baseline otu table
otu.tab <- otu_tab[, rep(1:nSubject, each=nTime)]
colnames(otu.tab) <- paste0(paste0('R',1:nTime), ':',rep(paste0('S',1:nSubject), each = nTime))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/baseline_rep.pdf", width = 3, height = 3)
ComplexHeatmap::Heatmap(otu.tab, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(otu.tab),
                        column_order = colnames(otu.tab),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


#[have you tested witch alpha is reasonable]
prop.max <- otu.tab[which.max(rowMeans(otu.tab)),]
adj.para <- t(apply(otu.tab,1,function(x) (prop.max / x)^alpha))

## Generate meta data 
nSam <- ncol(otu.tab)
SubjectID <- as.numeric(gsub('.*S','',colnames(otu.tab)))
time <- c(0.5,0.5,0.5,0.8,0.8,0.8)#as.numeric(as.factor(gsub(':.*','',colnames(otu.tab))))-1 
if(matched.pair){
  X <- time
}else{
  X <- as.vector(ifelse(SubjectID <= quantile(unique(SubjectID), grp.ratio / (1 + grp.ratio)), 0, 1))
}
rho <- sqrt(conf.cov.cor ^ 2 / (1 - conf.cov.cor ^ 2))
Z <- c(0.6,0.6,0.6,-0.6,-0.6,-0.6)#as.vector(rho * scale(X) + rnorm(nSam))
meta <- as.data.frame(cbind(X, Z, SubjectID, time))
rownames(meta) <- colnames(otu.tab)
# meta <- meta %>% dplyr::select(-Z)
# meta$Z <- paste0('z',1:nrow(meta))
meta$Z <- round(meta$Z,2)

library(gtable)
library(gridExtra)
# dev.off()
# pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/meta.pdf")
# g <- tableGrob(meta)
# tg <- grid.draw(g)
# dev.off()
# 
# pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/Xv.pdf")
# g <- tableGrob(meta[,'X', drop =F])
# tg <- grid.draw(g)
# dev.off()
# 
# 
# pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/Tv.pdf")
# g <- tableGrob(meta[,'time', drop =F])
# tg <- grid.draw(g)
# dev.off()
# 
# pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/SubjectIDv.pdf")
# g <- tableGrob(meta[,'SubjectID', drop =F])
# tg <- grid.draw(g)
# dev.off()
# 
# 
# pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/Zv.pdf")
# g <- tableGrob(meta[,'Z', drop =F])
# tg <- grid.draw(g)
# dev.off()



## X
X.m <- t(replicate(nOTU, X))
colnames(X.m) <- colnames(otu.tab)
rownames(X.m) <- (rownames(otu.tab))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/X.pdf")
g <- tableGrob(X.m)
tg <- grid.draw(g)
dev.off()
## time
time.m <- t(replicate(nOTU, time))
colnames(time.m) <- colnames(otu.tab)
rownames(time.m) <- (rownames(otu.tab))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/time.pdf")
g <- tableGrob(time.m)
tg <- grid.draw(g)
dev.off()


Z.m <- t(replicate(nOTU, Z))
colnames(Z.m) <- colnames(otu.tab)
rownames(Z.m) <- (rownames(otu.tab))
Z.m <- round(Z.m, 2)
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/Z.pdf")
g <- tableGrob(Z.m)
tg <- grid.draw(g)
dev.off()



# XT.m <- t(replicate(nOTU, X *time))
# colnames(XT.m) <- colnames(otu.tab)
# rownames(XT.m) <- rev(rownames(otu.tab))
# XT.m <- round(XT.m, 2)
# pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/XT.pdf")
# g <- tableGrob(XT.m)
# tg <- grid.draw(g)
# dev.off()



## SubjectID
# Subject.m <- t(replicate(nOTU, SubjectID))
# colnames(Subject.m) <- colnames(otu.tab)
# rownames(Subject.m) <- rownames(otu.tab)
# pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/Subject.pdf")
# g <- tableGrob(Subject.m)
# tg <- grid.draw(g)
# dev.off()





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
#interaction.diff.otu.ind <- sample(unique(c(time.diff.otu.ind, X.diff.otu.ind)), interaction.diff.otu.num)
interaction.diff.otu.ind <- c(interaction.diff.otu.ind, sample(otu.ord, interaction.diff.otu.num))

## Generate coefficient for X
if(balanced.X){
  coef.X <- sample(c(rnorm(floor(nOTU / 2), mean = -MgX, sd = SgX), rnorm(nOTU - floor(nOTU / 2), mean = MgX, sd = SgX)))
}else{
  coef.X <- rnorm(nOTU, MgX, SgX)
}

coef.X[setdiff(otu.ord, X.diff.otu.ind)] <- 0 

coef.X.m <- replicate(nSam, coef.X)
rownames(coef.X.m) <- rownames(otu.tab)
colnames(coef.X.m) <- rev(colnames(otu.tab))
dev.off()
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/coef.X.pdf")
# coef.X.m[coef.X.m !=0] <- 'x'
g <- tableGrob(coef.X.m)
tg <- grid.draw(g)
dev.off()



eta.X <- coef.X  %*% t(scale(X)) * adj.para 
eta.X <- round(eta.X,digits = 2)
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

# exp(abs(min(eta.X))+max(eta.X))


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
rownames(coef.T) <- rownames(otu.tab)
colnames(coef.T) <- colnames(otu.tab)
dev.off()
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/coef.T.pdf")
# coef.X.m[coef.X.m !=0] <- 'x'
g <- tableGrob(round(coef.T,2))
tg <- grid.draw(g)
dev.off()


# eta.time <- t(t(coef.T) * as.vector(scale(time))) * adj.para
eta.time <- t(t(coef.T) * as.vector(time)) * adj.para
col_fun = colorRamp2(c(min(eta.time), 0,max(eta.time)), c(brewer.pal(9,'Blues')[7],'white', brewer.pal(9,'Reds')[8]))
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
coef.Z <- sample(c(rnorm(floor(nOTU / 2), mean = -MgZ, sd = SgZ), rnorm(nOTU - floor(nOTU / 2), mean = MgZ, sd = SgZ)))
coef.Z[setdiff(otu.ord, Z.diff.otu.ind)] <- 0 

coef.Z.m <- replicate(nSam, coef.Z)
rownames(coef.Z.m) <- rownames(otu.tab)
colnames(coef.Z.m) <- colnames(otu.tab)
dev.off()
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/coef.Z.pdf")
g <- tableGrob(round(coef.Z.m,2))
tg <- grid.draw(g)
dev.off()




eta.Z <- coef.Z %*% t(scale(Z)) * adj.para
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



# ## Generate coefficient for interaction 
# if(balanced.XT){
#   coef.interaction <- sample(c(rnorm(floor(nOTU / 2), mean = -MgXT, sd = SgXT), rnorm(nOTU - floor(nOTU / 2), mean = MgXT, sd = SgXT)))
# }else{
#   coef.interaction  <- rnorm(nOTU, MgXT, SgXT)
# }
# coef.interaction[setdiff(otu.ord, interaction.diff.otu.ind)] <- 0 
# coef.interaction.Subject <- replicate(nSubject, rnorm(nOTU, coef.interaction, 0)) 
# coef.interaction <- matrix(data = apply(coef.interaction.Subject, 2, function(x) rep(x, nTime)), ncol = ncol(coef.interaction.Subject)*nTime)
# 
# rownames(coef.interaction) <- rownames(otu.tab)
# colnames(coef.interaction) <- colnames(otu.tab)
# dev.off()
# pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/coef.interaction.pdf")
# g <- tableGrob(round(coef.interaction,2))
# tg <- grid.draw(g)
# dev.off()


# eta.interaction <-  t(t(coef.interaction) * as.vector(scale(time) *scale(X))) * adj.para 
# m8 <- reshape2::melt(eta.interaction)
# colnames(m8) <- c('otu','sample','value')
# p8 <- ggplot(m8, aes(sample, otu, fill = value)) + 
#   geom_tile(color = 'black') + 
#   scale_fill_gradient2(mid = 'white', high = brewer.pal(9,'Greens')[9]) +
#   theme(axis.text.x = element_text(angle = 90, color = 'black', size = 14, vjust = 0.4),
#         axis.text.y = element_text(color = 'black', size = 14), 
#         legend.text = element_text(color = 'black', size = 14),
#         legend.title = element_text(color = 'black', size = 14)) + 
#   labs(x = '', y = '', fill = 'value')
# p8
# ggsave("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/eta.interaction.pdf", width = 3.5, height = 3)




## combine index matrix, coefficient matrix and design matrix
# eta.exp <- exp(t(eta.X + eta.Z + eta.time + eta.interaction + error))
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



eta.exp <- exp(t(eta.exp)) * t(otu.tab)
not.trans.eta.exp <- eta.exp

# Renormalize
otu.tab.prop <- eta.exp / rowSums(eta.exp)
otu.tab.prop <- t(otu.tab.prop) 
col_fun = colorRamp2(c(min(otu.tab.prop), median(otu.tab.prop),max(otu.tab.prop)), c(brewer.pal(9,'Blues')[1],brewer.pal(9,'Reds')[3], brewer.pal(9,'Reds')[8]))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/eta.pdf", width = 3.5, height = 3)
ComplexHeatmap::Heatmap(otu.tab.prop, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(otu.tab.prop),
                        column_order = colnames(otu.tab.prop),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


# Generate the sequence depth
nSeq <- rnegbin(ncol(otu.tab.prop), mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta)
otu.tab.sim <- sapply(1:ncol(otu.tab.prop), function (i) rmultinom(1, nSeq[i], otu.tab.prop[,i]))
colnames(otu.tab.sim) <- rownames(eta.exp)
rownames(otu.tab.sim) <- rownames(otu.tab)
col_fun = colorRamp2(c(min(otu.tab.sim), median(otu.tab.sim),max(otu.tab.sim)), c(brewer.pal(9,'Reds')[1],brewer.pal(9,'Reds')[3], brewer.pal(9,'Reds')[8]))
pdf("/research/bsi/projects/staff_analysis/m216453/correlated/Sim_steps/otu.tab.sim.pdf", width = 3.5, height = 3)
ComplexHeatmap::Heatmap(otu.tab.sim, col = col_fun, show_row_dend = F, show_column_dend = F, 
                        row_order = rownames(otu.tab.sim),
                        column_order = colnames(otu.tab.sim),
                        show_heatmap_legend = T, cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        border = TRUE,column_gap = unit(1, "mm"), rect_gp = gpar(col= "black"),
                        name = " ",show_row_names = T, show_column_names = T)
dev.off()


