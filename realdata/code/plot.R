## LDM  error: Error in check(sn, control = control, quietly = quietly) :Design must be balanced if permuting 'strata'.
## ZIBR does not support unbalanced design
## smoker;Nicolas2013;IBD2017
setwd('/research/bsi/projects/staff_analysis/m216453/correlated/Data/realshuffle/')
pttn <- c('smoker','Nicholas2013','IBD2017')

sig.dfs <- original.dfs <- sig.iterss <- list()
y.hlines <- c()
for(pttn in c('smoker','Nicholas2013','IBD2017')){
  if(pttn == 'IBD2017'){
    methods <- c('GLMMPQL','LinDA','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMB','MaAsLin2','glmernb')#
  }else{
    methods <- c('GLMMPQL','LinDA','NBMM','ZIGMM','ZINBMM','glmmadaptive','glmmTMB','MaAsLin2','glmernb','ZIBR','LDM')#
  }
  
  files <- list.files(pattern = paste0('*',pttn,'.res.Rdata'))
  for(method in methods){
    cat(method,': ')
    ld <- files[grep(paste0('^',method),files)]
    cat(length(ld),'\n')
  }
  
  sig <-  NULL
  for(method in methods){
    cat(method,'\n')
    ld <- files[grep(paste0('^',method),files)]
    for(f in ld){
      wrapper.obj <- NULL
      load(f)
      sig0 <- wrapper.obj$sig
      if(pttn == 'IBD2017'){
        names(sig0) <- paste0(f,'-',names(sig0))
      }else{
        names(sig0) <- f
      }
      sig <- c(sig,sig0)
    }
  }
  
  save(sig, file = paste0('/research/bsi/projects/staff_analysis/m216453/correlated/Data/sig_',pttn, '.Rdata'))
  
  if(pttn == 'IBD2017'){
    sig.df <- as.data.frame(cbind(sig)) %>% rownames_to_column('name') %>% 
      separate(name,into = c('method','iter','dataset'), sep ='_') %>% 
      separate(dataset,into = c('dataset','test'), sep ='-')%>% 
      mutate(dataset = gsub('.res.Rdata','', dataset))
    sig.df$n <- ifelse(sig.df$sig > 0,1,0)
    sig.iters <- aggregate(n ~ method + test, data = sig.df, function(x) mean(x)) %>% arrange(n)
    sig.df.sum <- aggregate(sig ~ method +test , data = sig.df, function(x) median(x)) %>% arrange(sig)
  }else{
    sig.df <- as.data.frame(cbind(sig)) %>% rownames_to_column('name') %>% 
      separate(name,into = c('method','iter','dataset'), sep ='_') %>% mutate(dataset = gsub('.res.Rdata','', dataset))
    sig.df$n <- ifelse(sig.df$sig >0,1,0)
    sig.iters <- aggregate(n ~ method, data = sig.df, function(x) mean(x)) %>% arrange(n)
    sig.df.sum <- aggregate(sig ~ method  , data = sig.df, function(x) median(x)) %>% arrange(sig)
  }
  
  sig.df <- within(sig.df, method <- factor(method, levels=unique(sig.df.sum$method)))
  
  real.files <- list.files('/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/',pattern = paste0(pttn,'.res.Rdata'))
  original <- c()
  for(real.file in real.files){
    wrapper.obj <- NULL
    load(paste0('/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/',real.file))
    sig0 <- wrapper.obj$sig
    if(pttn == 'IBD2017'){
      names(sig0) <- paste0(real.file,'-',names(sig0))
    }else{
      names(sig0) <- real.file
    }
    original <- c(original, sig0)
  }
  if(pttn == 'IBD2017'){
    original.df <- as.data.frame(original) %>% rownames_to_column('method') %>% 
      mutate(method = gsub(paste0('_',pttn,'.res.Rdata'),'',method)) %>% 
      separate(method,into = c('method','test'), sep ='-')
  }else{
    original.df <- as.data.frame(original) %>% rownames_to_column('method') %>% mutate(method = gsub(paste0('_',pttn,'.res.Rdata'),'',method))
  }
  y.hline <- nrow(wrapper.obj$res) * 0.05
  
  y.hlines <- c(y.hlines, y.hline)
  original.dfs[[pttn]] <- original.df
  sig.dfs[[pttn]] <- sig.df
  sig.iterss[[pttn]] <- sig.iters
}


# save(y.hlines,original.dfs, sig.dfs,sig.iterss, file= '/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/shuffle.sig1.Rdata')
load('/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/shuffle.sig1.Rdata')
## smoker:197; Nicholas: 110; IBD: 498
sig.df.all <- rbind(sig.dfs$smoker %>% mutate(grp = 'repeated-measure', test = 'X', total =197), 
                    sig.dfs$Nicholas2013 %>% mutate(grp = 'matched-pair', test = 'time', total = 110)) %>% 
  rbind(sig.dfs$IBD2017 %>% mutate(test = gsub('.padj','',test), total = 498) %>% mutate(grp = 'longitudinal') %>% dplyr::filter(test !='time'))  %>% 
  mutate(pct = sig/total)
head(sig.df.all)
sig.df.all$grp <- paste0(sig.df.all$grp,'(',sig.df.all$test,')')
sig.df.all <- within(sig.df.all, grp <- factor(grp, levels=c("repeated-measure(X)", "matched-pair(time)","longitudinal(X)","longitudinal(interaction)")))
library(RColorBrewer)
cols <- c(brewer.pal(8,'Set1')[c(1:5,7:9)],brewer.pal(3,'Set2'))
names(cols) <- as.character(unique(sig.df.all$method))
# save(cols, file = '../cols.Rdata')
sig.df.all$dataset <- gsub('^smoker$','Smoker2010',sig.df.all$dataset)
sig.df.all$dataset <- gsub('^Nicholas2013$','Nicholas2013',sig.df.all$dataset)
sig.df.all[sig.df.all$dataset=='IBD2017' & sig.df.all$test=='interaction','dataset'] <- 'IBD2017(interaction)'
sig.df.all[sig.df.all$dataset=='IBD2017' & sig.df.all$test=='X','dataset'] <- 'IBD2017'

sig.df.all1 <- sig.df.all %>% dplyr::filter(sig.df.all$grp == "repeated-measure(X)") %>% 
  group_by(dataset,method, test) %>% mutate(md = median(sig)) %>% 
  group_by(dataset,method, test) %>% arrange(md)
sig.df.all1 <- within(sig.df.all1, method <- factor(method, levels=as.character(unique(sig.df.all1$method))))

sig.df.all11 <- sig.df.all1%>% dplyr::filter(method %in% c('LDM','GLMMPQL','LinDA','MaAsLin2','NBMM','ZINBMM'))
summary(sig.df.all11$pct)
data_summary(data = sig.df.all11, formula = 'pct ~ method')
sig.df.all12 <- sig.df.all1%>% dplyr::filter(method %in% c('ZIGMM'))
summary(sig.df.all12$pct)
data_summary(data = sig.df.all12, formula = 'pct ~ method')
sig.df.all13 <- sig.df.all1%>% dplyr::filter(!(method %in% c('LDM','GLMMPQL','LinDA','MaAsLin2','NBMM','ZINBMM','ZIGMM')))
summary(sig.df.all13$pct)
data_summary(data = sig.df.all13, formula = 'pct ~ method')



sig.df.all2 <- sig.df.all %>% dplyr::filter(sig.df.all$grp == "matched-pair(time)") %>% 
  group_by(dataset,method, test) %>% mutate(md = median(sig)) %>% 
  group_by(dataset,method, test) %>% arrange(md)
sig.df.all2 <- within(sig.df.all2, method <- factor(method, levels=as.character(unique(sig.df.all2$method))))

sig.df.all21 <- sig.df.all2%>% dplyr::filter(method %in% c('LDM','LinDA'))
summary(sig.df.all21$pct)

sig.df.all3 <- sig.df.all %>% dplyr::filter(sig.df.all$grp == "longitudinal(X)") %>% 
  group_by(dataset,method, test) %>% mutate(md = median(sig)) %>% 
  group_by(dataset,method, test) %>% arrange(md)
sig.df.all3 <- within(sig.df.all3, method <- factor(method, levels=as.character(unique(sig.df.all3$method))))
sig.df.all31 <- sig.df.all3 %>% dplyr::filter(method %in% c('LinDA'))
summary(sig.df.all31$pct)


sig.df.all4 <- sig.df.all %>% dplyr::filter(sig.df.all$grp == "longitudinal(interaction)") %>% 
  group_by(dataset,method, test) %>% mutate(md = median(sig)) %>% 
  group_by(dataset,method, test) %>% arrange(md)
sig.df.all4 <- within(sig.df.all4, method <- factor(method, levels=as.character(unique(sig.df.all4$method))))
sig.df.all41 <- sig.df.all4 %>% dplyr::filter(method %in% c('LinDA'))
summary(sig.df.all41$pct)


p1 <- ggplot(sig.df.all1) + 
  geom_boxplot(aes(x = method, y = sig, fill = method), outlier.size = 0.5) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~./1.97, name="% of differential taxa from permuted data"))+
  # facet_wrap(.~grp, scales = 'free') + 
  scale_fill_manual(values = cols) + 
  facet_wrap(.~ dataset, scales = 'free') + 
  theme_bw() + 
  labs(x = '', y = '# of differential taxa from permuted data', fill = '') + 
  theme(axis.text.x = element_text(color = 'black', size = 14, angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none')
p2 <- ggplot(sig.df.all2) + 
  geom_boxplot(aes(x = method, y = sig, fill = method), outlier.size = 0.5) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~./1.1, name="% of differential taxa from permuted data"))+
  scale_fill_manual(values = cols) + 
  facet_wrap(.~ dataset, scales = 'free') + 
  theme_bw() + 
  labs(x = '', y = '# of differential taxa from permuted data', fill = '') + 
  theme(axis.text.x = element_text(color = 'black', size = 14, angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none')
p3 <- ggplot(sig.df.all3) + 
  geom_boxplot(aes(x = method, y = sig, fill = method), outlier.size = 0.5) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~./4.98, name=""))+
  scale_fill_manual(values = cols) + 
  facet_wrap(.~ dataset, scales = 'free') + 
  theme_bw() + 
  labs(x = '', y = '# of differential taxa from permuted data', fill = '') + 
  theme(axis.text.x = element_text(color = 'black', size = 14, angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none')
p4 <- ggplot(sig.df.all4) + 
  geom_boxplot(aes(x = method, y = sig, fill = method), outlier.size = 0.5) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~./4.98, name="% of differential taxa from permuted data"))+
  scale_fill_manual(values = cols) + 
  facet_wrap(.~ dataset, scales = 'free') + 
  theme_bw() + 
  labs(x = '', y = '', fill = '') + 
  theme(axis.text.x = element_text(color = 'black', size = 14, angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none')
require(grid)   # for the textGrob() function
require(ggpubr)
pp <- ggpubr::ggarrange(p1 + rremove("ylab") , p2 + rremove("ylab"), p3, p4, # remove axis labels from plots
                    labels = NULL, ncol = 2, nrow = 2,align = "hv")



sig.iters <- rbind(sig.iterss$smoker %>% mutate(grp = 'repeated-measure', test = 'X'), sig.iterss$Nicholas2013 %>% mutate(grp = 'matched-pair', test = 'time')) %>% 
  rbind(sig.iterss$IBD2017 %>% mutate(test = gsub('.padj','',test))%>% mutate(grp = 'longitudinal') %>% dplyr::filter(test !='time')) 
sig.iters$grp <- paste0(sig.iters$grp,'(',sig.iters$test,')')
sig.iters$grp <- gsub('repeated-measure\\(X\\)','Smoker2010',sig.iters$grp)
sig.iters$grp <- gsub('matched-pair\\(time\\)','Nicholas2013',sig.iters$grp)
sig.iters$grp <- gsub('longitudinal\\(X\\)','IBD2017',sig.iters$grp)
sig.iters$grp <- gsub('longitudinal\\(interaction\\)','IBD2017(interaction)',sig.iters$grp)
# sig.iters <- within(sig.iters, grp <- factor(grp, levels=c("repeated-measure(X)", "matched-pair(time)","longitudinal(X)","longitudinal(time)")))
sig.iters <- within(sig.iters, grp <- factor(grp, levels=c("Smoker2010", "Nicholas2013","IBD2017","IBD2017(interaction)")))

pp2 <- ggplot(sig.iters) + geom_bar(aes(x = reorder(method,n), y = n*100, fill = method), stat = 'identity') + theme_bw() + 
  facet_wrap(.~ grp, scales = 'free') + 
  geom_hline(yintercept = 5, linetype = 3, color = 'red') + 
  scale_fill_manual(values = cols) + 
  theme_bw() + 
  guides(fill=  F) + 
  labs(x = '', y = '% of permuted datasets with positive findings', fill = '') + 
  theme(axis.text.x = element_text(color = 'black', size = 14, angle = 90, hjust = 1),
        axis.text = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none')
library(ggpubr)
figure <- ggarrange(pp, pp2, nrow = 1, common.legend = F)
ggsave(file= paste0('/research/bsi/projects/staff_analysis/m216453/correlated/plot/combined.pdf'), width = 15, height = 8)



RES <- list()
for(pttn in c('smoker','Nicholas2013','IBD2017')){
  cat('[ ',pttn,' ]\n')
  real.files <- list.files('/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/',pattern = paste0(pttn,'.res.Rdata'))
  real.files <- real.files[!(real.files %in% c("IFAA_smoker.res.Rdata","IFAA_Nicholas2013.res.Rdata","IFAA_IBD2017.res.Rdata"))]
  wrapper.obj <- NULL
  load(paste0('/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/',real.files[1]))
  res00 <- wrapper.obj$res %>% rownames_to_column('taxa')
  if(pttn == 'IBD2017'){
    res00 <-  res00 %>% dplyr::select(taxa, X.padj, interaction.padj)
    colnames(res00)[2:3] <- paste0(gsub('.res.Rdata','',real.files[1]),'-',colnames(res00)[2:3])
  }else{
    colnames(res00)[2] <- gsub('.res.Rdata','',real.files[1])
  }
  
  for(real.file in real.files){
    cat(real.file,'\n')
    wrapper.obj <- NULL
    load(paste0('/research/bsi/projects/staff_analysis/m216453/correlated/Data/real/',real.file))
    if(length(grep('MaAsLin2',real.file))==1){
      res0 <- wrapper.obj$res %>% rownames_to_column('taxa')
      res0$taxa <- gsub('X','',res0$taxa)
    }else{
      res0 <- wrapper.obj$res %>% rownames_to_column('taxa')
    }
    
    if(pttn == 'IBD2017'){
      res0 <- res0 %>% dplyr::select(taxa, X.padj, interaction.padj)
      colnames(res0)[2:3] <- paste0(gsub('.res.Rdata','',real.file),'-',colnames(res0)[2:3])
    }else{
      colnames(res0)[2] <- gsub('.res.Rdata','',real.file)
    }
    res00 <- full_join(res00, res0)
  }
  RES[[pttn]] <- res00 %>% column_to_rownames('taxa')
}


### line plot for # of discoveries on the raw dataset
res00 <- RES$smoker
res00[is.na(res00)] <- 1
colnames(res00) <- gsub('_smoker','',colnames(res00))
smoker.fdr <- list()
for(fdr in c(0,0.05, 0.1,0.15, 0.2)){
  smoker.fdr[[paste0(fdr)]] <- apply(res00, 2, function(x) sum(x<=fdr))
}

res00 <- RES$Nicholas2013
res00[is.na(res00)] <- 1
colnames(res00) <- gsub('_Nicholas2013','',colnames(res00))
Nicholas.fdr <- list()
for(fdr in c(0,0.05, 0.1,0.15, 0.2)){
  Nicholas.fdr[[paste0(fdr)]] <- apply(res00, 2, function(x) sum(x<=fdr))
}

res00 <- RES$IBD2017
res00[is.na(res00)] <- 1
colnames(res00) <- gsub('_IBD2017','',colnames(res00))
IBD2017.fdr <- list()
for(fdr in c(0,0.05, 0.1,0.15, 0.2)){
  IBD2017.fdr[[paste0(fdr)]] <- apply(res00, 2, function(x) sum(x<=fdr))
}

smoker.X = plyr::ldply(smoker.fdr, rbind) %>% dplyr::rename(fdr = '.id')
Nicholas.T = plyr::ldply(Nicholas.fdr, rbind)%>% dplyr::rename(fdr = '.id')
IBD2017.fdr = plyr::ldply(IBD2017.fdr, rbind)%>% dplyr::rename(fdr = '.id')
IBD2017.X <- IBD2017.fdr[,c(1,grep('-X.padj',colnames(IBD2017.fdr)))]
IBD2017.interaction <- IBD2017.fdr[,c(1,grep('-interaction.padj',colnames(IBD2017.fdr)))]
colnames(IBD2017.interaction) <- gsub('-interaction.padj','',colnames(IBD2017.interaction))
colnames(IBD2017.X) <- gsub('-X.padj','',colnames(IBD2017.X))


# smoker.X.m <- melt(smoker.X) %>% mutate(grp = 'Smoker2010(X)')
# Nicholas.T.m <- melt(Nicholas.T) %>% mutate(grp = 'Nicholas2013(time)')
# IBD2017.X.m <- melt(IBD2017.X) %>% mutate(grp = 'IBD2017(X)')
# IBD2017.interaction.m <- melt(IBD2017.interaction) %>% mutate(grp = 'IBD2017(interaction)')
# df <- rbind(smoker.X.m, Nicholas.T.m, IBD2017.X.m, IBD2017.interaction.m)
# df$fdr <- as.numeric(df$fdr)
# df <- within(df, grp <- factor(grp, levels = c('Smoker2010(X)','Nicholas2013(time)','IBD2017(X)','IBD2017(interaction)')))
# ggplot(df, aes(x = fdr, y = value, group = variable, color = variable)) + 
#   geom_point()+
#   geom_line() + 
#   facet_wrap(~grp, scales = 'free_y') + 
#   scale_color_manual(values = cols) + 
#   # scale_x_continuous(expand = c(0, 0)) + 
#   theme_bw()+labs(x = '', y = 'Number of Discoveries', color = '') + 
#   theme(axis.text.x = element_text(color = 'black', size = 14),
#         axis.text = element_text(color = 'black', size = 14),
#         strip.text = element_text(color = 'black', size = 14),
#         axis.title = element_text(color = 'black', size = 14),
#         legend.title = element_text(color = 'black', size = 14),
#         legend.text = element_text(color = 'black', size = 14),
#         legend.position = 'right') 
# ggsave(file = '/research/bsi/projects/staff_analysis/m216453/correlated/plot/barplot_discoveries.pdf',
#        width = 8, height = 6)


### Upset plot
library(ComplexHeatmap)
res00 <- RES$smoker
res00[is.na(res00)] <- 1
colnames(res00) <- gsub('_smoker','',colnames(res00))
sum(apply(res00, 1, function(x) sum(x<=0.05)) >0)/nrow(res00)
sum(apply(res00, 1, function(x) sum(x<=0.1)) >0)/nrow(res00)
sum(apply(res00, 1, function(x) sum(x<=0.2)) >0)/nrow(res00)
res00.lt.05 <- apply(res00, 2, function(x) names(x[x<=0.05]))
res00.lt.1 <- apply(res00, 2, function(x) names(x[x<=0.1]))
res00.lt.2 <- apply(res00, 2, function(x) names(x[x<=0.2]))
pdf('/research/bsi/projects/staff_analysis/m216453/correlated/plot/smoker.pdf', width = 6, height = 5)
res00.lt.2.m1 = make_comb_mat(res00.lt.05)
UpSet(res00.lt.2.m1)
res00.lt.2.m1 = make_comb_mat(res00.lt.1)
UpSet(res00.lt.2.m1)
res00.lt.2.m1 = make_comb_mat(res00.lt.2)
UpSet(res00.lt.2.m1)
dev.off()

pdf('/research/bsi/projects/staff_analysis/m216453/correlated/plot/Nicholas2013.pdf', width = 8, height = 6)
res00 <- RES$Nicholas2013
res00[is.na(res00)] <- 1
colnames(res00) <- gsub('_Nicholas2013','',colnames(res00))
sum(apply(res00, 1, function(x) sum(x<=0.05)) >0)/nrow(res00)
sum(apply(res00, 1, function(x) sum(x<=0.1)) >0)/nrow(res00)
sum(apply(res00, 1, function(x) sum(x<=0.2)) >0)/nrow(res00)

res00.lt.05 <- apply(res00, 2, function(x) names(x[x<=0.05]))
res00.lt.1 <- apply(res00, 2, function(x) names(x[x<=0.1]))
res00.lt.2 <- apply(res00, 2, function(x) names(x[x<=0.2]))
res00.lt.2.m1 = make_comb_mat(res00.lt.05)
UpSet(res00.lt.2.m1)
res00.lt.2.m1 = make_comb_mat(res00.lt.1)
UpSet(res00.lt.2.m1)
res00.lt.2.m1 = make_comb_mat(res00.lt.2)
UpSet(res00.lt.2.m1)
dev.off()


pdf('/research/bsi/projects/staff_analysis/m216453/correlated/plot/IBD2017_X.pdf', width = 12, height = 5)
res00 <- RES$IBD2017[,c(colnames(RES$IBD2017)[grep('X.padj',colnames(RES$IBD2017))])]
res00[is.na(res00)] <- 1
colnames(res00) <- gsub('_IBD2017.*','',colnames(res00))
sum(apply(res00, 1, function(x) sum(x<=0.05)) >0)/nrow(res00)
sum(apply(res00, 1, function(x) sum(x<=0.1)) >0)/nrow(res00)
sum(apply(res00, 1, function(x) sum(x<=0.2)) >0)/nrow(res00)
res00.lt.05 <- apply(res00, 2, function(x) names(x[x<=0.05]))
res00.lt.1 <- apply(res00, 2, function(x) names(x[x<=0.1]))
res00.lt.2 <- apply(res00, 2, function(x) names(x[x<=0.2]))
res00.lt.2.m1 = make_comb_mat(res00.lt.05)
UpSet(res00.lt.2.m1,lwd = 1,pt_size = unit(1.5, "mm"))
res00.lt.2.m1 = make_comb_mat(res00.lt.1)
UpSet(res00.lt.2.m1,lwd = 1,pt_size = unit(1.5, "mm"))
res00.lt.2.m1 = make_comb_mat(res00.lt.2)
UpSet(res00.lt.2.m1,lwd = 1,pt_size = unit(1.5, "mm"))
dev.off()


pdf('/research/bsi/projects/staff_analysis/m216453/correlated/plot/IBD2017_interaction.pdf', width = 6, height = 4)
res00 <- RES$IBD2017[,c(colnames(RES$IBD2017)[grep('interaction.padj',colnames(RES$IBD2017))])]
res00[is.na(res00)] <- 1
colnames(res00) <- gsub('_IBD2017.*','',colnames(res00))
sum(apply(res00, 1, function(x) sum(x<=0.05)) >0)/nrow(res00)
sum(apply(res00, 1, function(x) sum(x<=0.1)) >0)/nrow(res00)
sum(apply(res00, 1, function(x) sum(x<=0.2)) >0)/nrow(res00)
res00.lt.05 <- apply(res00, 2, function(x) names(x[x<=0.05]))
res00.lt.1 <- apply(res00, 2, function(x) names(x[x<=0.1]))
res00.lt.2 <- apply(res00, 2, function(x) names(x[x<=0.2]))
res00.lt.2.m1 = make_comb_mat(res00.lt.05)
UpSet(res00.lt.2.m1)
res00.lt.2.m1 = make_comb_mat(res00.lt.1)
UpSet(res00.lt.2.m1)
res00.lt.2.m1 = make_comb_mat(res00.lt.2)
UpSet(res00.lt.2.m1)
dev.off()




