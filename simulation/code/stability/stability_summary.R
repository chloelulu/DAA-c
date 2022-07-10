## I tried repeated measure balanced(er1) and unbalanced (er1 and er4). er1 is enough.

setwd('/research/bsi/projects/staff_analysis/m216453/correlated/')
library(dplyr)
library(tibble)
## ---- Balanced 
er = 'er1'
nofilter <- list.files(paste0('stability/',er,'/'), pattern = paste0(er,'_summarymatrix*'))
# filter <- list.files('filter/er1/', pattern = 'er1_summarymatrix*')
filter <- list.files(paste0('filter40/',er,'/'), pattern = paste0(er,'_summarymatrix*'))
sum(nofilter==filter)

rhos <- NULL
for(file in nofilter){
  res_seqs <- NULL
  load(paste0('stability/',er,'/',file))
  res_seqs1 <- res_seqs
  res_seqs <- NULL
  load(paste0('filter40/',er,'/',file))
  res_seqs2 <- res_seqs
  if(length(res_seqs1)!=0 & length(res_seqs2)!=0){
    sub <- names(res_seqs1)
    for(i in 1:length(sub)){
      f1 <- res_seqs1[[sub[i]]] %>% dplyr::select(X.padj) %>% rownames_to_column('taxa') %>% dplyr::rename(nofilter = X.padj)
      f2 <- res_seqs2[[sub[i]]]  %>% dplyr::select(X.padj) %>% rownames_to_column('taxa') %>% dplyr::rename(filter = X.padj)
      f <- inner_join(f1, f2)
      rho <- cor.test(f[,2],f[,3], method = 'spearman')$estimate
      names(rho) <- gsub('.*summarymatrix','',gsub('.Rdata','',file))
    }
    rhos <- c(rhos, rho)
  }
}
save(rhos, file= paste0('/research/bsi/projects/staff_analysis/m216453/correlated/rho_',er,'_B.Rdata'))

## Unbalanced
nofilter <- list.files(paste0('stability_U/',er,'/'), pattern = paste0(er,'_summarymatrix*'))
# filter <- list.files('filter/er1/', pattern = 'er1_summarymatrix*')
filter <- list.files(paste0('filter40_U/',er,'/'), pattern = paste0(er,'_summarymatrix*'))
sum(nofilter==filter)

rhos <- NULL
for(file in nofilter){
  res_seqs <- NULL
  load(paste0('stability_U/',er,'/',file))
  res_seqs1 <- res_seqs
  res_seqs <- NULL
  load(paste0('filter40_U/',er,'/',file))
  res_seqs2 <- res_seqs
  if(length(res_seqs1)!=0 & length(res_seqs2)!=0){
    sub <- names(res_seqs1)
    for(i in 1:length(sub)){
      f1 <- res_seqs1[[sub[i]]] %>% dplyr::select(X.padj) %>% rownames_to_column('taxa') %>% dplyr::rename(nofilter = X.padj)
      f2 <- res_seqs2[[sub[i]]]  %>% dplyr::select(X.padj) %>% rownames_to_column('taxa') %>% dplyr::rename(filter = X.padj)
      f <- inner_join(f1, f2)
      rho <- cor.test(f[,2],f[,3], method = 'spearman')$estimate
      names(rho) <- gsub('.*summarymatrix','',gsub('.Rdata','',file))
    }
    rhos <- c(rhos, rho)
  }
}
save(rhos, file= paste0('/research/bsi/projects/staff_analysis/m216453/correlated/rho_',er,'_U.Rdata'))

load(paste0('/research/bsi/projects/staff_analysis/m216453/correlated/rho_',er,'_U.Rdata'))
df1 = as.data.frame(rhos) %>% rownames_to_column('method')
load(paste0('/research/bsi/projects/staff_analysis/m216453/correlated/rho_',er,'_B.Rdata'))
df2 = as.data.frame(rhos) %>% rownames_to_column('method')

df <- rbind(df1 %>% mutate(grp = 'U'), df2 %>% mutate(grp = 'B'))
df$methods <- gsub('.*-','',df$method)
df <- df %>% dplyr::select(-method)
head(df)
ord <- (aggregate(rhos ~ methods, data = df, function(x) median(x)) %>% arrange(rhos))[,1]
df <- within(df,methods <- factor(methods, levels=as.character(ord)))
load('/research/bsi/projects/staff_analysis/m216453/correlated/Data/cols.Rdata')
ggplot(df, aes(x = methods, y = rhos, fill = methods)) + 
  geom_boxplot(stat ='boxplot', outlier.size = 0.5) + 
  scale_fill_manual(values = cols) +
  # facet_grid(.~ grp) + 
  scale_y_continuous(limits = c(min(df$rhos),1),expand = c(0.1, 0.0, 0.1, 0)) +
  theme_bw() + 
  guides(fill = F) + 
  theme(text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, color = 'black'),
        axis.text.y = element_text(color = 'black'),
        legend.text = element_text(size = 16, color = "black"),
        axis.title =  element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(),
        legend.position = 'none',
        panel.grid.minor = element_blank()) + labs(x= '') 
# ggsave('/research/bsi/projects/staff_analysis/m216453/correlated/rhos.pdf', width = 6, height = 6)
ggsave('/research/bsi/projects/staff_analysis/m216453/correlated/plot/rhos20VS0.pdf', width = 6, height = 6)

df_mean <- aggregate(rhos ~ methods, df, function(x) mean(x))
df_mean <- df_mean[order(-df_mean$rhos),]
# save(df_mean, df, file = '/research/bsi/projects/staff_analysis/m216453/correlated/rho.Rdata')
save(df_mean, df, file = '/research/bsi/projects/staff_analysis/m216453/correlated/rho20.Rdata')
