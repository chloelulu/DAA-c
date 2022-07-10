source('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_project/Code/DailyCode.R')
pkg = c('ggforce','vegan','dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringr','scales','tibble','kableExtra','formattable','reactable','htmltools',"htmlwidgets","webshot2")
suppressPackageStartupMessages(sapply(pkg, require, character = T))

clean <- function(res){
  df <- reshape::melt(res)
  colnames(df) <- c('balanced.X','balanced.T','balanced.XT','nSubject','time.n','error.sd','MgX','SgX','X.diff.otu.pct','MgT','SgT','SbT','time.diff.otu.pct',
                    'MgXT','SgXT','interaction.diff.otu.pct','measure','method','value','iter')
  df$MgX <- paste0('MgX=',df$MgX)
  df$SgX <- paste0('SgX=',df$SgX)
  df$MgT <- paste0('MgT=',df$MgT)
  df$SgT <- paste0('SgT=',df$SgT)
  df$MgXT <- paste0('MgXT=',df$MgXT)
  df$SgXT <- paste0('SgXT=',df$SgXT)
  df$error.sd <- paste0('error.sd=',df$error.sd)
  df = df %>% dplyr::filter(!(method %in% c('metamicrobiomeR','glmmadaptiveP','lme','lmer')))#,'NBMM','ZINBMM','ZIGMM'
  return(df)
}


dir = 'repeatedB/';
data1=paste0(dir,'er1_res.Rdata');
load(data1)
f1 = 'MgT'
df = clean(res)%>% mutate(error.sd ='error.sd=1')# %>% filter(!(MgX %in% c('MgX=1','MgX=1.5')))
data <- aggregate(as.formula(paste0('value ~ measure + error.sd +', f1, '+ method')), na.omit(df[df$error.sd=='error.sd=1',]),function(x) mean(x[!is.na(x)])) %>% filter(measure %in% c('FDR','TPR'))
data[data$measure=='FDR',]
data[data$measure=='TPR',]


data2=paste0(dir,'er4_res.Rdata');
f1 = 'MgX'
load(data2)
df = clean(res)%>% mutate(error.sd ='error.sd=4')# %>% filter(!(MgX %in% c('MgX=1','MgX=1.5')))
data <- aggregate(as.formula(paste0('value ~ measure + error.sd +', f1, '+ method')), na.omit(df[df$error.sd=='error.sd=4',]),function(x) mean(x[!is.na(x)])) %>% filter(measure %in% c('FDR','TPR'))
data[data$measure=='FDR',]
data[data$measure=='TPR',]


col <- c('LinDA' = brewer.pal(12, "Set3")[1], 'lme'= brewer.pal(8, "Dark2")[2], 'GLMMPQL' = brewer.pal(12, "Set3")[3],
         'LDM' = brewer.pal(12, "Set3")[4], 'glmmadaptive'= brewer.pal(12, "Set3")[5], 'glmmadaptiveP' = brewer.pal(12, "Set3")[6],
         'glmmTMB' = brewer.pal(12, "Set3")[7], 'ZIBR'= brewer.pal(12, "Set3")[8], 'ZIGMM' = brewer.pal(12, "Set3")[9],
         'ZINBMM' = brewer.pal(12, "Set3")[10], 'NBMM'= brewer.pal(12, "Set3")[11], 'glmernb' = brewer.pal(12, "Set3")[12],
         'lmer' = brewer.pal(8, "Dark2")[1], 'MaAsLin2' = brewer.pal(8, "Dark2")[2])
myplot <- function(data1='repeatedB_er1_res.Rdata', 
                   data2='repeatedB_er4_res.Rdata', 
                   wd = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/res/',
                   sub.methods= c("glmernb","GLMMPQL","ZIBR","LDM","NBMM","ZIGMM","ZINBMM","glmmadaptive","glmmTMB","LinDA","MaAsLin2"),
                   para = 'MgX',level1 = c('MgX=1','MgX=2'),level2 = c('MgX=2.5','MgX=3.5'),
                   f1 = 'MgX', name = 'Repeated_Balanced_Stool', data.type='Within-subject correlation', 
                   output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/'){
  setwd(wd)
  res <- NULL
  load(data1)
  df1 <- clean(res) %>% mutate(error.sd ='error.sd=1') %>% dplyr::filter(!!as.name(para) %in% level1) %>% dplyr::filter(method %in% sub.methods)
  df1[,para] <- as.character(df1[,para])
  df1[,para] <- sub(paste0('^',unique(df1[,para])[1],'$'),'+',df1[,para])
  df1[,para] <- sub(paste0('^',unique(df1[,para])[2],'$'),'++',df1[,para])
  
  res <- NULL
  load(data2)
  df2 <- clean(res) %>% mutate(error.sd ='error.sd=4') %>% dplyr::filter(!!as.name(para) %in% level2) %>% dplyr::filter(method %in% sub.methods)
  df2[,para] <- as.character(df2[,para])
  df2[,para] <- sub(paste0('^',unique(df2[,para])[1],'$'),'+',df2[,para])
  df2[,para] <- sub(paste0('^',unique(df2[,para])[2],'$'),'++',df2[,para])
  
  res <- rbind(df1,df2)  
  res$method <- gsub('MaAsLin2TMM','MaAsLin2(TMM)',res$method)
  res$method <- gsub('MaAsLin2CSS','MaAsLin2(CSS)',res$method)
  res$method <- gsub('MaAsLin2GMPR','MaAsLin2(GMPR)',res$method)
  if(length(grep('MaAsLin2norm$',name))>0){res$method <- gsub('^MaAsLin2$','MaAsLin2(TSS)',res$method)}
  
  res.df2 <- data_summary(res %>% dplyr::filter(measure %in% c('FDR','TPR')),paste0('value ~ measure + error.sd +', f1, '+ method'))
  
  fdr = res.df2 %>% dplyr::filter(measure =='FDR') %>% dplyr::select(c(f1, 'method','error.sd','value'))
  fdr[is.na(fdr)] = 100
  fdr = fdr %>% unite('grp',c(f1,'error.sd')) %>% spread('grp', value) %>% column_to_rownames('method') %>% as.data.frame()
  fdr <- fdr[,c(colnames(fdr)[grep('error.sd=1',colnames(fdr))],colnames(fdr)[grep('error.sd=4',colnames(fdr))]),drop =F]
  fdr = fdr[,c('+_error.sd=1', '++_error.sd=1','+_error.sd=4','++_error.sd=4'),drop = F]
  
  # Lower CI: fdr.est - 1.96 * sqrt(fdr.est * (1 - fdr.est)/ num.iter))
  fdr.CI = res.df2 %>% dplyr::filter(measure =='FDR') %>% dplyr::select(c(f1, 'method','error.sd','ymin')) %>% unite('grp',c(f1,'error.sd')) %>% spread('grp', ymin)
  fdr.CI = fdr.CI[,c('method','+_error.sd=1', '++_error.sd=1','+_error.sd=4','++_error.sd=4'),drop = F]
  
  fdr.CI[is.na(fdr.CI)] = 100
  fdr.CI = fdr.CI %>% column_to_rownames('method')
  fdr.CI = fdr.CI[rownames(fdr),,drop =F]
  test <- function(v1, v2) ifelse(v1<=0.05,"***", ifelse(v2<=0.1, "**", ifelse(v2<=0.2, "*", 'x')))
  label.fdr <- mapply(test,v1 = fdr.CI, v2 = fdr)
  colnames(label.fdr) = colnames(fdr);rownames(label.fdr) = rownames(fdr)
  rank = apply(label.fdr, 1, function(x) str_count(x, "\\*") %>% sum())
  label.fdr <- as.data.frame(label.fdr)
  label.fdr$score <- 2*rank(rank)
  
  
  ## TPR
  tpr = res.df2 %>% dplyr::filter(measure =='TPR') %>% dplyr::select(c(f1, 'error.sd','method','value'))
  tpr = tpr %>% unite('grp',c(f1,'error.sd')) %>% spread('grp', value) %>% column_to_rownames('method') %>% as.data.frame()
  tpr = tpr[rownames(fdr),,drop =F]
  tpr = tpr[,c('+_error.sd=1', '++_error.sd=1','+_error.sd=4','++_error.sd=4')]
  
  tpr.rank = apply(tpr, 2, rank)
  save(fdr, tpr.rank, file = paste0('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/heatmap/',name,'.rdata'))
  
  # add score
  tpr[,'FDR rank'] = label.fdr[,ncol(label.fdr)]
  tpr[,'TPR rank'] = rank(apply(tpr.rank, 1, function(x) round(mean(x), digits = 1)))
  
  # reorder the columns
  tpr = tpr  %>% rownames_to_column('method')
  
  ord <- apply((tpr[,c(1,(ncol(tpr)-1):ncol(tpr)),drop=F] %>% column_to_rownames('method')), 1, function(x) sum(x)) %>% sort() %>% rev() %>% names()
  
  tpr = (tpr %>% column_to_rownames('method'))[ord,,drop = F] %>% rownames_to_column('method')
  tpr[,"FDR rank"] <- tpr[,"FDR rank"]/2
  fdr.dt = label.fdr[ord,,drop =F]
  
  colnames(tpr)[2:5] <- colnames(fdr.dt)[1:4] <-  c('+','++','+ ','++ ')
  colnames(tpr)[c(1,6:7)] = c('Effect size','FDR','TPR')
  
  
  
  for(i in c(2:(ncol(tpr)-2))){tpr[,i] = format(round(tpr[,i], digits = 2), nsmall = 2)}
  
  bar_chart <- function(label, width = "100%", height = "15px", fill = "forestgreen", background = NULL) {
    bar <- div(style = list(background = fill, width = width, height = height))
    chart <- div(style = list(flexGrow = 1, marginLeft = "1px", background = background), bar)
    div(style = list(display = "flex", alignItems = "center"), label, chart)
  }
  
  get_color <- function(fdr){
    if(fdr =='***'){
      col = brewer.pal(8,'Dark2')[5]
    }else if(fdr =='**'){
      col = brewer.pal(8,'Set2')[6]
    }else if(fdr =='*'){
      col = brewer.pal(8,'RdGy')[2]
    }else{
      col = 'grey'
    }
    return(col)
  }
  
  # making table
  NM <- colnames(tpr)[grep('\\+',colnames(tpr))]
  
  ## define color by fdr.dt matrix
  coldefs.list <- function(numcols){
    coldefs_list = NULL
    color_list = list()
    for (idx in 1:length(numcols)){
      fdr_col = fdr.dt[, idx] # column 8 is fdr score column, which does not match tpr table, need to be deleted; +1: fdr first column is methods
      colors <- sapply(fdr_col, function(x) get_color(x), USE.NAMES = F)
      name = numcols[idx]
      color_list[[name]] = colors
      
      cell.func <- function(value, index, name) {
        width <- value
        bar_chart(format(round(value, digits = 2), nsmall = 2), width = width*60, background = brewer.pal(12,'Set3')[9],
                  fill = color_list[[name]][index])
      }
      
      style.func <- function(value, index, name) {
        color <- color_list[[name]][index]
        list(fontWeight = 600, fontSize = 24,color = color)}
      
      coldefs <- list(reactable::colDef(style = style.func, cell=cell.func, name = NM[idx], align = 'center'))
      coldefs_list = c(coldefs_list,coldefs)
    }
    # change names
    names(coldefs_list) <- numcols
    return (list(coldefs_list, color_list))
  }
  
  numcols <- colnames(tpr)[c(2:5)]
  
  coldefs_list <- coldefs.list(numcols)[[1]]
  color_list <- coldefs.list(numcols)[[2]]
  
  minW = 50;fontsize = 24;fontweight = 600
  coldefs_list[[colnames(tpr)[1]]] = colDef(minWidth = 230, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[1]
    name <- tpr[,1][index]
    tagList(div(style = list(fontWeight = fontweight, fontSize = fontsize), name))})
  coldefs_list[[colnames(tpr)[6]]] = colDef(name = "FDR",minWidth = minW, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[6]
    name <- tpr[,6][index]
    tagList(div(style = list(fontWeight = fontweight, fontSize = fontsize), name))})
  coldefs_list[[colnames(tpr)[7]]] = colDef(name = "TPR",minWidth = minW, align = 'center', cell = function(value, index) {
    name = colnames(tpr)[7]
    name <- tpr[,7][index]
    tagList(div(style = list(fontWeight = fontweight, fontSize = fontsize), name))})
  
  tpr[,c(2:5)] <- apply(tpr[,c(2:5)], 2, function(x) as.numeric(x))
  
  html_file <- paste0(output,name,'.html')
  tb <- reactable(tpr, pagination=FALSE,resizable = FALSE, wrap = FALSE, bordered = F,
                  style = list(fontFamily = "Helvetica", fontSize = "24px",fontWeight= 100),
                  theme = reactableTheme(headerStyle = list("&:hover[aria-sort]" = list(background = "#f7f7f8"),
                                                            "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 98%)"),
                                                            borderColor = "#555")),
                  columnGroups = list(
                    colGroup(name = data.type, columns = colnames(tpr)[1]),
                    colGroup(name = 'High', columns = colnames(tpr)[2:3]),
                    colGroup(name = 'Low', columns = colnames(tpr)[4:5]),
                    colGroup(name = "Rank", columns = colnames(tpr)[6:7])
                  ),
                  columns = coldefs_list
  )
  saveWidget(widget = tb, file = html_file, selfcontained = TRUE)
  
  
  img_file <- paste0(output,'png/',name,'.png')
  webshot2::webshot(url = html_file, file = img_file, vwidth = 1100, vheight = 600, zoom =3) # change zoom to a higher number
  
}









setwd("/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/res/")
## ------- reapeated measure
## Balanced
myplot(data1='repeatedB_er1_res.Rdata', data2='repeatedB_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=1','MgX=2'),level2 = c('MgX=2.5','MgX=3.5'),
       f1 = 'MgX', name = 'Repeated_Balanced_Stool', 
       output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced
myplot(data1='repeatedU_er1_res.Rdata', data2='repeatedU_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=1','MgX=2'),level2 = c('MgX=2.5','MgX=3.5'),
       f1 = 'MgX', name = 'Repeated_UnBalanced_Stool',output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## ------- matched case
## Balanced
myplot(data1='matchedB_er1_res.Rdata', data2='matchedB_er4_res.Rdata',
       para = 'MgT',level1 = c('MgT=0.5','MgT=0.8'),level2 = c('MgT=2.5','MgT=3.5'),
       f1 = 'MgT', name = 'Matched_Balanced_Stool', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced
myplot(data1='matchedU_er1_res.Rdata', data2='matchedU_er4_res.Rdata',
       para = 'MgT',level1 =c('MgT=0.55','MgT=0.8'),level2 = c('MgT=4','MgT=5'),
       f1 = 'MgT', name = 'Matched_UnBalanced_Stool',output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## ------- longitudinal: X
## Balanced
myplot(data1='longXB_er1_res.Rdata', data2='longXB_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=2','MgX=3'),level2 = c('MgX=3','MgX=4.5'),
       sub.methods= c("glmernb","GLMMPQL","ZIBR","NBMM","ZIGMM","ZINBMM","glmmadaptive","glmmTMB","LinDA","MaAsLin2"),
       f1 = 'MgX', name = 'longX_Balanced_Stool', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced
myplot(data1='longXU_er1_res.Rdata', data2='longXU_er4_res.Rdata',
       para = 'MgX',level1 = c('MgX=2','MgX=3'),level2 = c('MgX=3','MgX=4'),
       sub.methods= c("glmernb","GLMMPQL","ZIBR","NBMM","ZIGMM","ZINBMM","glmmadaptive","glmmTMB","LinDA","MaAsLin2"),
       f1 = 'MgX', name = 'longX_UnBalanced_Stool', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## ------- longitudinal: T
## Balanced
myplot(data1='longTB_er1_res.Rdata', data2='longTB_er4_res.Rdata',
       para = 'MgT',level1 = c('MgT=0.6','MgT=0.8'),level2 = c('MgT=2.5','MgT=3'),
       sub.methods= c("glmernb","GLMMPQL","ZIBR","NBMM","ZIGMM","ZINBMM","glmmadaptive","glmmTMB","LinDA","MaAsLin2"),
       f1 = 'MgT', name = 'longT_Balanced_Stool', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced
myplot(data1='longTU_er1_res.Rdata', data2='longTU_er4_res.Rdata',
       para = 'MgT',level1 = c('MgT=0.6','MgT=0.8'),level2 = c('MgT=2.5','MgT=3'),
       sub.methods= c("glmernb","GLMMPQL","ZIBR","NBMM","ZIGMM","ZINBMM","glmmadaptive","glmmTMB","LinDA","MaAsLin2"),
       f1 = 'MgT', name = 'longT_UnBalanced_Stool', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')



## Balanced:sample sample
myplot(data1='smallsampleB_er1_res.Rdata', data2='smallsampleB_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=4','MgX=6'),level2 = c('MgX=8','MgX=10'),
       f1 = 'MgX', name = 'smallsample_Balanced_Stool', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced:sample sample
myplot(data1='smallsampleU_er1_res.Rdata', data2='smallsampleU_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=4','MgX=6'),level2 = c('MgX=8','MgX=10'),
       f1 = 'MgX', name = 'smallsample_UnBalanced_Stool',output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Balanced:sample otu
myplot(data1='smallotuB_er1_res.Rdata', data2='smallotuB_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=1','MgX=2'),level2 = c('MgX=1.8','MgX=2.3'),
       f1 = 'MgX', name = 'smallotu_Balanced_Stool', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced:sample otu
myplot(data1='smallotuU_er1_res.Rdata', data2='smallotuU_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=1','MgX=2'),level2 = c('MgX=2.5','MgX=3.5'),
       f1 = 'MgX', name = 'smallotu_UnBalanced_Stool', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')







### vaginal
## ------- reapeated measure
## Balanced
myplot(data1='repeatedB_v_er1_res.Rdata', data2='repeatedB_v_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=2','MgX=3'),level2 = c('MgX=2.5','MgX=3.5'),
       f1 = 'MgX',name = 'Repeated_Balanced_Vaginal', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced
myplot(data1='repeatedU_v_er1_res.Rdata', data2='repeatedU_v_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=2','MgX=3'),level2 = c('MgX=2.5','MgX=3.5'),
       f1 = 'MgX',name = 'Repeated_UnBalanced_Vaginal', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## ------- matched case
## Balanced
myplot(data1='matchedB_v_er1_res.Rdata', data2='matchedB_v_er4_res.Rdata',
       para = 'MgT',level1 = c('MgT=0.9','MgT=1.8'),level2 = c('MgT=2.5','MgT=4.5'),
       f1 = 'MgT',name = 'Matched_Balanced_Vaginal', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced
myplot(data1='matchedU_v_er1_res.Rdata', data2='matchedU_v_er4_res.Rdata',
       para = 'MgT',level1 =c('MgT=1','MgT=1.8'),level2 = c('MgT=4.5','MgT=6.5'),
       f1 = 'MgT',name = 'Matched_UnBalanced_Vaginal', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## ------- longitudinal: X
## Balanced
myplot(data1='longXB_v_er1_res.Rdata', data2='longXB_v_er4_res.Rdata', 
       para = 'MgX',level1 = c('MgX=5','MgX=9'),level2 = c('MgX=4','MgX=8'),
       sub.methods= c("glmernb","GLMMPQL","ZIBR","NBMM","ZIGMM","ZINBMM","glmmadaptive","glmmTMB","LinDA","MaAsLin2"),
       f1 = 'MgX',name = 'longX_Balanced_Vaginal', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced
myplot(data1='longXU_v_er1_res.Rdata', data2='longXU_v_er4_res.Rdata',
       para = 'MgX',level1 = c('MgX=3.5','MgX=5.5'),level2 = c('MgX=4','MgX=5'),
       sub.methods= c("glmernb","GLMMPQL","ZIBR","NBMM","ZIGMM","ZINBMM","glmmadaptive","glmmTMB","LinDA","MaAsLin2"),
       f1 = 'MgX',name = 'longX_UnBalanced_Vaginal',output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')


## ------- longitudinal: T
## Balanced
myplot(data1='longTB_v_er1_res.Rdata', data2='longTB_v_er4_res.Rdata',
       para = 'MgT',level1 = c('MgT=0.8','MgT=1.2'),level2 = c('MgT=3','MgT=4'),
       sub.methods= c("glmernb","GLMMPQL","ZIBR","NBMM","ZIGMM","ZINBMM","glmmadaptive","glmmTMB","LinDA","MaAsLin2"),
       f1 = 'MgT',name = 'longT_Balanced_Vaginal',output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

## Unbalanced
myplot(data1='longTU_v_er1_res.Rdata', data2='longTU_v_er4_res.Rdata',
       para = 'MgT',level1 = c('MgT=0.8','MgT=1.2'),level2 = c('MgT=2.5','MgT=3.5'),
       sub.methods= c("glmernb","GLMMPQL","ZIBR","NBMM","ZIGMM","ZINBMM","glmmadaptive","glmmTMB","LinDA","MaAsLin2"),
       f1 = 'MgT',name = 'longT_UnBalanced_Vaginal', output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')










## NULL
load('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/res/U_res.Rdata')
df1 = clean(res) %>% filter(measure=='FP') %>% mutate(balanced = 'Vaginal')
load('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/B_res.Rdata')
df2 = clean(res) %>% filter(measure=='FP')  %>% mutate(balanced = 'Stool') 
df <- rbind(df1, df2) %>% droplevels();head(df)
df[is.na(df$value),'value'] <- 0
df$ct <- 0
df[df$value >0,'ct'] <- 1
# df = data_summary(df, ct ~ measure + error.sd + method + balanced) %>% dplyr::select(measure, error.sd, method, balanced, ct)
res <- aggregate(ct ~ error.sd + method + balanced, data = df, function(x) sum(x)/1000)
head(res)
res <- res %>% unite('grp',c('balanced','error.sd')) %>% spread('grp', ct) %>% column_to_rownames('method') %>% as.data.frame()
# save(res,file = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/heatmap/null.Rdata')
load('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/heatmap/null.Rdata')
ratings <- apply(res, 2, function(x) ifelse(x <=0.05,"3", ifelse(x<=0.1, "2", ifelse(x<=0.2, "1", '0'))))
ratings1 <- apply(ratings, 2, function(x) as.numeric(x))
rownames(ratings1) <- rownames(ratings)
ord <-  names(sort(-rowSums(ratings1)))
ratings1 <- ratings1[ord,,drop =F]

star <- function(rating, max_rating = 3) {
  star_icon <- function(max_rating, empty = FALSE) {
    col_func <- ifelse(max_rating =='3',brewer.pal(8, 'Set1')[3], ifelse(max_rating =='2', brewer.pal(8, 'Set2')[6], brewer.pal(8, 'Set1')[1]))
    tagAppendAttributes(shiny::icon("star"),
                        style = paste("color:", if (empty) "#edf0f2" else col_func),
                        "aria-hidden" = "true")
  }
  rounded_rating <- rating  # always round up
  stars <- lapply(seq_len(max_rating), function(i) {
    if (i <= rounded_rating) star_icon(rating) else star_icon(rating,empty = TRUE)
  })
  label <- sprintf("%s out of %s stars", rating, max_rating)
  div(title = label, role = "img", stars)
}


colnames(ratings1) <- c("error.sd=1","error.sd=4","error.sd=1 ","error.sd=4 ")




output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/'
html_file <- paste0(output,'Null.html')
colnames(ratings1) <- c('High','Low','High ','Low ')
tb <- reactable(ratings1, pagination=FALSE,resizable = FALSE, wrap = FALSE, bordered = F,
                style = list(fontFamily = "Helvetica", fontSize = "24px",fontWeight= 600), 
                theme = reactableTheme(headerStyle = list("&:hover[aria-sort]" = list(background = "#f7f7f8"),
                                                          "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 98%)"),
                                                          borderColor = "#555")),
                columns = list(
                  `High` = colDef(cell = function(value) star(value)),
                  `Low` = colDef(cell = function(value) star(value)),
                  `High `  = colDef(cell = function(value) star(value)),
                  `Low ` = colDef(cell = function(value) star(value))),
                columnGroups = list(
                  colGroup(name = 'Balanced', columns = colnames(ratings1)[1:2]),
                  colGroup(name = 'Unbalanced', columns = colnames(ratings1)[3:4])
                )
)
saveWidget(widget = tb, file = html_file, selfcontained = TRUE)
img_file <- paste0(output,'Null.png')
webshot(url = html_file, file = img_file, vwidth = 1000, vheight = 600, zoom =3)



df <- rbind(df1, df2)
df[is.na(df$value),'value'] <- 0
df <- aggregate(value ~ error.sd + method + balanced, data = df, function(x) (mean(x)/500))
df$error.sd <- gsub('error.sd=1','High',df$error.sd)
df$error.sd <- gsub('error.sd=4','Low',df$error.sd)
p2 <- ggplot(df, aes(x = method, y = value*100, fill = method)) +
  geom_bar(stat = 'identity') + facet_grid(error.sd ~ balanced, scales = 'free') + theme_bw() +
  scale_fill_manual(values = col) +
  theme(text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, color = 'black'),
        axis.text.y = element_text(color = 'black'),
        legend.text = element_text(size = 16, color = "black"),
        axis.title =  element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(),
        legend.position = 'none',
        panel.grid.minor = element_blank()) + labs(y = '% of differential taxa', x= '') 
p2
ggsave('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/FigureS2.pdf', width =7, height = 7, dpi = 100)


do.call(file.remove, list(list.files("/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/", pattern = '*html',full.names = TRUE)))









########################## Time ########################## 
files = list.files('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/res/', pattern = 'Rdata$')
files = files[grep('long|matched|repeat',files)]
setwd('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/res/')
# files <- files[grep('repeated',files)]
DF <- dff <- NULL
for(i in 1:length(files)){
  df <- NULL
  load(files[i])
  df <- clean(res) %>% dplyr::filter(measure == 'time') %>% dplyr::select(method, iter, value) 
  if(sum(is.na(df$value))!=nrow(df)){cat(files[i],'\n')}
  dff <- rbind(dff, df)
  df <- df %>% group_by(method) %>% summarise(mean = mean(value))
  DF <- rbind(DF, df)
}

DF1 <- DF %>% na.omit()
dff1 <- dff %>% na.omit()
ggplot(dff1, aes(x = reorder(method, value), y = value/60)) + geom_boxplot() +
  theme_bw() + 
  guides(fill=  F) + 
  labs(x = '', y = 'time(mins)', fill = '') + 
  theme(axis.text.x = element_text(color = 'black', size = 14, angle = 90, hjust = 1),
        axis.text = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none')
ggsave('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/time.pdf', width = 6, height = 5)

DF1 <- aggregate(mean ~ method, DF1, function(x) mean(x))
time <- DF1[order(-DF1$mean),]
rownames(time) <- NULL
time1 <- as.vector(t(time %>% column_to_rownames('method')))
names(time1) <- as.vector(time$method)
time1/60
# save(time1, file = 'long_time.rdata')

########################## Heatmap ##########################
files <- list.files('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/heatmap/')
files <- files[grep('Vaginal|Stool', files)]
files <- c("Repeated_Balanced_Vaginal.rdata",files[!(files %in% "Repeated_Balanced_Vaginal.rdata")]) ## LDM does not apply to general longitudinal data
load(paste0('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/heatmap/Repeated_Balanced_Vaginal.rdata'))
names <- names(rowMeans(fdr))
FDR <- TPR <- nm <- NULL
for(file in files){
  load(paste0('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/heatmap/',file))
  fdr1 <- rowMeans(fdr)
  fdr1 <- fdr1[names]
  tpr1 <- rowMeans(tpr.rank)
  tpr1 <- tpr1[names]
  nm <- c(nm, gsub('.rdata','',file))
  cat(names(fdr1), '\n')
  cat(names(tpr1), '\n')
  FDR <- rbind(FDR, fdr1)
  TPR <- rbind(TPR,tpr1)
}

rownames(FDR) <- rownames(TPR) <- nm

FDR <- apply(FDR, 2, function(x) ifelse(x < 0.05, 'Good', ifelse(x>0.2,'Poor','Intermediate'))) %>% as.data.frame()
TPR <- apply(TPR, 2, function(x) ifelse(x > (11*0.7), 'Good', ifelse(x<(11* 0.3),'Poor','Intermediate'))) %>% as.data.frame()
time1 <- apply(time %>% column_to_rownames('method'), 2, function(x) ifelse(x < 60, 'Good', ifelse(x>600,'Poor','Intermediate'))) %>% t()
time1 <- time1[,colnames(FDR), drop =F]
rownames(time1) = 'Speed'
time1 = cbind(time1, class = 'yOther')
colnames(FDR)==colnames(TPR)
colnames(FDR)==colnames(time1)[1:(ncol(time1)-1)]

load('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/heatmap/rho.Rdata')
rownames(df_mean) <- NULL
rhos <- as.data.frame(df_mean) %>% column_to_rownames('methods') %>% t()
rhos[rhos > .8] <- 'Good'
rhos <- cbind(rhos, class = 'yOther')
rhos <- rhos[,colnames(time1), drop =F]
rownames(rhos) <- 'Stability'

load('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/heatmap/null.Rdata')
null = res %>% mutate(Stool = (`Stool_error.sd=1`+`Stool_error.sd=4`)/2,
                      Vaginal = (`Vaginal_error.sd=1`+`Vaginal_error.sd=4`)/2) %>% dplyr::select(Stool, Vaginal) %>% t()
null = apply(null, 2, function(x) ifelse(x < 0.05, 'Good', ifelse(x>0.2,'Poor','Intermediate'))) %>% as.data.frame()

null <- cbind(null, class = c('1Null','1Null'))
null <- null[,colnames(time1), drop =F]

d1 <- rbind(as.data.frame(FDR) %>% mutate(class = 'FDR'), as.data.frame(TPR) %>% mutate(class = 'TPR')) %>% rbind(time1) %>% rbind(rhos) %>% rbind(null)
d1 <- d1 %>% rownames_to_column('type')
rownames(d1) <- d1$type
d1$type <- gsub('_Balanced_',' ',d1$type)
d1$type <- gsub('_UnBalanced_',' ',d1$type)
d1$type <- gsub('1','',d1$type)

# d1[is.na(d1)] <- 'Good'
d1$type[grep('smallsample',rownames(d1))] <- paste0('Small sample size')
d1$type[grep('smallotu',rownames(d1))] <- paste0('Small taxa number')
d1$comp = 'Other'
d1$comp[grep('_Balanced_',rownames(d1))] <- 'Balanced'
d1$comp[grep('_UnBalanced_',rownames(d1))] <- 'Unbalanced'

idx.fdr.stool1 <- intersect(intersect(which(d1$class=='FDR'),grep('Stool',d1$type)), grep('Balanced', d1$comp))
idx.fdr.vaginal1 <- intersect(intersect(which(d1$class=='FDR'),grep('Vaginal',d1$type)), grep('Balanced', d1$comp))

idx.fdr.stool2 <- intersect(intersect(which(d1$class=='FDR'),grep('Stool',d1$type)), grep('Unbalanced', d1$comp))
idx.fdr.vaginal2 <- intersect(intersect(which(d1$class=='FDR'),grep('Vaginal',d1$type)), grep('Unbalanced', d1$comp))

idx.fdr.smallotu <- which(d1$type=='Small taxa' & d1$class =='FDR')
idx.fdr.smallsample <- which(d1$type=='Small sample'  & d1$class =='FDR')

idx.tpr.smallotu <- which(d1$type=='Small taxa' & d1$class =='TPR')
idx.tpr.smallsample <- which(d1$type=='Small sample'  & d1$class =='TPR')

idx.tpr.stool1 <- intersect(intersect(which(d1$class=='TPR'),grep('Stool',d1$type)), grep('Balanced', d1$comp))
idx.tpr.vaginal1 <- intersect(intersect(which(d1$class=='TPR'),grep('Vaginal',d1$type)), grep('Balanced', d1$comp))

idx.tpr.stool2 <- intersect(intersect(which(d1$class=='TPR'),grep('Stool',d1$type)), grep('Unbalanced', d1$comp))
idx.tpr.vaginal2 <- intersect(intersect(which(d1$class=='TPR'),grep('Vaginal',d1$type)), grep('Unbalanced', d1$comp))

idx.other <- which(d1$class=='Other')
idx.null <- grep('Null',d1$class)

length(idx.fdr.stool1) + length(idx.fdr.stool1) + length(idx.tpr.stool1) + length(idx.tpr.stool1) + 
  length(idx.fdr.stool2) + length(idx.fdr.stool2) + length(idx.tpr.stool2) + length(idx.tpr.stool2) +
  length(idx.other) + length(idx.tpr.smallotu) + length(idx.tpr.smallotu)+ length(idx.fdr.smallotu) + length(idx.fdr.smallotu)
dim(d1)

d1 <- d1[c(idx.null,idx.fdr.stool1,idx.fdr.stool2, idx.fdr.vaginal1,idx.fdr.vaginal2, idx.fdr.smallotu, idx.fdr.smallsample, idx.tpr.stool1,idx.tpr.stool2, idx.tpr.vaginal1,idx.tpr.vaginal2, idx.tpr.smallotu, idx.tpr.smallsample, idx.other),,drop =F]
head(d1)

d1$type1 <- gsub(' Stool| Vaginal','',d1$type)
d1$type[grep('Stool',d1$type)] <- 'Stool'
d1$type[grep('Vaginal',d1$type)] <- 'Vaginal'
d1$type[grep('Small',d1$type)] <- 'Stool'

d1$type1[grep('longX', d1$type1)] = 'Longitudinal(X)'
d1$type1[grep('longT', d1$type1)] = 'Longitudinal(T)'
d1$type1[grep('Matched', d1$type1)] = 'Matched-pair'
d1$type1[grep('Repeated', d1$type1)] = 'Replicate sampling'
d1$type1[grep('Stool|Vaginal',d1$type1)] <- 'Null'

cols = c("Poor"=brewer.pal(8,'Set1')[1],"Good"=brewer.pal(8,'Set1')[3],"Intermediate"=brewer.pal(8,'Set2')[6])
rowsplit = d1$class
rowann = d1[,c('type','comp','type1'),drop =F]
ann_name <- d1$type1
library(ComplexHeatmap)
row_ha = rowAnnotation(df = rowann[,c(1,2)], foo = anno_text(ann_name, location =0, just = "left"),simple_anno_size = unit(0.3, "cm"), 
                       col = list(type = c("Speed" = brewer.pal(8,'Dark2')[3],"Stability" = brewer.pal(8,'Dark2')[3], 
                                           "Stool" = brewer.pal(8,'Dark2')[2],'Vaginal' = brewer.pal(8,'Dark2')[1]),
                                  comp = c("Unbalanced" = brewer.pal(8,'Pastel2')[1], "Balanced" = brewer.pal(8,'Pastel2')[2],"Other" = brewer.pal(8,'Pastel2')[3])))

pdf('/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/heatmap_f.pdf', width =6, height = 9)
idx1 = colnames(d1)[grep('MaAsLin2|LinDA|LDM',colnames(d1))]
idx2 = colnames(d1)[-grep('LDM|MaAsLin2|LinDA',colnames(d1))]
d1 = d1[,c(idx1,idx2)]
ht = Heatmap(as.matrix(d1 %>% dplyr::select(-c('type','type1','class','comp'))),
             col = cols,
             rect_gp = gpar(col = "white", lwd = 0.5), 
             show_row_names =F,cluster_rows = FALSE, cluster_columns = F,
             border = TRUE,row_split = rowsplit,
             right_annotation = row_ha)
draw(ht, padding = unit(c(1, 1, 1, 1), "cm"))
dev.off()







## Add 05102022: Compare MaAsLin2 different normalization methods
myplot(data1='repeatedU_er1_res.Rdata', data2='repeatedU_er4_res.Rdata', 
       wd = "/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Data/res05032022/",
       sub.methods= c('MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR','MaAsLin2','LinDA'),
       para = 'MgX',level1 = c('MgX=1','MgX=2'),level2 = c('MgX=2.5','MgX=3.5'),
       f1 = 'MgX', name = 'Repeated_UnBalanced_Stool_MaAsLin2norm',
       output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

myplot(data1='matchedU_er1_res.Rdata', data2='matchedU_er4_res.Rdata',
       wd = "/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Data/res05032022/",
       sub.methods= c('MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR','MaAsLin2','LinDA'),
       para = 'MgT',level1 =c('MgT=0.55','MgT=0.8'),level2 = c('MgT=4','MgT=5'),
       f1 = 'MgT', name = 'Matched_UnBalanced_Stool_MaAsLin2norm', 
       output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

myplot(data1='longXU_er1_res.Rdata', data2='longXU_er4_res.Rdata',
       wd = "/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Data/res05032022/",
       sub.methods= c('MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR','MaAsLin2','LinDA'),
       para = 'MgX',level1 = c('MgX=2','MgX=3'),level2 = c('MgX=3','MgX=4'),
       f1 = 'MgX', name = 'longX_UnBalanced_Stool_MaAsLin2norm', 
       output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

myplot(data1='longTU_er1_res.Rdata', data2='longTU_er4_res.Rdata',
       wd = "/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Data/res05032022/",
       sub.methods= c('MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR','MaAsLin2','LinDA'),
       para = 'MgT',level1 = c('MgT=0.6','MgT=0.8'),level2 = c('MgT=2.5','MgT=3'),
       f1 = 'MgT', name = 'longT_UnBalanced_Stool_MaAsLin2norm',
       output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

myplot(data1='repeatedU_v_er1_res.Rdata', data2='repeatedU_v_er4_res.Rdata', 
       wd = "/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Data/res05032022/",
       sub.methods= c('MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR','MaAsLin2','LinDA'),
       para = 'MgX',level1 = c('MgX=2','MgX=3'),level2 = c('MgX=2.5','MgX=3.5'),
       f1 = 'MgX',name = 'Repeated_UnBalanced_Vaginal_MaAsLin2norm', 
       output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

myplot(data1='matchedU_v_er1_res.Rdata', data2='matchedU_v_er4_res.Rdata',
       wd = "/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Data/res05032022/",
       sub.methods= c('MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR','MaAsLin2','LinDA'),
       para = 'MgT',level1 =c('MgT=1','MgT=1.8'),level2 = c('MgT=4.5','MgT=6.5'),
       f1 = 'MgT',name = 'Matched_UnBalanced_Vaginal_MaAsLin2norm', 
       output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

myplot(data1='longXU_v_er1_res.Rdata', data2='longXU_v_er4_res.Rdata',
       wd = "/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Data/res05032022/",
       sub.methods= c('MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR','MaAsLin2','LinDA'),
       para = 'MgX',level1 = c('MgX=3.5','MgX=5.5'),level2 = c('MgX=4','MgX=5'),
       f1 = 'MgX',name = 'longX_UnBalanced_Vaginal_MaAsLin2norm',
       output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')

myplot(data1='longTU_v_er1_res.Rdata', data2='longTU_v_er4_res.Rdata',
       wd = "/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Data/res05032022/",
       sub.methods= c('MaAsLin2CSS','MaAsLin2TMM','MaAsLin2GMPR','MaAsLin2','LinDA'),
       para = 'MgT',level1 = c('MgT=0.8','MgT=1.2'),level2 = c('MgT=2.5','MgT=3.5'),
       f1 = 'MgT',name = 'longT_UnBalanced_Vaginal_MaAsLin2norm', 
       output = '/Users/m216453/Library/Mobile Documents/com~apple~CloudDocs/Documents/Mayo_Research/2021_03_29_CorrelatedData/Result/')
