## ======= TPR analysis code======
## ---- combo ----
setwd('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/')
load("combo.bact.fungi.gen.ffq.adj.RData") #"bact.gen.abund" "fungi.gen.abund" "ffq.adj" "bact.ann" "fungi.ann"
dim(ffq.adj);dim(fungi.gen.abund);dim(bact.gen.abund)
rownames(ffq.adj) <- paste0('S',rownames(ffq.adj))
colnames(fungi.gen.abund) <- paste0('S',colnames(fungi.gen.abund))
colnames(bact.gen.abund) <- paste0('S',colnames(bact.gen.abund))
idx <- rowSums(bact.gen.abund==0) > (ncol(bact.gen.abund) * 0)
comm <- bact.gen.abund[idx,]
meta <- ffq.adj
# association test
set.seed(123)
pValues <- NULL
for(i in colnames(meta)){
  cat('.')
  meta.tmp <- meta[,i, drop =F]
  fit <- GUniFrac::ZicoSeq(meta.dat = meta.tmp, feature.dat = comm, feature.dat.type = 'count',
                           prev.filter =0.1, grp.name = i, return.feature.dat = T)
  pValues <- cbind(pValues, fit$p.raw)
}
colnames(pValues) <- colnames(meta)
P.df <- reshape2::melt(pValues)
pvals <- P.df[,3]
covariate <- P.df[,c(1,2),drop =F];colnames(covariate) <- c('gene','outcome')
covariate$gene <- as.factor(covariate$gene)
covariate$outcome <- as.factor(covariate$outcome)
dat <- list(pValues = pValues, pvals = pvals, covariate =covariate)
save(dat, file = '/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/dat/P.combo.ZicoSeq.RData')

# FDR adjustment
load('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/dat/P.combo.ZicoSeq.RData')
source('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/code/func.R')
methods_funs <- list('AdaptiveGBHo'='run_AdaptiveGBHo',
                     'AdaptiveGBHg'='run_AdaptiveGBHg',
                     'AdaptiveGBHo_storey'='run_AdaptiveGBHo_st',
                     'AdaptiveGBHg_storey'='run_AdaptiveGBHg_st',
                     'BH'='run_BH',
                     'ST'='run_ST',
                     'stratBHo'='run_stratBHo',
                     'stratSTo'='run_stratSTo',
                     'stratBHg'='run_stratBHg',
                     'stratSTg'='run_stratSTg',
                     'StoreyNew0'='Storey_New_0',
                     'StoreyNew05'='Storey_New_0.05',
                     'StoreyNew1'='Storey_New_0.1',
                     'StoreyNew2'='Storey_New_0.2',
                     'StoreyNew4'='Storey_New_0.4',
                     'StoreyGeo1'='Storey_Geo_0.1',
                     'StoreyAri1'='Storey_Ari_0.1',
                     'StoreyNew1N'='Storey_New_0.1_norm',
                     'StoreyNew1S'='Storey_New_0.1_S',
                     'tdgbh3F'='run_tdgbh3F',
                     'tdgbh3T'='run_tdgbh3T',
                     'swfdrBH'='run_swfdr.BH',
                     'swfdr'='run_swfdr',
                     'CAMT'='run_CAMT',
                     'FDRregT'='run_FDRregT',
                     'FDRregE'='run_FDRregE',
                     'adaptMT'='run_adaptMT')
# tdgbh3F tdgbh3T swfdr CAMT FDRregT FDRregE adaptMT
for(m in c('tdgbh3F','tdgbh3T','swfdr','CAMT','FDRregT','FDRregE','adaptMT')){
  # for(m in c('stratBHg','stratSTg','BH','ST','stratBHo','stratSTo','AdaptiveGBHo','AdaptiveGBHg','StoreyNew1','AdaptiveGBHo_storey','AdaptiveGBHg_storey')){
  cat(m,'\n')
  tryCatch({
    res <- time <- NULL
    wrapper <- match.fun(methods_funs[[m]])
    out <- wrapper(dat)
    res <- sum(out$fdr<=0.05)
    save(out, res, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/res/P.combo.ZicoSeq_',m,'.Rdata'))
  })
}


## ---- adenoma ----
setwd('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/')
load("adenoma.gen.mtb.RData") #"bact.gen.abund"    "mtb.subpath.abund" "meta.dat"          "mtb.ann"           "bact.ann"
meta <- t(mtb.subpath.abund)
colnames(meta)[1] <- 'pathway1'
comm <- bact.gen.abund
idx <- rowSums(comm==0) > (ncol(bact.gen.abund) * 0)
comm <- comm[idx,];dim(comm)

set.seed(123)
pValues <- NULL
for(i in colnames(meta)){
  cat('[.]')
  meta.tmp <- as.data.frame(meta[,i, drop =F])
  suppressMessages(fit <- GUniFrac::ZicoSeq(meta.dat = meta.tmp, feature.dat = comm, feature.dat.type = 'count',
                                            prev.filter =0.1, grp.name = i, return.feature.dat = T))
  pValues <- cbind(pValues, fit$p.raw)
}
colnames(pValues) <- colnames(meta)
P.df <- reshape2::melt(pValues)
pvals <- P.df[,3]
covariate <- P.df[,c(1,2),drop =F];colnames(covariate) <- c('gene','outcome')
covariate$gene <- as.factor(covariate$gene)
covariate$outcome <- as.factor(covariate$outcome)
dat <- list(pValues = pValues, pvals = pvals, covariate =covariate)
save(dat, file = paste0("/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/dat/P.adenoma.ZicoSeq.RData"))

# FDR adjustment
load("/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/dat/P.adenoma.ZicoSeq.RData")
for(m in c('tdgbh3F','tdgbh3T','swfdr','CAMT','FDRregT','FDRregE','adaptMT')){#,'stratBHg','stratSTg','BH','ST','stratBHo','stratSTo','AdaptiveGBHo','AdaptiveGBHg','StoreyNew1','AdaptiveGBHo_storey','AdaptiveGBHg_storey'
  cat(m,'\n')
  tryCatch({
    res <- time <- NULL
    wrapper <- match.fun(methods_funs[[m]])
    out <- wrapper(dat)
    res <- sum(out$fdr<=0.05)
    save(out, res, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/res/P.adenoma.ZicoSeq_',m,'.Rdata'))
  })
}

## ---- autism (submit to cluster by calculateP_noshuffle.sh)----
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  file = args[1]
  k = args[2]
}

library(Seurat)
library(plyr)
library(dplyr)
library(tibble)
library(reshape2)
library(edgeR, quietly = T)
library(DESeq2)

norm ='GMPR'
source('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/code/func.R')
setwd('/research/bsi/projects/staff_analysis/m216453/scRNA_pool/2group/')
file = 'Autism.RData'
load(file)
folder <- paste0("/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/noshuffle/P/",gsub('.RData','',file))
if(!dir.exists(folder)){dir.create(folder)}
clusters <- colnames(seurat.obj@meta.data)[grep('^k',colnames(seurat.obj@meta.data))]
celltype <- which(apply(cbind(table(seurat.obj$disease,seurat.obj@meta.data[,clusters[length(clusters)]])), 2, function(x) sum(x)>100))
idx <- names(which(rowSums(seurat.obj@assays$RNA@counts)>0)) # Filter: genes with counts > 0 in all tested samples
expr.raw.ct <- seurat.obj@assays[['RNA']]@counts
meta.cells <- seurat.obj@meta.data
meta <- unique(meta.cells[,c('subject',"disease")])
rownames(meta) <- NULL
meta <- meta %>% column_to_rownames('subject')
cell.column <- clusters[length(clusters)]
# subset the cell type
meta.cells.sub <- meta.cells[meta.cells[,cell.column] %in% unname(k), , drop =F]
expr.raw.ct.sub <- expr.raw.ct[,rownames(meta.cells.sub), drop =F]

# genes with expression in each cell type > 100
x1 <- rowSums(expr.raw.ct.sub)
Y <- as(expr.raw.ct.sub, "lgCMatrix") #should be more efficient than X != 0
Y@x[] <- TRUE #set all values to TRUE
ind <- x1 > 0

m1 <- meta.cells.sub[meta.cells.sub$disease ==unique(meta.cells.sub$disease)[1],]
expr.raw.ct.sub1 <- expr.raw.ct.sub[,rownames(m1)]
Y1 <- as(expr.raw.ct.sub1, "lgCMatrix") #should be more efficient than X != 0
Y1@x[] <- TRUE #set all values to TRUE
ind1 <- x1 > 10

m2 <- meta.cells.sub[meta.cells.sub$disease ==unique(meta.cells.sub$disease)[2],]
expr.raw.ct.sub2 <- expr.raw.ct.sub[,rownames(m2)]
Y2 <- as(expr.raw.ct.sub2, "lgCMatrix") #should be more efficient than X != 0
Y2@x[] <- TRUE #set all values to TRUE
ind2 <- x1 > 10

feature.names <- unique(c(names(which(ind1)),names(which(ind2))))

expr.raw.ct.sel <- expr.raw.ct.sub[feature.names, ]
summary(rowSums(expr.raw.ct.sel));dim(expr.raw.ct.sel)
sam.ids <- meta.cells.sub$subject

# summary the cells
for(feature.name in feature.names){
  measurement.list1 <- tapply(expr.raw.ct.sel[feature.name, ], factor(sam.ids), function (x) x)
  sum.measurement <- sapply(measurement.list1, function (x) sum(x))
}
capture.list1 <- lapply(feature.names, global.function1)
sum.measurement <- list_to_df(capture.list1)
rownames(sum.measurement) <- feature.names
sum.measurement <- sum.measurement[,colSums(sum.measurement)>0]
conditions <- meta[colnames(sum.measurement),]

keep <- rowSums(sum.measurement) > 0
y <- sum.measurement[keep,, drop =F]
if(norm =='RLE'){
  suppressMessages(dds <- DESeqDataSetFromMatrix(round(y), DataFrame(conditions), ~conditions))
  suppressMessages(dds <- estimateSizeFactors(dds))
  s <- sizeFactors(dds)
}
if(norm =='GMPR'){s <- GUniFrac::GMPR(y)}
if(norm =='TMM'){s <- calcNormFactors(y, method = 'TMM')} # I think this function is not correctly written
if(norm =='TSS'){s <- colSums(y)}
count_norm <- y/s
pvalues <- sapply(1:nrow(count_norm),function(i){
  data <- cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p <- wilcox.test(gene~conditions, data)$p.value
  return(p)
})
names(pvalues) <- rownames(count_norm)
p <- pvalues
folder.data <- paste0('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/noshuffle/data/',gsub('.RData','',file),'/')
if(!dir.exists(folder.data)){dir.create(folder.data)}
save(conditions, sum.measurement, file = paste0(folder.data,gsub('.RData','',file),'_ident',k,'.FS.Rdata'))
ident <- paste0('ident',unname(k))
save(p,ident, file = paste0(folder,"/P.",gsub('.RData','',file),'_',k,".Wilcox.RData"))



## convert independent p to dat for later FDR fit
files <- list.files('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/P/Autism/', pattern = '^P')
for(f in files){
  cat(f,'\n')
  p <- NULL
  try({
    load(paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/P/',f))
    p <- P[[1]]
    colnames(p)[2] <- 'ident1'
    for(i in 2:length(P)){
      tmp <- P[[i]]
      colnames(tmp)[2] <- paste0('ident',i)
      p <- full_join(p, tmp)
    }})
  
  if(!is.null(p)){
    pValues <- p
    P.df <- reshape2::melt(pValues)
    pValues <- as.matrix(pValues %>% tibble::column_to_rownames('gene'))
    pvals <- P.df[,3]
    covariate <- P.df[,c(1,2),drop =F];colnames(covariate) <- c('gene','outcome')
    covariate$gene <- as.factor(covariate$gene)
    covariate$outcome <- as.factor(covariate$outcome)
    dat <- list(pValues = pValues, pvals = pvals, covariate =covariate)
    save(dat, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/dat/',f))
  }
}

files <- list.files('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/P/Autism/', pattern = '^P')
setwd('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/P/Autism/')
P <- list()
for(file in files){
  p <- ident <- NULL
  load(file)
  P[[ident]] <- p
}

df <- as.data.frame(p) %>% rownames_to_column('gene')
colnames(df)[2] <- names(P)[1]
for(i in 2:length(P)){
  df.tmp <- as.data.frame(P[[i]]) %>% rownames_to_column('gene')
  colnames(df.tmp)[2] <- names(P)[i]
  df <- full_join(df, df.tmp)
}
pValues <- df
P.df <- reshape2::melt(pValues)
pValues <- as.matrix(pValues %>% tibble::column_to_rownames('gene'))
pvals <- P.df[,3]
names(pvals) <- paste0(P.df[,1], '_',P.df[,2])
covariate <- P.df[,c(1,2),drop =F];colnames(covariate) <- c('gene','outcome')
covariate$gene <- as.factor(covariate$gene)
covariate$outcome <- as.factor(covariate$outcome)
dat <- list(pValues = pValues, pvals = pvals, covariate =covariate)
save(dat, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/dat/Autism.ZicoSeq.RData'))



## fit FDR adjust methods together with run_real_noshuffle.sh 
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  f = args[1]
  m = args[2]
}
source('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/code/func.R')
load(paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/dat/',f))
methods_funs <- list('AdaptiveGBHo'='run_AdaptiveGBHo',
                     'AdaptiveGBHg'='run_AdaptiveGBHg',
                     'AdaptiveGBHo_storey'='run_AdaptiveGBHo_st',
                     'AdaptiveGBHg_storey'='run_AdaptiveGBHg_st',
                     'BH'='run_BH',
                     'ST'='run_ST',
                     'stratBHo'='run_stratBHo',
                     'stratSTo'='run_stratSTo',
                     'stratBHg'='run_stratBHg',
                     'stratSTg'='run_stratSTg',
                     'StoreyNew0'='Storey_New_0',
                     'StoreyNew05'='Storey_New_0.05',
                     'StoreyNew1'='Storey_New_0.1',
                     'StoreyNew2'='Storey_New_0.2',
                     'StoreyNew4'='Storey_New_0.4',
                     'StoreyGeo1'='Storey_Geo_0.1',
                     'StoreyAri1'='Storey_Ari_0.1',
                     'StoreyNew1N'='Storey_New_0.1_norm',
                     'StoreyNew1S'='Storey_New_0.1_S',
                     'tdgbh3F'='run_tdgbh3F',
                     'tdgbh3T'='run_tdgbh3T',
                     'swfdrBH'='run_swfdr.BH',
                     'swfdr'='run_swfdr',
                     'CAMT'='run_CAMT',
                     'FDRregT'='run_FDRregT',
                     'FDRregE'='run_FDRregE',
                     'adaptMT'='run_adaptMT')
tryCatch({
  res <- time <- NULL
  wrapper <- match.fun(methods_funs[[m]])
  out <- wrapper(dat)
  res <- sum(out$fdr<=0.05, na.rm = T)
  save(out, res, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/res/', gsub('.RData','',f), '_',m,'.Rdata'))
})


## ==== TPR plot code ======
library(ComplexHeatmap,lib.loc = '/usr/local/biotools/rpackages/R-4.2.2-2023-02-01')
library(ggplot2);library(dplyr);library(tidyr);library(reshape2);library(tibble);library(RColorBrewer)
setwd('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/res/')
col <- c('2dGBH-geo'="#FCBBA1",'2dGBH-ari'="forestgreen",'2dGBH'="#EF3B2C",'lslStoreyNew1'='green','tstStoreyNew1'='blue','lsl_New1'='yellow',
         '2dGBH-norm'="#EFEDF5",'2dGBH-P'="#B3DE69",
         '2dGBH-0'="#FFF5EB",'2dGBH-0.05'="#FEE6CE",'2dGBH-0.1'="#EF3B2C",
         '2dGBH-0.2'="#FDAE6B", '2dGBH-0.4'="#F16913", 
         "ST"="#C6DBEF","BH"="#6BAED6", 
         "stratBH_o"="#C7E9C0", "stratST_o"="#41AB5D",
         "stratBH_g"="#238443", "stratST_g"="#F7FCB9",
         "AdaptiveGBH_o"="#6A51A3", "AdaptiveGBH_g"="#9E9AC8",
         'NSC_1'="#E6AB02", 'NSC_2'="#A6761D", 
         'swfdr'="#B3E2CD", 'swfdrBH'="#FDCDAC", 'CAMT'="#DECBE4", 
         'FDRregT'="#969696", 'FDRregE'="#525252", 'adaptMT'="#FBB4AE")
setwd('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/noshuffle/res/')
line.p <- list()
for(pt in c('combo','adenoma','Autism')){
  files <- list.files(pattern = pt)
  files <- files[!(files %in% c("P.combo.ZicoSeq_AdaptiveGBHg.Rdata","P.combo.ZicoSeq_AdaptiveGBHo.Rdata",
                                "P.adenoma.ZicoSeq_AdaptiveGBHg.Rdata","P.adenoma.ZicoSeq_AdaptiveGBHo.Rdata"))]
  if(pt =='Autism'){
    files <- c("Autism_AdaptiveGBHg_storey.Rdata","Autism_AdaptiveGBHo_storey.Rdata",
               "Autism_BH.Rdata", "Autism_ST.Rdata","Autism_StoreyNew1.Rdata",
               "Autism_stratBHo.Rdata","Autism_stratSTo.Rdata")
  }
  
  P <- c()
  for(file in files){
    cat(file,': ')
    out <- NULL
    load(file)
    p <- as.vector(out$fdr)
    cat(length(p))
    cat('\n')
    P <- cbind(P, p)
  }
  rownames(P) <- paste0('g',1:nrow(P))
  colnames(P) <- gsub('Autism_|.Rdata|P.combo.ZicoSeq_|P.adenoma.ZicoSeq_|Autism.ZicoSeq_|P.adenoma.ZicoSeq_','',files)
  colnames(P)[colnames(P)=='StoreyGeo1'] <- '2dGBH-geo'
  colnames(P)[colnames(P)=='StoreyNew1'] <- '2dGBH'
  colnames(P)[colnames(P)=='stratBHo'] <- 'stratBH_o'
  colnames(P)[colnames(P)=='stratBHg'] <- 'stratBH_g'
  colnames(P)[colnames(P)=='stratSTo'] <- 'stratST_o'
  colnames(P)[colnames(P)=='stratSTg'] <- 'stratST_g'
  colnames(P)[colnames(P)=='AdaptiveGBHo_storey'] <- 'AdaptiveGBH_o'
  colnames(P)[colnames(P)=='AdaptiveGBHg_storey'] <- 'AdaptiveGBH_g'
  # colnames(P)[colnames(P)=='AdaptiveGBHo'] <- 'AdaptiveGBH_o'
  # colnames(P)[colnames(P)=='AdaptiveGBHg'] <- 'AdaptiveGBH_g'
  P <- P[, c("BH","ST","2dGBH","AdaptiveGBH_g","AdaptiveGBH_o","stratBH_o","stratST_o" )]
  sort(colSums(P<0.2, na.rm = T))
  sort(colSums(P<0.15, na.rm = T))
  sort(colSums(P<0.1, na.rm = T))
  sort(colSums(P<0.05, na.rm = T))
  df <- reshape2::melt(cbind.data.frame(`20%`=colSums(P<0.2, na.rm = T), `15%`=colSums(P<0.15, na.rm = T), `10%`=colSums(P<0.1, na.rm = T),`5%`=colSums(P<0.05, na.rm = T)) %>% rownames_to_column('method'))
  df <- within(df, variable <- factor(variable, levels = c('5%','10%','15%','20%')))
  df <- within(df, method <- factor(method, levels = c("BH","ST","2dGBH","AdaptiveGBH_g","AdaptiveGBH_o","stratBH_o","stratST_o" )))
  line.p[[pt]] <- ggplot(df, aes(x = variable, y = value)) + 
    geom_line(aes(group = method,color = method)) + 
    geom_point(aes(color = method), size = 2, shape = 1) + 
    scale_color_manual(values = col) + 
    theme_bw()+
    theme(axis.ticks.x = element_blank(),
          axis.text = element_text(color = 'black', size = 14),
          axis.title = element_text(color = 'black', size = 14),
          legend.text = element_text(color = 'black', size = 14),
          panel.grid = element_blank()) + 
    scale_color_manual(values = col) + 
    labs(color = '', x='FDR', y = '# of findings')
  # ggsave(paste0('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/result/(Figure4C)',pt,'_TPR_line_rev.pdf'), width = 6, height = 5)
  
  res00 <- P
  res00[is.na(res00)] <- 1
  sum(apply(res00, 1, function(x) sum(x<=0.05)) >0)/nrow(res00)
  sum(apply(res00, 1, function(x) sum(x<=0.1)) >0)/nrow(res00)
  sum(apply(res00, 1, function(x) sum(x<=0.15)) >0)/nrow(res00)
  sum(apply(res00, 1, function(x) sum(x<=0.2)) >0)/nrow(res00)
  res00.lt.05 <- apply(res00, 2, function(x) names(x[x<=0.05]))
  res00.lt.1 <- apply(res00, 2, function(x) names(x[x<=0.1]))
  res00.lt.15 <- apply(res00, 2, function(x) names(x[x<=0.15]))
  res00.lt.2 <- apply(res00, 2, function(x) names(x[x<=0.2]))
  pdf(paste0('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/result/(Figure4D)',pt,'_UpSetPlot_rev.pdf'), width = 6, height = 5)
  res00.lt.2.m1 = make_comb_mat(res00.lt.05)
  UpSet(res00.lt.2.m1)

  res00.lt.2.m2 = make_comb_mat(res00.lt.1)
  UpSet(res00.lt.2.m2)

  res00.lt.2.m3 = make_comb_mat(res00.lt.15)
  UpSet(res00.lt.2.m3)

  res00.lt.2.m4 = make_comb_mat(res00.lt.2)
  UpSet(res00.lt.2.m4)
  dev.off()
}
library(ggpubr)
ggarrange(line.p$combo, line.p$adenoma, line.p$Autism, nrow = 1, common.legend = T)
ggsave(paste0('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/result/(Figure4C)TPR_line_rev.pdf'), width = 10, height = 5)


## ==== FDR analysis code ======
## ---- combo + adenoma (together with calculateP_shuffle_microbiome.sh)-----
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  k = args[1]
}
setwd('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/')
load("combo.bact.fungi.gen.ffq.adj.RData") #"bact.gen.abund" "fungi.gen.abund" "ffq.adj" "bact.ann" "fungi.ann"
dim(ffq.adj);dim(fungi.gen.abund);dim(bact.gen.abund)
rownames(ffq.adj) <- paste0('S',rownames(ffq.adj))
colnames(fungi.gen.abund) <- paste0('S',colnames(fungi.gen.abund))
colnames(bact.gen.abund) <- paste0('S',colnames(bact.gen.abund))
# prevalence filter
idx <- rowSums(bact.gen.abund==0) > (ncol(bact.gen.abund) * 0)
comm <- bact.gen.abund[idx,]
meta <- ffq.adj
# association test
set.seed(k)
pValues <- NULL
for(i in colnames(meta)){
  meta.tmp <- meta[,i, drop =F]
  meta.tmp1 <- as.data.frame(meta.tmp[sample(nrow(meta.tmp)),]) ## shuffle data
  rownames(meta.tmp1) <- rownames(meta.tmp)
  colnames(meta.tmp1) <- colnames(meta.tmp)
  fit <- GUniFrac::ZicoSeq(meta.dat = meta.tmp1, feature.dat = comm, feature.dat.type = 'count',
                           prev.filter =0.1, grp.name = i, return.feature.dat = T)
  pValues <- cbind(pValues, fit$p.raw)
}
colnames(pValues) <- colnames(meta)
P.df <- reshape2::melt(pValues)
pvals <- P.df[,3]
covariate <- P.df[,c(1,2),drop =F];colnames(covariate) <- c('gene','outcome')
covariate$gene <- as.factor(covariate$gene)
covariate$outcome <- as.factor(covariate$outcome)
dat <- list(pValues = pValues, pvals = pvals, covariate =covariate)
save(dat, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/shuffle/dat/',"P.combo_",k,".ZicoSeq.RData" ))



setwd('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/')
load("adenoma.gen.mtb.RData") #"bact.gen.abund"    "mtb.subpath.abund" "meta.dat"          "mtb.ann"           "bact.ann"
meta <- t(mtb.subpath.abund)
colnames(meta)[1] <- 'pathway1'
comm <- bact.gen.abund
idx <- rowSums(comm==0) > (ncol(bact.gen.abund) * 0)
comm <- comm[idx,];dim(comm)
set.seed(k)
pValues <- NULL
for(i in colnames(meta)){
  cat('[------',i,' ------] \n')
  meta.tmp <- as.data.frame(meta[,i, drop =F])
  meta.tmp[,i] <- unname(sample(meta.tmp[,i]))
  suppressMessages(fit <- GUniFrac::ZicoSeq(meta.dat = meta.tmp, feature.dat = comm, feature.dat.type = 'count',
                                            prev.filter =0.1, grp.name = i, return.feature.dat = T))
  pValues <- cbind(pValues, fit$p.raw)
}
colnames(pValues) <- colnames(meta)

P.df <- reshape2::melt(pValues)
pvals <- P.df[,3]
covariate <- P.df[,c(1,2),drop =F];colnames(covariate) <- c('gene','outcome')
covariate$gene <- as.factor(covariate$gene)
covariate$outcome <- as.factor(covariate$outcome)
dat <- list(pValues = pValues, pvals = pvals, covariate =covariate)
save(dat, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/shuffle/dat/',"P.adenoma_",k,".ZicoSeq.RData" ))


## ---- Autism (together with calculateP_shuffle.sh, run_real_shuffle.sh)-----
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  dir = args[1]
  file = args[2]
  i = as.integer(args[3])
}

library(plyr)
library(dplyr)
library(tibble)
library(reshape2)
library(GUniFrac)

set.seed(i)
setwd(paste0('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/noshuffle/data/',dir))
norm ='GMPR'; method='wilocx'
cat(file,'\n')
load(file)
y <- sum.measurement
if(!is.null(norm)){
  if(norm =='RLE'){
    suppressMessages(dds <- DESeqDataSetFromMatrix(round(y), DataFrame(conditions), ~conditions))
    suppressMessages(dds <- estimateSizeFactors(dds))
    s <- sizeFactors(dds)
  }
  if(norm =='GMPR'){s <- GUniFrac::GMPR(y)}
  if(norm =='TMM'){s <- calcNormFactors(y, method = 'TMM')} # I think this function is not correctly written
  if(norm =='TSS'){s <- colSums(y)}
  count_norm <- y/s
}
if(method =='wilocx'){
  conditions <- sample(conditions)
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data <- cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    pvalues <- wilcox.test(gene~conditions, data)$p.value
    return(pvalues)
  })
  names(pvalues) <- rownames(count_norm)
}

if(method=='ZicoSeq'){
  meta.dat <- cbind.data.frame(conditions=sample(conditions))
  fit1 <- GUniFrac::ZicoSeq(meta.dat = meta.dat,feature.dat = as.matrix(sum.measurement), grp.name = 'conditions')
  pvalues <- fit1$p.raw
}
p <- pvalues
setwd('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/shuffle/P/')
if(!dir.exists(dir)){dir.create(dir)}
save(p, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/shuffle/P/',dir, '/P.',gsub('.Rdata','',file),'-iter',i,'.',method,'.RData'))



## change p to dat for later FDR methods fit
setwd('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/shuffle/P/Autism/')
files <- list.files()
idents <- unique(gsub('P.Autism_|-iter.*','',files))
for(iter in 1:100){
  cat(iter,'\n')
  P <- NULL
  load(paste0("P.Autism_",ident,"-iter",iter,".wilcox.RData"))
  P <- as.data.frame(p) %>% rownames_to_column('gene')
  colnames(P)[2] <- ident
  
  for(ident in idents[-1]){
    p <- NULL
    load(paste0("P.Autism_",ident,"-iter",iter,".wilcox.RData"))
    p <- as.data.frame(p) %>% rownames_to_column('gene')
    colnames(p)[2] <- ident
    P <- full_join(P, p)
  }
  pValues <- P
  suppressMessages(P.df <- reshape2::melt(pValues))
  pValues <- as.matrix(pValues %>% column_to_rownames('gene'))
  pvals <- P.df[,3]
  covariate <- P.df[,c(1,2),drop =F];colnames(covariate) <- c('gene','outcome')
  covariate$gene <- as.factor(covariate$gene)
  covariate$outcome <- as.factor(covariate$outcome)
  dat <- list(pValues = pValues, pvals = pvals, covariate =covariate)
  save(dat, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/shuffle/dat/P.Autism_',iter,'.wilcox.RData'))
  cat('\n')
}


# setwd('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/shuffle/dat/')
# files <- list.files(pattern = 'wilcox.RData')
# n <- c()
# for(file in files){
#   load(file)
#   n <- c(n, mean(dat$pValues<0.05, na.rm = T))
# }
# summary(n)



## calculate FDR (with run_real_shuffle.sh)
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  f = args[1]
  m = args[2]
}
source('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/code/func.R')
load(paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/shuffle/dat/',f))
methods_funs <- list('AdaptiveGBHo'='run_AdaptiveGBHo',
                     'AdaptiveGBHg'='run_AdaptiveGBHg',
                     'AdaptiveGBHo_storey'='run_AdaptiveGBHo_st',
                     'AdaptiveGBHg_storey'='run_AdaptiveGBHg_st',
                     'BH'='run_BH',
                     'ST'='run_ST',
                     'stratBHo'='run_stratBHo',
                     'stratSTo'='run_stratSTo',
                     'stratBHg'='run_stratBHg',
                     'stratSTg'='run_stratSTg',
                     'StoreyNew0'='Storey_New_0',
                     'StoreyNew05'='Storey_New_0.05',
                     'StoreyNew1'='Storey_New_0.1',
                     'StoreyNew2'='Storey_New_0.2',
                     'StoreyNew4'='Storey_New_0.4',
                     'StoreyGeo1'='Storey_Geo_0.1',
                     'StoreyAri1'='Storey_Ari_0.1',
                     'StoreyNew1N'='Storey_New_0.1_norm',
                     'StoreyNew1S'='Storey_New_0.1_S',
                     'tdgbh3F'='run_tdgbh3F',
                     'tdgbh3T'='run_tdgbh3T',
                     'swfdr'='run_swfdr',
                     'CAMT'='run_CAMT',
                     'FDRregT'='run_FDRregT',
                     'FDRregE'='run_FDRregE',
                     'adaptMT'='run_adaptMT')
tryCatch({
  res <- time <- NULL
  wrapper <- match.fun(methods_funs[[m]])
  out <- wrapper(dat)
  res <- sum(out$fdr<=0.05, na.rm = T)
  save(out, res, file = paste0('/home/mayo/m216453/lu/pairwiseFDR/data/realdata/shuffle/res/', gsub('.RData','',f), '_',m,'.Wilcox.Rdata'))
})


## ==== FDR plot code ======
setwd('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/shuffle/res/')
files <- list.files(pattern = 'Rdata$')
files1 <- files[grep('P\\.adenoma',files)];files1 <- files1[-grep('AdaptiveGBHg.ZicoSeq.Rdata|AdaptiveGBHo.ZicoSeq.Rdata|stratSTg|stratBHg',files1)]
files2 <- files[grep('P\\.combo',files)];files2 <- files2[-grep('AdaptiveGBHg.ZicoSeq.Rdata|AdaptiveGBHo.ZicoSeq.Rdata|stratSTg|stratBHg',files2)]
files3 <- files[grep('P\\.Autism',files)];files3 <- files3[-grep('AdaptiveGBHg.ZicoSeq.Rdata|AdaptiveGBHo.Wilcox.Rdata|stratSTg|stratBHg|AdaptiveGBHg.Wilcox.Rdata',files3)]
table(gsub('.*\\.wilcox\\_|.Wilcox.Rdata','',files3))
table(gsub('.*\\.wilcox\\_|.Wilcox.Rdata|.ZicoSeq.Rdata|.*\\.ZicoSeq\\_','',files2))
table(gsub('.*\\.wilcox\\_|.Wilcox.Rdata|.ZicoSeq.Rdata|.*\\.ZicoSeq\\_','',files1))

files <- c(files1, files2, files3)
ct0.01 <- ct0.05 <- ct0.1 <- ct0.15 <- ct0.2 <- c()
for(file in files){
  out<-NULL
  cat('.')
  load(file)
  ct0.01 <- c(ct0.01,sum(out$fdr<=0.01, na.rm = T))
  ct0.05 <- c(ct0.05,sum(out$fdr<=0.05, na.rm = T))
  # mn0.05 <- c(mn0.05,mean(out$fdr<=0.05, na.rm = T))
  ct0.1 <- c(ct0.1,sum(out$fdr<=0.1, na.rm = T))
  ct0.15 <- c(ct0.15,sum(out$fdr<=0.15, na.rm = T))
  ct0.2 <- c(ct0.2,sum(out$fdr<=0.2, na.rm = T))
}
names(ct0.01) <- names(ct0.05) <- names(ct0.1) <- names(ct0.15) <- names(ct0.2) <- files
df <- cbind(ct0.01,ct0.05,ct0.1,ct0.15,ct0.2)
colnames(df) <- c('FDR(0.01)','FDR(0.05)','FDR(0.1)','FDR(0.15)','FDR(0.2)')
df <- as.data.frame(df) %>% rownames_to_column('file')
df$dataset <- gsub('_.*|P.','',df$file)
df$method <- gsub('.Wilcox.Rdata|.ZicoSeq.Rdata|.*ZicoSeq_|.*Wilcox_|.*wilcox_','',df$file)
df$method[df$method=='StoreyNew1'] <- '2dGBH'
df$method[df$method=='stratBHo'] <- 'stratBH_o'
df$method[df$method=='stratBHg'] <- 'stratBH_g'
df$method[df$method=='stratSTo'] <- 'stratST_o'
df$method[df$method=='stratSTg'] <- 'stratST_g'
df$method[df$method=='AdaptiveGBHo_storey'] <- 'AdaptiveGBH_o'
df$method[df$method=='AdaptiveGBHg_storey'] <- 'AdaptiveGBH_g'
# save(df, file = '/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/adenomaZ_comboZ_autismZ(shuffle)_rev.Rdata')


# load('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/realdata/adenomaZ_comboZ_autismZ(shuffle)_rev.Rdata')
df <- melt(df)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
method.o <- c('2dGBH','BH','ST','AdaptiveGBH_g','AdaptiveGBH_o','stratST_o','stratBH_o')
df$dataset <- firstup(df$dataset)
df$iter <- ifelse(df$value==0,0,1)
fdr <- aggregate(iter ~ dataset + method + variable, df, function(x) mean(x))
mn <- aggregate(iter ~ method, fdr, function(x) mean(x))
mn <- mn[order(mn$iter),]
fdr <- within(fdr,  method<- factor(method, levels =mn$method))
fdr3 <- fdr[(fdr$method %in% method.o),]
fdr3$iter <- fdr3$iter *100
fdr3 <- within(fdr3,  dataset<- factor(dataset, levels =c("Combo","Adenoma","Autism")))

fdr3.sub <- fdr3[fdr3$dataset=='Adenoma' &fdr3$variable=='FDR(0.05)',]
ord <- levels(unique(fdr3.sub[order(fdr3.sub$iter),'method']))
p1 <- ggplot(fdr3[fdr3$variable=='FDR(0.05)' & fdr3$dataset %in% c("Autism","Adenoma","Combo"),], aes(x = reorder(method, iter), y = iter, fill = method)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values = col) + 
  facet_wrap(.~ dataset) + 
  geom_hline(yintercept = 5, color = brewer.pal(9,'BrBG')[1], linetype='dashed') +
  theme_bw() + theme(axis.text.x = element_text(color = 'black', angle = 90, hjust = 1, size = 14),
                     axis.text = element_text(color = 'black', size = 14),
                     axis.title = element_text(color = 'black', size = 14),
                     strip.text = element_text(color = 'black', size = 14),
                     legend.text = element_text(color = 'black', size = 14),
                     panel.grid = element_blank(),
                     legend.position = 'none') + 
  labs(x = '',y = '% of permuted dataset with positive findings', fill ='') 
p1
ggsave(file='/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/result/Figure4A_rev.pdf',width = 10, height =  6)


df3 <- df[df$variable %in% c('FDR(0.05)') & (df$method %in% method.o),]
df3 <- within(df3,  method<- factor(method, levels =ord))
df3 <- within(df3,  dataset<- factor(dataset, levels =c("Combo","Adenoma","Autism")))
combo.df <- df3[df3$dataset=="Combo",] %>% mutate(pct = value/79.18)
Adenoma.df <- df3[df3$dataset=="Adenoma",] %>% mutate(pct = value/70.84)
summary(Adenoma.df$pct[Adenoma.df$method =='AdaptiveGBH_o'])
summary(Adenoma.df$pct[Adenoma.df$method %in% c('stratBH_o','stratST_o')])
Autism.df <- df3[df3$dataset=="Autism",] %>% mutate(pct = value/1581.57)
summary(Autism.df$pct[Autism.df$method %in% c('AdaptiveGBH_g','stratST_o')])
summary(Autism.df$pct[Autism.df$method %in% c('BH','ST')])
summary(Autism.df$pct[Autism.df$method %in% c('2dGBH')])
summary(Autism.df$pct[Autism.df$method %in% c('AdaptiveGBH_o')])

p21 <- ggplot(df3[df3$dataset=="Combo",], aes(x = method, y = value, fill = method)) + 
  geom_boxplot(outlier.size = 0.4) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~./79.18 , name=""))+
  scale_fill_manual(values = col) + 
  facet_wrap(~ dataset, scales = 'free_y') + 
  theme_bw() + 
  theme(axis.text.x = element_text(color = 'black', angle = 90, hjust = 1, vjust = 0.25, size = 14), 
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none') + 
  labs(x = '',y = '# of positive findings in permuted datasets', fill ='') 
p22 <- ggplot(df3[df3$dataset=="Adenoma",], aes(x = method, y = value, fill = method)) + 
  geom_boxplot(outlier.size = 0.4) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~./70.84 , name=""))+
  scale_fill_manual(values = col) + 
  facet_wrap(~ dataset, scales = 'free_y') + 
  theme_bw() + 
  theme(axis.text.x = element_text(color = 'black', angle = 90, hjust = 1, vjust = 0.25, size = 14), 
        axis.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        panel.grid = element_blank(),
        legend.position = 'none') + 
  labs(x = '',y = '', fill ='') 
p23 <- ggplot(df3[df3$dataset=="Autism",], aes(x = method, y = value, fill = method)) + 
  geom_boxplot(outlier.size = 0.4) + 
  scale_y_continuous(sec.axis = sec_axis(trans=~./1581.57, name="% of positive findings from permuted data"))+#442.65
  scale_fill_manual(values = col) + 
  facet_wrap(~ dataset, scales = 'free_y') + 
  theme_bw() + 
  theme(axis.text.x = element_text(color = 'black', angle = 90, hjust = 1, vjust = 0.25, size = 14), 
        axis.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        strip.text = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        panel.grid = element_blank(),
        legend.position = 'none') + 
  labs(x = '',y = '', fill ='') 

library(ggpubr)
p2 <- ggpubr::ggarrange(p21 , p22 , p23, # remove axis labels from plots
                        labels = NULL, ncol = 3, nrow = 1,align = "hv")
p2
ggsave('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/result/Figure4B_rev.pdf',width = 10, height = 4)




