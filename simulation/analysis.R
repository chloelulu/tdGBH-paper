args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  f = args[1]
  m = args[2]
}


library(reshape2)
source('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/code/func.R')
dir <- gsub('.*\\/','',getwd())
output <- paste0("/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/result/",dir,'/')
load(f)
pValues <- sim$p.mat
truths <- sim$truth.mat
pvals <- melt(pValues)[,3]
covariate <- melt(pValues)[,c(1,2),drop =F];colnames(covariate) <- c('gene','outcome')
covariate$gene <- as.factor(covariate$gene)
covariate$outcome <- as.factor(covariate$outcome)
truth <- melt(truths)[,3]
dat <- list(truth = truth, pValues = pValues, pvals = pvals, covariate =covariate)
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
                     'FDRregE'='run_FDRregE')
tryCatch({
  wrapper <- match.fun(methods_funs[[m]])
  out <- wrapper(dat)
  res <- cal.fdr.tpr(fdr = out$fdr, truth = dat$truth, cutoffs=c(0.05))
  save(out, res, file = paste0(output, gsub('.RData','',f), '_',m,'.Rdata'))
})



