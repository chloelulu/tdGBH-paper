## === 2dGBH ====
require(structSSI)
weight_BH2 <- function(p.mat, pi0.method = 'storey', #c('lsl','tst', 'storey')
                       global.pi0.method = 'storey', shrink = 0.1){#c('lsl','tst', 'storey')
  
  pi0.method <- match.arg(pi0.method)
  global.pi0.method  <- match.arg( global.pi0.method )
  
  pi0 <- structSSI::estimate.pi0(as.vector(p.mat), method =  global.pi0.method)
  
  pi0.o <- apply(p.mat, 2, function(x) structSSI::estimate.pi0(x, method = pi0.method))
  pi0.g <- apply(p.mat, 1, function(x) structSSI::estimate.pi0(x, method = pi0.method))
  
  pi0.o <- (1 - shrink) * pi0.o + shrink * pi0
  pi0.g <- (1 - shrink)  * pi0.g + shrink * pi0
  
  pi0.o.mat <- t(matrix(pi0.o, nrow = length(pi0.o), ncol = length(pi0.g)))
  pi0.g.mat <- matrix(pi0.g, nrow = length(pi0.g), ncol = length(pi0.o))
  
  sd.r <- (sd(pi0.o) / sqrt(length(pi0.o)))  / (sd(pi0.o) / sqrt(length(pi0.o)) + sd(pi0.g) / sqrt(length(pi0.g)))
  pi0.og.mat <- sqrt((pi0.o.mat^(2 * sd.r)) * (pi0.g.mat^(2 * (1 - sd.r))))

  pi0 <- mean(pi0.og.mat)
  ws.og.mat <- (1 - pi0.og.mat) / pi0.og.mat
  p.ws.mat <- p.mat / ws.og.mat * (1 - pi0)
  p.adj <- matrix(p.adjust(as.vector(p.ws.mat), 'BH'), length(pi0.g), length(pi0.o))
  
  return(list(p.adj = p.adj, pi0 = pi0, pi0.o = pi0.o.mat, pi0.g = pi0.g.mat, pi0.og = pi0.og.mat))
}


## ===== simulation func =====
sim.NULL <- function(genes, signal.strength, feature.density, signal.density){
  
  H0 <- 1-pnorm(rnorm(round(genes * (1-signal.density)))) # one side pvalues
  H1 <- 1-pnorm(rnorm(round(genes * signal.density), signal.strength, sd = 1))
  H <- c(H0, H1)
  
  # select TRUE and FALSE, together consumes feature_density portion
  if (signal.density > feature.density){
    cat('Please provide a large feature density!')
  }else{
    # i.e.: feature.density = 0.5, signal.density = 0.1, genes = 100. 10 H1 genes randomly exists in 50 genes. Other 50 genes are 
    idx <- sample(length(H0))
    idx1 <- idx[1:((feature.density - signal.density) * genes)] # select Ps in H0 by feature density(i.e., feature.density=0.7 means signals exist in specific 70% features)
    idx2 <- idx[!(idx %in% idx1)] #edit 11/30/2022
    sub.P <- c(H1, H0[idx1]) # combine signals in H1 and selected H0[exists in specific genes].
    sub.truth <- c(rep(TRUE,length(H1)), rep(FALSE,length(idx1)))
    sub.idx <- sample(1:length(sub.P)) # shuffle the enriched genes. Let %signals[signal density] randomly exists in % features[feature density]
    P <- c(sub.P[sub.idx], H0[idx2])
    truth <- c(sub.truth[sub.idx], rep(FALSE,length(idx2)))
    return(list(P=P, truth = truth))
  }
}

## Within block correlation is 0.7
sim.AR2 <- function(genes, blocks, cor.rho, signal.strength, feature.density, signal.density){
  # simulate one outcome with #gene into #blocks
  sub.genes <- genes/blocks
  
  # prepare means for signals and non-signals
  mu.H1 <- rep(signal.strength, round(blocks*signal.density))
  mu.H0.1 <- rep(0, round(blocks*(feature.density-signal.density)))
  mu.H0.2 <- rep(0, round(blocks*(1-feature.density)))
  mu <- c(sample(c(mu.H1, mu.H0.1), length(c(mu.H1, mu.H0.1))), mu.H0.2)
  
  # prepare truth
  truth <- rep(TRUE, length(mu))
  truth[mu==0] <- FALSE
  
  # simulate data
  sim <- numeric()
  Sigma <- cor.rho ^ abs(outer(1:sub.genes, 1:sub.genes, "-"))
  sim <- c(mapply(function(mu_i) {
    MASS::mvrnorm(n=1, mu=rep(mu_i, sub.genes), Sigma=Sigma)
  }, mu_i = mu))
  
  
  # truth repeated to match genes
  truth <- rep(truth, each=sub.genes)
  P <- 1-pnorm(sim)
  return(list(P=P, truth = truth))
}

sim.ARI2 <- function(genes, blocks, cor.rho, signal.strength, feature.density, signal.density){
  # simulate one outcome with #gene into #blocks
  sub.genes <- genes/blocks
  
  # simulate data
  Sigma <- cor.rho ^ abs(outer(1:sub.genes, 1:sub.genes, "-"))
  
  # Generate multiple multivariate normal vectors at once
  sim <- c(t(replicate(blocks, MASS::mvrnorm(n=1, mu=rep(0, sub.genes), Sigma=Sigma))))
  z1 <- c(t(replicate(blocks, MASS::mvrnorm(n=1, mu=rep(2, sub.genes), Sigma=Sigma))))
  
  # prepare truth
  truth <- rep(FALSE, genes)
  idx.H1 <- sample(round(feature.density * genes), round(signal.density * genes)) # sample signals from %features(control by feature.density)
  sim[idx.H1] <- z1[idx.H1]
  truth[idx.H1] <- TRUE
  
  P <- 1-pnorm(sim)
  
  return(list(P=P, truth = truth))
}

sim.Block2 <- function(genes, blocks, cor.rho, signal.strength, feature.density, signal.density){
  # simulate one outcome with #gene into #blocks
  sub.genes <- genes/blocks
  
  # prepare means for signals and non-signals
  mu.H1 <- rep(signal.strength, round(blocks*signal.density))
  mu.H0.1 <- rep(0, round(blocks*(feature.density-signal.density)))
  mu.H0.2 <- rep(0, round(blocks*(1-feature.density)))
  mu <- c(sample(c(mu.H1, mu.H0.1), length(c(mu.H1, mu.H0.1))), mu.H0.2)
  
  # prepare truth
  truth <- rep(TRUE, length(mu))
  truth[mu==0] <- FALSE
  
  # simulate data
  sim <- numeric()
  Sigma <- diag(sub.genes)
  Sigma[, ] <- cor.rho
  diag(Sigma) <- 1
  sim <- c(mapply(function(mu_i) {
    MASS::mvrnorm(n=1, mu=rep(mu_i, sub.genes), Sigma=Sigma)
  }, mu_i = mu))
  
  # truth repeated to match genes
  truth <- rep(truth, each=sub.genes)
  P <- 1-pnorm(sim)
  return(list(P=P, truth = truth))
}

sim.BlockI2 <- function(genes, blocks, cor.rho, signal.strength, feature.density, signal.density){
  # simulate one outcome with #gene into #blocks
  sub.genes <- genes/blocks
  
  # simulate data
  Sigma <- diag(sub.genes)
  Sigma[, ] <- cor.rho
  diag(Sigma) <- 1
  
  # Generate multiple multivariate normal vectors at once
  sim <- c(t(replicate(blocks, MASS::mvrnorm(n=1, mu=rep(0, sub.genes), Sigma=Sigma))))
  z1 <- c(t(replicate(blocks, MASS::mvrnorm(n=1, mu=rep(2, sub.genes), Sigma=Sigma))))
  
  # prepare truth
  truth <- rep(FALSE, genes)
  idx.H1 <- sample(round(feature.density * genes), round(signal.density * genes)) # sample signals from %features(control by feature.density)
  sim[idx.H1] <- z1[idx.H1]
  truth[idx.H1] <- TRUE
  
  P <- 1-pnorm(sim)
  
  return(list(P=P, truth = truth))
}

# signal density:pct of genes rejects null hypotheses in each outcome.
# Thus if I fix signal.density=1, all blocks have signals in this outcome; if signal.density = 0.2, 20% blocks(genes) have signals.
# If I do not sample the index of this outcome, 20% signals will distributed in the first 20% genes.
# If I sample the index of this outcome, 20% signals will be randomly distributed in this outcome.
# If each column has the same signal, then signal.density = 0.2 means 20% genes distributed in 20% genes
simulata.data <- function(outcome,          # an integer, the number of outcomes to be simulated
                          feature.no,       # an integer, total hypothesis to be simulated
                          blocks = 100,
                          cor.struct = c('None','AR','Block','BlockI','ARI'),
                          cor.rho = 0.7,
                          signal.strength,  # a number, signify the strength of signal, the larger number the more small pvalues
                          signal.densities, # a vector with values between 0-1, represents the proportion of features are significant. 
                          # Each value in the vector is the signal density of each outcome. 
                          # For example, 0.3 in the vector means 30% features are significant(H01) and 70% are not significant(H1).
                          feature.density=1 # a number between 0-1, signals resides in a proportion of features. Should be larger than signal.density.
                          # For example, 0.7 means signal exists in a specific 70% of the features
){
  
  genes <- feature.no/outcome
  
  res <- lapply(signal.densities, function(signal.density) {
    if(cor.struct =='None'){
      return(sim.NULL(genes, signal.strength, feature.density, signal.density))
    }
    
    if(cor.struct == 'AR'){
      return(sim.AR2(genes, blocks, cor.rho, signal.strength, feature.density, signal.density))
    }
    
    if(cor.struct == 'ARI'){
      return(sim.ARI2(genes, blocks, cor.rho, signal.strength, feature.density, signal.density))
    }
    
    if(cor.struct == 'Block'){
      return(sim.Block2(genes, blocks, cor.rho, signal.strength, feature.density, signal.density))
    }
    
    if(cor.struct == 'BlockI'){
      return(sim.BlockI2(genes, blocks, cor.rho, signal.strength, feature.density, signal.density))
    }
    
  })
  
  pValues <- do.call(cbind, lapply(res, function(x) x$P))
  truths <- do.call(cbind, lapply(res, function(x) x$truth))
  
  return(list(p.mat = pValues, truth.mat = truths))
}





## ===== competing methods func =====
require(qvalue)
require(structSSI)
run_BH <-  function(dat){
  pvals <- dat$pvals
  t1 <- Sys.time()
  fdr <- p.adjust(pvals, method = 'BH')
  time <- difftime(Sys.time(),t1)
  return(list(fdr=fdr, time = time))
}

run_stratBHo <-  function(dat){
  pValues <- dat$pValues
  t1 <- Sys.time()
  fdr <- apply(pValues, 2, function(x) p.adjust(x, 'BH'))
  time <- difftime(Sys.time(),t1)
  return(list(fdr=fdr, time = time))
}

run_stratBHg <-  function(dat){
  pValues <- dat$pValues
  t1 <- Sys.time()
  fdr <- apply(pValues, 1, function(x) p.adjust(x, 'BH'))
  time <- difftime(Sys.time(),t1)
  return(list(fdr=fdr, time = time))
}

run_ST <- function(dat){
  pvals <- dat$pvals
  require(qvalue)
  t1 <- Sys.time()
  fdr <- qvalue::qvalue(pvals)$qvalues
  time <- difftime(Sys.time(),t1)
  return(list(fdr=fdr, time = time))
}

run_stratSTo <-  function(dat){
  pValues <- dat$pValues
  require(qvalue)
  t1 <- Sys.time()
  #  https://github.com/StoreyLab/qvalue/issues/20
  fdr <- apply(pValues, 2, function(x){
    tryCatch({qvalue(x)$qvalue}, error = function(err){return(qvalue(x, pi0=1)$qvalue)})
  })  
  time <- difftime(Sys.time(),t1)
  return(list(fdr=fdr, time = time))
}

run_stratSTg <-  function(dat){
  pValues <- dat$pValues
  # idx <- (apply(pValues, 1, function(x) length(x[!is.na(x)])) >1)
  # pValues <- pValues[idx,]
  require(qvalue)
  t1 <- Sys.time()
  #  https://github.com/StoreyLab/qvalue/issues/20
  fdr <- apply(pValues, 1, function(x){
    tryCatch({qvalue(x)$qvalue}, error = function(err){return(qvalue(x, pi0=1)$qvalue)})
  })  
  time <- difftime(Sys.time(),t1)
  return(list(fdr=fdr, time = time))
}

run_AdaptiveGBHo <-  function(dat){
  pvals <- dat$pvals
  covariate <- dat$covariate$outcome
  idx <- which(!is.na(pvals))
  pvals <- pvals[idx]
  covariate <- droplevels(covariate[idx])
  require(structSSI)
  t1 <- Sys.time()
  res <- structSSI::Adaptive.GBH(unadj.p=pvals, group.index=covariate)
  idx <- order(res@p.vals$hypothesisIndex)
  res <- res@p.vals[idx,]
  
  fdr <- dat$pvals
  fdr[which(!is.na(pvals))] <- res$adjp
  
  time <- difftime(Sys.time(),t1)
  return(list(fdr=fdr, time = time))
}

run_AdaptiveGBHg <-  function(dat){
  pvals <- dat$pvals
  covariate <- dat$covariate$gene
  idx <- which(!is.na(pvals))
  pvals <- pvals[idx]
  covariate <- droplevels(covariate[idx])
  require(structSSI)
  t1 <- Sys.time()
  res <- structSSI::Adaptive.GBH(unadj.p=pvals, group.index=covariate)
  idx <- order(res@p.vals$hypothesisIndex)
  res <- res@p.vals[idx,]
  
  fdr <- dat$pvals
  fdr[which(!is.na(pvals))] <- res$adjp
  
  time <- difftime(Sys.time(),t1)
  return(list(fdr=fdr, time = time))
}

Storey_New_0.1 <- function(dat, shrink = 0.1){
  t1 <- Sys.time()
  fdr <- weight_BH2(p.mat = dat$pValues, pi0.method = 'storey', 
                    global.pi0.method = 'storey', shrink = shrink)
  time <- difftime(Sys.time(),t1)
  return(list(fdr=as.vector(fdr$p.adj), time = time))
}




## ===== Other func =====
cal.fdr.tpr <- function(fdr, truth, cutoffs){
  FDR <- TPR <- c()
  for(cutoff in cutoffs){
    tp <- sum((fdr <= cutoff) & truth, na.rm = T)
    fn <- sum((fdr > cutoff)  & truth, na.rm = T)
    fp <- sum((fdr <= cutoff) & !truth, na.rm = T)
    TPR <- c(TPR, ifelse(!tp, 0, tp/(tp+fn)))
    FDR <- c(FDR, ifelse(!fp, 0, fp/(tp+fp)))
  }
  names(FDR) <- names(TPR) <- cutoffs
  return(list(fdr = FDR, tpr = TPR))
}


