setwd('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/data/simulation/')
source('/research/bsi/projects/staff_analysis/m216453/pairwiseFDR/code/func.R')

### ==== NULL Setting ====
## independent
for(iter in 1:1000){
  for(affected.outcome in c(0.2)){ # c(0.1, 0.2, 0.4)
    for(affected.feature in c(0.2)){
      for(signal.density in c(0)){ #0.001, 0.01, 0.02
        for(outcome in c(20, 50, 200, 500)){#c(50, 200, 500)
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density
            space <- affected.feature * affected.outcome * outcome * feature
            cat(iter,':(',signal.density, '*',outcome,'*',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,'spaces)\n')
            affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature #*affected.feature
            outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
            sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='None',signal.strength = 2,
                                 feature.density=affected.feature, signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting0LMb_indep/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}

## AR(1)
for(iter in 1:1000){
  for(affected.outcome in c(0.2)){ # c(0.1, 0.2, 0.4)
    for(affected.feature in c(0.2)){
      for(signal.density in c(0)){ #0.001, 0.01, 0.02
        for(outcome in c(20, 50, 200, 500)){#c(50, 200, 500)
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density
            space <- affected.feature * affected.outcome * outcome * feature
            cat(iter,':(',signal.density, '*',outcome,'*',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,'spaces)\n')
            affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature #*affected.feature
            outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
            sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='ARI',signal.strength = 2,
                                 feature.density=affected.feature, signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting0LMb_ARI/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}

## Block
for(iter in 701:1000){
  gc()
  for(affected.outcome in c(0.2)){ # c(0.1, 0.2, 0.4)
    for(affected.feature in c(0.2)){
      for(signal.density in c(0)){ #0.001, 0.01, 0.02
        for(outcome in c(20, 50, 200, 500)){#c(50, 200, 500)
          for(feature in c(1000)){
            
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density
            space <- affected.feature * affected.outcome * outcome * feature
            cat(iter,':(',signal.density, '*',outcome,'*',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,'spaces)\n')
            affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature #*affected.feature
            outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
            sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='Block',signal.strength = 2,
                                 feature.density=affected.feature, signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting0LMb_BlockI/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}





### ---- independent-20 & 500 ----
## outcome
for(iter in 1:1000){
  for(affected.outcome in 0.2){
    for(affected.feature in 1){ 
      for(signal.density in c(0.01, 0.02, 0.05, 0.1)){
        for(outcome in c(20, 500)){
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density; space = feature.no * affected.feature * affected.outcome
            cat('signal density(',signal.density, ')* outcome(',outcome,')* feature(',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,')\n')
            affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature# *affected.feature
            cat('=> signal density for affected outcome is:',affected.outcome.signal.density,'\n')
            outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
            sim <- simulata.data(outcome =outcome, feature.no = feature.no, cor.struct ='None',signal.strength = 2,
                                 feature.density=affected.feature,signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting01LMo/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}
gc()

## gene
for(iter in 1:100){
  for(affected.outcome in c(1)){
    for(affected.feature in c(0.2)){
      for(signal.density in c(0.01, 0.02, 0.05, 0.1)){
        for(outcome in c(500)){#c(20, 500)
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density; space = feature.no * affected.feature * affected.outcome
            cat('signal density(',signal.density, ')* outcome(',outcome,')* feature(',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,')\n')
            affected.outcome.signal.density <- signal.density
            cat('=> signal density for affected outcome is:',affected.outcome.signal.density,'\n')
            outcome.signal.densities <- rep(signal.density, outcome * affected.outcome)
            sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='None',signal.strength = 2,
                                 feature.density=affected.feature, signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting01LMg/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}
gc()
## both-fix affected feature = 0.2, affected outcome = 0.2, check signal trend
for(iter in 1:100){
  for(affected.outcome in 0.2){
    for(affected.feature in 0.2){
      for(signal.density in c(0.005,0.01,0.02,0.04)){ #c(0.005,0.01,0.02,0.04)
        for(outcome in c(20)){#c(20, 500)
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density; space = feature.no * affected.feature * affected.outcome
            cat('signal density(',signal.density, ')* outcome(',outcome,')* feature(',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,')\n')
            affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature #*affected.feature
            cat('=> signal density for affected outcome is:',affected.outcome.signal.density,'\n')
            outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
            sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='None',signal.strength = 2,
                                 feature.density=affected.feature, signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting01LMb/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}
gc()



### ---- AR(1)-20 & 500 ----
## outcome
for(iter in 1:100){
  for(affected.outcome in 0.2){
    for(affected.feature in 1){ 
      for(signal.density in c(0.01, 0.02, 0.05, 0.1)){
        for(outcome in c(20, 500)){
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density; space = feature.no * affected.feature * affected.outcome
            cat('signal density(',signal.density, ')* outcome(',outcome,')* feature(',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,')\n')
            affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature# *affected.feature
            cat('=> signal density for affected outcome is:',affected.outcome.signal.density,'\n')
            outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
            sim <- simulata.data(outcome =outcome, feature.no = feature.no, cor.struct ='ARI',signal.strength = 2,
                                 feature.density=affected.feature,signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting03LMo/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}
gc()
## gene
for(iter in 1:100){
  for(affected.outcome in c(1)){
    for(affected.feature in c(0.2)){
      for(signal.density in c(0.01, 0.02, 0.05, 0.1)){
        for(outcome in c(20, 500)){
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density; space = feature.no * affected.feature * affected.outcome
            cat('signal density(',signal.density, ')* outcome(',outcome,')* feature(',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,')\n')
            affected.outcome.signal.density <- signal.density
            cat('=> signal density for affected outcome is:',affected.outcome.signal.density,'\n')
            outcome.signal.densities <- rep(signal.density, outcome * affected.outcome)
            sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='ARI',signal.strength = 2,
                                 feature.density=affected.feature, signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting03LMg/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}
gc()
## both-fix affected feature = 0.2, affected outcome = 0.2, check signal trend
for(iter in 1:100){
  for(affected.outcome in 0.2){
    for(affected.feature in 0.2){
      for(signal.density in c(0.005, 0.01,0.02,0.04)){ 
        for(outcome in c(20, 500)){
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density; space = feature.no * affected.feature * affected.outcome
            cat('signal density(',signal.density, ')* outcome(',outcome,')* feature(',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,')\n')
            affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature #*affected.feature
            cat('=> signal density for affected outcome is:',affected.outcome.signal.density,'\n')
            outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
            sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='ARI',signal.strength = 2,
                                 feature.density=affected.feature, signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting03LMb/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}
gc()

### ---- Block-20 & 500 ----
## outcome
for(iter in 1:100){
  for(affected.outcome in 0.2){
    for(affected.feature in 1){ 
      for(signal.density in c(0.01, 0.02, 0.05, 0.1)){
        for(outcome in c(20, 500)){
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density; space = feature.no * affected.feature * affected.outcome
            cat('signal density(',signal.density, ')* outcome(',outcome,')* feature(',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,')\n')
            affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature# *affected.feature
            cat('=> signal density for affected outcome is:',affected.outcome.signal.density,'\n')
            outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
            sim <- simulata.data(outcome =outcome, feature.no = feature.no, cor.struct ='BlockI',signal.strength = 2,
                                 feature.density=affected.feature,signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting02LMo/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}
gc()
## gene
for(iter in 1:100){
  for(affected.outcome in c(1)){
    for(affected.feature in c(0.2)){
      for(signal.density in c(0.01, 0.02, 0.05, 0.1)){
        for(outcome in c(20, 500)){
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density; space = feature.no * affected.feature * affected.outcome
            cat('signal density(',signal.density, ')* outcome(',outcome,')* feature(',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,')\n')
            affected.outcome.signal.density <- signal.density
            cat('=> signal density for affected outcome is:',affected.outcome.signal.density,'\n')
            outcome.signal.densities <- rep(signal.density, outcome * affected.outcome)
            sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='BlockI',signal.strength = 2,
                                 feature.density=affected.feature, signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting02LMg/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}
gc()
## both-fix affected feature = 0.2, affected outcome = 0.2, check signal trend
for(iter in 1:100){
  for(affected.outcome in 0.2){
    for(affected.feature in 0.2){
      for(signal.density in c(0.005, 0.01,0.02,0.04)){ 
        for(outcome in c(20, 500)){
          for(feature in c(1000)){
            cat('---------------------\n')
            feature.no = feature * outcome; signals <- feature.no * signal.density; space = feature.no * affected.feature * affected.outcome
            cat('signal density(',signal.density, ')* outcome(',outcome,')* feature(',feature,')=',signals, 'signals distributed in ',affected.outcome,'outcomes and',affected.feature,'features(',space,')\n')
            affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature #*affected.feature
            cat('=> signal density for affected outcome is:',affected.outcome.signal.density,'\n')
            outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
            sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='BlockI',signal.strength = 2,
                                 feature.density=affected.feature, signal.densities = outcome.signal.densities)
            cat('[',mean(colSums(sim$truth.mat)>0),';');cat(mean(rowSums(sim$truth.mat)>0),';');cat(sum(sim$truth.mat),']\n')
            save(sim, file = paste0('setting02LMb/iter',iter,'_AO',affected.outcome,'_AF',affected.feature,'_SD',signal.density,'_o',outcome,'_f',feature,'.RData'))
          }
        }
      }
    }
  }
}
gc()
