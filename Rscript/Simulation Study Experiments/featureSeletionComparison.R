rm(list=ls())
gc()

### Load/Install Required Packages  ####
library(compositions)
library(data.table)
library(reshape2)
library(doParallel)
library(igraph)
library(caret)
library(tidyverse)
library(glmnet)
library(vegan)
library(keras)
library(PRROC)
library(energy)
library(diffCompVarRcpp)
library(selEnergyPermR)


## Start Cluster
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)



## Load Helper Functions
source(file = "Helper Functions/functions.R")


######===========================================================-
### Read External Arguments
######===========================================================-
args = c(10,1,100)
args = commandArgs(trailingOnly = TRUE)
seed_ = as.numeric(args[1]) # random seed selection
sampleDesign = as.numeric(args[2]); design = if_else(sampleDesign==1,"Balanced","Unbalanced")
dims = as.numeric(args[3])
######===========================================================-

system.time({
## Read Data
if(sampleDesign==1){
  g1 = 20
  g2 = 60
}else{
  g1 = 40
  g2 = 40
}


## Define Feature Data Frame
performance.df = data.frame()
trainModels  = T


for(s in 1:5){

  switch(s,
         {
           ## Scenario-1 = xxx
           dat = selEnergyPermR::scenario1(seed = seed_,n1 = g1,n2 = g2,dms_ = dims);
           fname = "scenario1"
         },

         {
           ## Scenario-2 = Sim from WGS Data
           dat = selEnergyPermR::scenario2(seed = seed_,n1 = g1,n2 = g2,dms_ = dims);
           fname = "scenario2"
         },

         {
           ## Scenario-3 = Sim from additive log normal
           dat = selEnergyPermR::scenario3(seed = seed_,n1 = g1,n2 = g2,dms_ = dims);
           fname = "scenario3"
         },

         {
           ## Scenario-4 = Sim from sparse dirichilet
           dat = selEnergyPermR::scenario4(seed = seed_,n1 = g1,n2 = g2,dms_ = dims);
           fname = "scenario4"
         },

         {
           ## Scenario-4 = Sim from sparse dirichilet
           dat = selEnergyPermR::scenario5(seed = seed_,n1 = g1,n2 = g2,dms_ = dims);
           fname = "scenario5"
         }
  )


  ## Sparisty Pre-Process
  procData = selEnergyPermR::processCompData(dat,minPrevalence = .9)
  dat = procData$processedData
  impFact = procData$impFactor
  minorityClass = procData$minClss
  majorityClass = procData$majClass
  impClass = minorityClass
  Scenario = fname

  ## Pre Process -  Removbe col with all 0's
  y = dat[,-1]
  bool = colSums(y)==0
  y = y[,!bool]
  dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact))

  ## Get logratios
  time_lr = system.time({
    dat.f = data.frame(Status = factor(dat[,1]), fastImputeZeroes( compositions::clo(dat[,-1]) , impFactor  = impFact ) )
    lrs = calcLogRatio(dat.f)
    class = dat[,1]
  })





  models = c("pls")
  total_Ratios = ncol(lrs[,-1])

  ## all features
  Perf = data.frame(method = "All",Dims = dims,totalRatios = total_Ratios,
                    Scenario = s,seed = seed_,numRatios = total_Ratios,
                    compTime_sec = time_lr[3])
  bm = benchMarkPerformance_features(tbl = lrs,ensemble = models)
  ph = bm$performance
  Perf.all = cbind(Perf,ph)




  ##-------------------------*
  ## Boruta
  ##--------------------------*
  compTime = system.time({
    b = Boruta::Boruta(x = lrs[,-1],y = class,doTrace = 2,maxRun = 100,getImp = Boruta::getImpExtraGini)
    dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
    keep = dec %>%
      filter(Decision!="Rejected")
    kr =as.character(keep$Ratio)
  })

  feats.df = subset(lrs,select = c("Status",kr))

  if(ncol(feats.df)>2){
    ## Performance
    Perf = data.frame(method = "Boruta",Dims = dims,totalRatios = total_Ratios,
                      Scenario = s,seed = seed_,numRatios = ncol(feats.df[,-1]),compTime_sec = time_lr[3]+compTime[3])
    bm = benchMarkPerformance_features(tbl = feats.df,ensemble = models)
    ph = bm$performance
    Perf = cbind(Perf,ph)
  }else{
    Perf = Perf.all
    Perf$method = "Boruta"
    Perf$numRatios = ncol(lrs[,-1])
  }
  performance.df = rbind(performance.df,Perf)


  ##-------------------------*
  ## LASSO
  ##-------------------------*
  compTime = system.time({
    cv.lasso <- cv.glmnet(as.matrix(lrs[,-1]), lrs[,1], family='binomial', alpha=1,
                          #nfolds = numFold_glm,nlambda = num_lambda,
                          parallel=TRUE, type.measure='auc')
    tt = glmnet::glmnet(as.matrix(lrs[,-1]), lrs[,1], family='binomial', alpha=1,lambda = cv.lasso$lambda.min)
    df_coef = as.matrix(tt$beta)
    impVar = df_coef[df_coef[, 1] != 0, ]
    impVar_names = names(impVar)
  })
  feats.df = subset(lrs,select = c("Status",impVar_names))

  if(ncol(feats.df)>2){
    ## Performance
    Perf = data.frame(method = "LASSO",Dims = dims,totalRatios = total_Ratios,
                      Scenario = s,seed = seed_,numRatios = ncol(feats.df[,-1]),compTime_sec = time_lr[3]+compTime[3])
    bm = benchMarkPerformance_features(tbl = feats.df,ensemble = models)
    ph = bm$performance
    Perf = cbind(Perf,ph)
  }else{
    Perf = Perf.all
    Perf$method = "LASSO"
    Perf$numRatios = ncol(lrs[,-1])
  }
  performance.df = rbind(performance.df,Perf)





  ##-------------------------*
  ## RFE
  ##-------------------------*
  # define the control using a random forest selection function
  compTime = system.time({
    control <- rfeControl(functions=rfFuncs, method="cv", number=5)
    mx = floor(log(ncol(lrs[,-1]),2))
    sets = 2^(1:mx)
    # run the RFE algorithm
    results <- rfe(lrs[,-1], lrs[,1], sizes=sets,rfeControl=control)
    #print(results)
    kr = predictors(results)
  })

  ## Select Features
  feats.df = subset(lrs,select = c("Status",kr))

  if(ncol(feats.df)>2){
    ## Performance
    Perf = data.frame(method = "RFE",Dims = dims,totalRatios = total_Ratios,
                      Scenario = s,seed = seed_,numRatios = ncol(feats.df[,-1]),compTime_sec = time_lr[3]+compTime[3])
    bm = benchMarkPerformance_features(tbl = feats.df,ensemble = models)
    ph = bm$performance
    Perf = cbind(Perf,ph)
  }else{
    Perf = Perf.all
    Perf$method = "RFE"
    Perf$numRatios = ncol(lrs[,-1])
  }

  performance.df = rbind(performance.df,Perf)



  ##-------------------------*
  ## Filter Method
  ##-------------------------*
  compTime = system.time({
    ig = FSelectorRcpp::information_gain(x = lrs[,-1],y = as.factor(lrs[,1]),type = "infogain")$importance
    weights = data.frame(Ratio = colnames(lrs)[-1],importance = ig)
    subset <- weights[weights$importance>0,1]
  })

  ## Select Features
  feats.df = subset(lrs,select = c("Status",subset))

  if(ncol(feats.df)>2){
    ## Performance
    Perf = data.frame(method = "InformationGain",Dims = dims,totalRatios = total_Ratios,
                      Scenario = s,seed = seed_,numRatios = ncol(feats.df[,-1]),compTime_sec = time_lr[3]+compTime[3])
    bm = benchMarkPerformance_features(tbl = feats.df,ensemble = models)
    ph = bm$performance
    Perf = cbind(Perf,ph)
  }else{
    Perf = Perf.all
    Perf$method = "InformationGain"
    Perf$numRatios = ncol(lrs[,-1])
  }

  performance.df = rbind(performance.df,Perf)

  ##-------------------------*
  ### SelPermEnergy
  ##---------------------------*
  compTime = system.time({
    sep.true = selectionEnergy.scaled(inputData = dat,nreps_energy = 5e3)
  })
  feats.df = sep.true$finalSubset

  ## Select Features
  if(ncol(feats.df)>2){
    ## Performance
    Perf = data.frame(method = "selEnergyPerm",Dims = dims,totalRatios = total_Ratios,
                      Scenario = s,seed = seed_,numRatios = ncol(feats.df[,-1]),compTime_sec = compTime[3])
    bm = benchMarkPerformance_features(tbl = feats.df,ensemble = models)
    ph = bm$performance
    Perf = cbind(Perf,ph)
  }else{
    Perf = Perf.all
    Perf$method = "selEnergyPerm"
    Perf$numRatios = ncol(lrs[,-1])
  }
  performance.df = rbind(performance.df,Perf)


}


})


## Write Results
setwd("Results/featureSelection/")
csv_name = paste0("seed_",seed_,"-",design,"-dims_",dims,"-.csv")
write_csv(x = performance.df,path = csv_name)
