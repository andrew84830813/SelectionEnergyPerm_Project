library(doParallel)
library(igraph)
library(caret)
library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(compositions)
library(diffCompVarRcpp)
library(selEnergyPermR)



## Start Local Cluster
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)


######===========================================================-
### Read External Arguments
######===========================================================-
args = c(1,8)
args = commandArgs(trailingOnly = TRUE)
seed_ = as.numeric(args[1]) # random seed selection
ds = as.numeric(args[2])
######===========================================================-




blocks = NULL
om = NULL
nfolds = 1
switch(ds,

       {
         dat = read_csv(file = "Data/deliveryMode_month1.csv");
         fname = "birthDelivery_month1";
         ## get sample names
         snames = dat$sample_name;
         dat = data.frame(Status = dat$delivery,dat[,5:ncol(dat)]);
         dat = data.frame(dat);
         rownames(dat) = snames;
         ## fix column names
         colnames(dat) = str_replace(colnames(dat),pattern  = "__\\.",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = "#");
         colnames(dat) = str_replace(colnames(dat),pattern  = "f__",replacement = "");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "#g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "g__",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = ".")

       }, #1 - birth study

       {
         dat = read_csv(file = "Data/deliveryMode_month2.csv");
         fname = "birthDelivery_month2";
         ## get sample names
         snames = dat$sample_name;
         dat = data.frame(Status = dat$delivery,dat[,5:ncol(dat)]);
         dat = data.frame(dat);
         rownames(dat) = snames;
         ## fix column names
         colnames(dat) = str_replace(colnames(dat),pattern  = "__\\.",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = "#");
         colnames(dat) = str_replace(colnames(dat),pattern  = "f__",replacement = "");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "#g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "g__",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = ".")
       }, #2

       {
         dat = read_csv(file = "Data/deliveryMode_month3.csv");
         fname = "birthDelivery_month3";
         ## get sample names
         snames = dat$sample_name;
         dat = data.frame(Status = dat$delivery,dat[,5:ncol(dat)]);
         dat = data.frame(dat);
         rownames(dat) = snames;
         ## fix column names
         colnames(dat) = str_replace(colnames(dat),pattern  = "__\\.",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = "#");
         colnames(dat) = str_replace(colnames(dat),pattern  = "f__",replacement = "");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "#g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "g__",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = ".")
       }, #3

       {
         dat = read_csv(file = "Data/deliveryMode_month4.csv");
         fname = "birthDelivery_month4";
         ## get sample names
         snames = dat$sample_name;
         dat = data.frame(Status = dat$delivery,dat[,5:ncol(dat)]);
         dat = data.frame(dat);
         rownames(dat) = snames;
         ## fix column names
         colnames(dat) = str_replace(colnames(dat),pattern  = "__\\.",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = "#");
         colnames(dat) = str_replace(colnames(dat),pattern  = "f__",replacement = "");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "#g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "g__",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = ".")
       }, #4

       {
         dat = read_csv(file = "Data/deliveryMode_month5.csv");
         fname = "birthDelivery_month5";
         ## get sample names
         snames = dat$sample_name;
         dat = data.frame(Status = dat$delivery,dat[,5:ncol(dat)]);
         dat = data.frame(dat);
         rownames(dat) = snames;
         ## fix column names
         colnames(dat) = str_replace(colnames(dat),pattern  = "__\\.",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = "#");
         colnames(dat) = str_replace(colnames(dat),pattern  = "f__",replacement = "");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "#g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "g__",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = ".")
       }, #5

       {
         dat = read_csv(file = "Data/deliveryMode_month6.csv");
         fname = "birthDelivery_month6";
         ## get sample names
         snames = dat$sample_name;
         dat = data.frame(Status = dat$delivery,dat[,5:ncol(dat)]);
         dat = data.frame(dat);
         rownames(dat) = snames;
         ## fix column names
         colnames(dat) = str_replace(colnames(dat),pattern  = "__\\.",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = "#");
         colnames(dat) = str_replace(colnames(dat),pattern  = "f__",replacement = "");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "#g__",replacement = ".");
         colnames(dat) = str_replace(colnames(dat),pattern  = "g__",replacement = "__");
         colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = ".")
       }, #6

       {
         ## ds  = 7
         dat = read_csv(file = "Data/Uguanda_PIH.csv");
         fname = "pih";
       }, #7

       {
         dat = read_csv(file = "Data/diabimmune_f6.csv");
         fname = "diabimmune_f6";
         md = read_csv(file = "Data/diabimmune_f6_metadata.csv");
         blocks = md$permutation_blocks
         nfolds =1
        }, #8 - diabimmune f6

       {
         dat = read_csv(file = "Data/diabimmune_6-12.csv");
         fname = "diabimmune_6-12";
         md = read_csv(file = "Data/diabimmune_6-12_metadata.csv");
         blocks = md$permutation_blocks
         nfolds = 1
       }, #9 - diabimmune 6-12

       {
         dat = read_csv(file = "Data/diabimmune_12-18.csv");
         fname = "diabimmune_12-18";
         md = read_csv(file = "Data/diabimmune_12-18_metadata.csv");
         blocks = md$permutation_blocks
         nfolds = 1
       }, #10 - diabimmune 12-18

       {
         md = read_csv(file = "Data/diabimmune_18-24_metadata.csv");
         fname = "diabimmune_18-24";
         dat = read_csv(file = "Data/diabimmune_18-24.csv");
         blocks = md$permutation_blocks
         nfolds = 1
       }, ##11 - diabimmune 18-24

       {
         dat = read_csv(file = "Data/fecalcalprotecin_UC.csv");
         fname = "fecalCalProtectin_UC";
         md = read_csv(file = "Data/fecalcalprotecin_UC_metadata.csv");
         blocks = md$visit_num ;
         om = "combinedF"
       }, #12 - fecal calprotecin UC

       {
         dat = read_csv(file = "Data/fecalcalprotecin_CD.csv");
         fname = "fecalCalProtectin_CD";
         md = read_csv(file = "Data/fecalcalprotecin_CD_metadata.csv");
         blocks = md$visit_num ;
         om = "combinedF"
       }, #13 - fecal calprotecin CD

       {
         dat = read_csv(file = "Data/fecalcalprotecin_ALL.csv");
         fname = "fecalCalProtectin_ALL";
         md = read_csv(file = "Data/fecalcalprotecin_ALL_metadata.csv");
         blocks = md$perm_blocks;
         om = "combinedF"
       } #14 - fecal calprotecin All

)



## PRe Process Data
dat = data.frame(dat)
table(dat[,1])
maxSparisty = .9
procData = processCompData(dat,minPrevalence = maxSparisty)
dat = procData$processedData
impFact = procData$impFactor
minorityClass = procData$minClss
majorityClass = procData$majClass
impClass = "pos"


## Pre Process -  Removbe col with all 0's
y = dat[,-1]
bool = colSums(y)==0
y = y[,!bool]
dat = data.frame(Status = dat[,1],
                 fastImputeZeroes( compositions::clo(y),impFactor = impFact))
message("preProcess Complete")




# SelEnergyPerm -----------------------------------------------------------

## Set Common Permutations
permRep = 1000
library(permute)
if(is.null(blocks)){
  set.seed(08272008)
  seeds = shuffleSet(nrow(dat),nset = permRep)
  seeds = data.frame(t(seeds))
}else{
  set.seed(08272008)
  h1 <- how(blocks = blocks)
  ## Set Common Permutations
  seeds = shuffleSet(nrow(dat),nset = permRep,control = h1)
  seeds = data.frame(t(seeds))
}
b = 100*seed_-99
ee = 100*seed_
seeds = seeds[,b:ee]
permRep = ncol(seeds)

## run SelEnergyPerm
px = 500
sep.true = selEnergyPermR::selectionEnergy.scaled(inputData = dat,
                                  nreps_energy = 1000,
                                  dcv_nfold = nfolds,
                                  optimizationMetric = om,
                                  patience = px)
x = sep.true$retainedDCV
## Get final features ####
feature.df = sep.true$finalSubset
## Retiecve DCV ####
lrs_dcv = sep.true$retainedDCV
rxDCV = sep.true$rawDCV[,2:6]
## compute Metric on Final Subet ####
set.seed(08272008)
spEPerf = selEnergyPermR::featureSlectionPerformance(tbl = feature.df,Method_Name = "spE",
                                                     plot_ = F,nreps_energy = 1e6)
spEPerf.df = spEPerf$performance
targetRatios = spEPerf.df$NumRatios
optMetric = sep.true$optimization_Metric
if(optMetric=="combinedF"){
  testStat = spEPerf.df$combinedF
}else{
  testStat = spEPerf.df$EnergyF
}
testStat.df =  data.frame(Type = "empirical",
                          Seed = seed_,
                          fname,rep = 1,
                          num_ratios = ncol(feature.df),
                          tstat = testStat,
                          optMetric)


## Null Results ####
nullPerf = data.frame()
nreps = 10
nreps = permRep
null.testStat = data.frame()
system.time({
  for(r in 1:nreps){

    suppressMessages(suppressWarnings({

      ## Run selection Energy Routine on permuted data
      null_data = dat
      null_data[,1] = null_data[seeds[,r],1]
      sep.null = selEnergyPermR::selectionEnergy.scaled(inputData = null_data,dcv_nfold = nfolds,
                                        targetFeats = targetRatios,
                                        optimizationMetric = optMetric,
                                        nreps_energy = 1000)

      ## Compute Null Metrics
      nullSubset = sep.null$finalSubset
      spEPerf.null = featureSlectionPerformance(tbl = nullSubset,Method_Name = "spE",plot_ = F)
      spEPerf.null.df = spEPerf.null$performance

      #store results
      nullPerf = rbind(nullPerf,spEPerf.null.df)

      if(optMetric=="combinedF"){
        null_tstat = spEPerf.null.df$combinedF
      }else{
        null_tstat = spEPerf.null.df$EnergyF
      }

      testStat.df = rbind(testStat.df,
                          data.frame(Type = "null",
                                     Seed = seed_,
                                     fname,rep = r,
                                     num_ratios = ncol(nullSubset),
                                     tstat =null_tstat,
                                     optMetric)
                          )


      ## set output path here
      pp = paste0("Results/Case_Study/","NewExp_Seed",seed_,"__",fname,".csv")
      write_csv(x = testStat.df,pp)
    }))

    message("selPermEnergy ",r," of ",permRep)

  }
})

