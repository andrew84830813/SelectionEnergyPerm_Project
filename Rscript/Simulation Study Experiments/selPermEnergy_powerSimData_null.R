rm(list=ls())
gc()

### Load Required Packages  ####
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


## Notes
## update path for your machine


######===========================================================-
### Read External Arguments
######===========================================================-
args = c(1,2,150,2)
args = commandArgs(trailingOnly = TRUE)
seed_ = as.numeric(args[1]) # random seed selection
sampleDesign = as.numeric(args[2]); design = if_else(sampleDesign==1,"Balanced","Unbalanced")
dims = as.numeric(args[3])
scenario = as.numeric(args[4])
######===========================================================-



###----------------------------------------*
## Power  Analysis Simulation Scenarios
###----------------------------------------*
if(sampleDesign==1){
  g1 = 40
  g2 = 40
}else{
  g1 = 20
  g2 = 60
}

set.seed(seed_)
switch(scenario,
       {
         ## Scenario-1 = difference in varaince via dirichlet distrubition
         dat = scenario1(seed = seed_,n1 = g1+g2,n2 = g2,dms_ = dims);
         dat = dat[dat$Status=="S1",]
         dat$Status = sample(c(rep("S1",g1),rep("S2",g2)))
         fname = "scenario1";
         path_ = "/nas/longleaf/home/andrew84/selEnergyPerm2/Results/powerSimData_Null/s1/"
       },

       {
         ## Scenario-2 = single feature difference among sparse noisy signal
         dat = scenario2(seed = seed_,n1 = g1+g2,n2 = g2,dms_ = dims);
         dat = dat[dat$Status=="S1",]
         dat$Status = sample(c(rep("S1",g1),rep("S2",g2)))
         fname = "scenario2";
         path_ = "/nas/longleaf/home/andrew84/selEnergyPerm2/Results/powerSimData_Null/s2/"
       },

       {
         ## Scenario-3 = Small Signal and Difference in Sparsity
         dat = scenario3(seed = seed_,n1 = g1+g2,n2 = g2,dms_ = dims);
         dat = dat[dat$Status=="S1",]
         dat$Status = sample(c(rep("S1",g1),rep("S2",g2)))
         fname = "scenario3";
         path_ = "/nas/longleaf/home/andrew84/selEnergyPerm2/Results/powerSimData_Null/s3/"
       },

       {
         ## Scenario-4 = Difference in Mean same covariance additive log normal
         dat = scenario4(seed = seed_,n1 = g1+g2,n2 = g2,dms_ = dims);
         dat = dat[dat$Status=="S1",]
         dat$Status = sample(c(rep("S1",g1),rep("S2",g2)))
         fname = "scenario4";
         path_ = "/nas/longleaf/home/andrew84/selEnergyPerm2/Results/powerSimData_Null/s4/"
       },

       {
         ## Scenario-5 = Difference in Mean same covariance additive log normal
         dat = scenario5(seed = seed_,n1 = g1+g2,n2 = g2,dms_ = dims);
         dat = dat[dat$Status=="S1",]
         dat$Status = sample(c(rep("S1",g1),rep("S2",g2)))
         fname = "scenario5";
         path_ = "/nas/longleaf/home/andrew84/selEnergyPerm2/Results/powerSimData_Null/s5/"
       }


)


message("seed",seed_,"-",design,"-dims",dims,"-",fname)


## Sparisty process
procData = processCompData(dat,minPrevalence = .9)
dat = procData$processedData
impFact = procData$impFactor
minorityClass = procData$minClss
majorityClass = procData$majClass
impClass = minorityClass
Scenario = fname
message("sparisty process Complete")


## Pre Process -  Removbe col with all 0's
y = dat[,-1]
bool = colSums(y)==0
y = y[,!bool]
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact))
message("preProcess Complete")

## Compute Ratios
dat.f = data.frame(Status = factor(dat[,1]), fastImputeZeroes( compositions::clo(dat[,-1]) , impFactor  = impFact ) )
lrs = calcLogRatio(dat.f)
message("logratios computed")




######===========================================================-
### Distance Matrix
######===========================================================-
lrs = lrs %>%
  arrange(desc(Status))
levs = data.frame(Labels = unique(lrs[,1]))
labels = lrs[,1]
tb = data.frame((table(lrs[,1])))
colnames(tb)[1] = "Labels"
sz = left_join(levs,tb)
#Distance Matrix
distMatTime = system.time({
  d.compdat = parallelDist::parDist(as.matrix(lrs[,-1]),method = "euclidean")
})
message("distance matrix computed")


## Set Common Permutations
permRep = 150
seeds = foreach(x = 1:permRep,.combine = cbind)%do%{
  set.seed(x)
  sample(1:nrow(dat))
}


##-------------------------*
### ANOSIM
##-------------------------*
tt = system.time({
  ano = vegan::anosim(x = d.compdat,grouping = lrs[,1],permutations = 0)$statistic
  null.Stats = foreach (x =1:permRep,.combine = c) %do%{
    message("ANOSIM ",x," of ",permRep)
    ano.perm = vegan::anosim(x = d.compdat,grouping = labels[seeds[,x]],permutations = 0)$statistic
  }
})
p_ = (sum(ano < null.Stats ,na.rm = T) + 1) / (permRep+1)
anosim_aitch = data.frame(numRatios = ncol(lrs[,-1]),compTime = distMatTime[3]+tt[3] ,
                          Statistic = "R",
                          Value = ano,
                          Seed = seed_,
                          permNumber = permRep,
                          type = "Empirical", test = "ANOSIM",
                          p = p_)
message("anosim complete")



##-------------------------*
### Permanova
##-------------------------*
tt = system.time({
  f.PERMANOVA_aitch = adonis2(d.compdat ~ Type,
                              data =  data.frame(Type = lrs[,1]),
                              permutations = 0)$F[1]
  null.Stats = foreach (x =1:permRep,.combine = c) %do%{
    message("PERMANOVA ",x," of ",permRep)
    f.PERMANOVA_aitch = adonis2(d.compdat ~ Type,
                                data =  data.frame(Type = labels[seeds[,x]]),
                                permutations = 0)$F[1]
  }
})

p_ = (sum(f.PERMANOVA_aitch < null.Stats ,na.rm = T) + 1) / (permRep+1)
permanova_aitch = data.frame(numRatios = ncol(lrs[,-1]),compTime = distMatTime[3]+tt[3] ,
                             Statistic = "permanovaF",
                             Value = f.PERMANOVA_aitch,
                             Seed = seed_,
                             permNumber = permRep,
                             type = "Empirical", test = "PERMANOVA",
                             p = p_)
message("permanova complete")


##-------------------------*
### Energy Test
##-------------------------*
tt = system.time({
  energy_aitch =  energy::eqdist.etest(x = d.compdat,sizes = sz$Freq,distance = T,R = 0)$statistic
  null.Stats = foreach (x =1:permRep,.combine = c) %do%{
    message("Energy test ",x," of ",permRep)
    ph = data.frame(Labels = labels[seeds[,x]] , lrs[,-1])
    ph = ph %>%
      arrange(desc(Labels))
    dperm = parallelDist::parDist(as.matrix(ph[,-1]),method = "euclidean")
    energy::eqdist.etest(x = dperm,sizes = sz$Freq,distance = T,R = 0)$statistic
  }
})
p_ = (sum(energy_aitch < null.Stats ,na.rm = T) + 1) / (permRep+1)
energy_aitch = data.frame(numRatios = ncol(lrs[,-1]),compTime = distMatTime[3]+tt[3],
                          Statistic = "estat",
                          Value = as.numeric(energy_aitch),
                          Seed = seed_,
                          permNumber = permRep,
                          type = "Empirical", test = "Energy_Aitch_ALL",
                          p = p_)
message("energy test complete")


##-------------------------*
### Dipsersion Test
##-------------------------*
tt = system.time({
  mod = vegan::betadisper(d.compdat,group = lrs[,1])
  md = as.numeric(permutest(mod,permutations = 1)$statistic)
  null.Stats = foreach (x =1:permRep,.combine = c) %do%{
    message("permDisp2 ",x," of ",permRep)
    mod = vegan::betadisper(d.compdat,group = labels[seeds[,x]] )
    as.numeric(permutest(mod,permutations = 1)$statistic)
  }
})
p_ = (sum(md < null.Stats ,na.rm = T) + 1) / (permRep+1)
permDisp_Aitch = data.frame(numRatios = ncol(lrs[,-1]),compTime = distMatTime[3]+tt[3],
                            Statistic = "permDispF",
                            Value = md,
                            Seed = seed_,
                            permNumber = permRep,
                            type = "Empirical", test = "PermDisp2",
                            p = p_
)
message("dispersion test complete")



##-------------------------*
### Selperm ENergy
##-------------------------*
tt=system.time({

  ##-------------------------------------------*
  ## Run selection Energy Routine ####
  ##-------------------------------------------*
  sep.true = selectionEnergy.scaled(inputData = dat,nreps_energy = 5e3)
  ## Get final features ####
  features.df = sep.true$finalSubset
  ## Retiecve DCV ####
  lrs_dcv = sep.true$retainedDCV
  rxDCV = sep.true$rawDCV[,2:6]
  ## compute Metric on Final Subet ####
  spEPerf = featureSlectionPerformance(tbl = features.df,Method_Name = "spE",cvRepeats = 1,plot_ = F,nreps_energy = 5e4)
  spEPerf.df = spEPerf$performance
  targetRatios = spEPerf.df$NumRatios
  optMetric = sep.true$optimization_Metric
  if(optMetric=="combinedF"){
    testStat = spEPerf.df$combinedF
  }else{
    testStat = spEPerf.df$EnergyF
  }

  ##-------------------------------------------*
  ## Null Results ####
  ##-------------------------------------------*
  nullPerf = data.frame()
  nreps = 10
  testStat.null = c()
  nreps = permRep
  system.time({
    for(r in 1:nreps){

      suppressMessages(suppressWarnings({

        ## Run selection Energy Routine on permuted data
        null_data = dat
        null_data[,1] = null_data[seeds[,r],1]
        sep.null = selectionEnergy.scaled(inputData = null_data,targetFeats = targetRatios,optimizationMetric = optMetric,nreps_energy = 5e3)

        ## Compute Null Metrics
        nullSubset = sep.null$finalSubset
        spEPerf.null = featureSlectionPerformance(tbl = nullSubset,Method_Name = "spE",cvRepeats = 1,plot_ = F,nreps_energy = 5e4)
        spEPerf.null.df = spEPerf.null$performance

        if(optMetric=="combinedF"){
          testStat.null[r] = spEPerf.null.df$combinedF
        }else{
          testStat.null[r] = spEPerf.null.df$EnergyF
        }


        #store results
        nullPerf = rbind(nullPerf,spEPerf.null.df)

      }))

      message("selPermEnergy ",r," of ",permRep)

    }
  })

})
## compute SPE pvals
trueMetrics = spEPerf.df[,5:(ncol(spEPerf.df)-1)]
trueMetrics$testStat = testStat
nullMetrics = nullPerf[,5:(ncol(nullPerf)-1)]
nullMetrics$testStat = testStat.null
discard = c("model_AUC","mds_rawEnergy")
pvals = sapply(1:ncol(trueMetrics), function(x) ( sum(trueMetrics[1,x] < nullMetrics[,x] ,na.rm = T) + 1) / (nrow(nullMetrics)+1))
stats_ = c(colnames(trueMetrics))
statValues = c(as.numeric(trueMetrics[1,]))
spp.df = data.frame(numRatios = ncol(features.df[,-1]),compTime = tt[3],
                    Statistic = stats_,
                    Value = statValues,
                    Seed = seed_,
                    permNumber = permRep,
                    type = "Empirical", test = "selProPerm",
                    p = pvals)
spp.df = spp.df[!spp.df$Statistic%in%discard,]
spp.df=  data.frame(spp.df)

######===========================================================-
# Comnined Output
######===========================================================-
message("Bind all Pvals")
out = rbind.data.frame(anosim_aitch,permanova_aitch,energy_aitch,permDisp_Aitch,spp.df)
out$scenario = fname
out$sampleDesign = design
out$dimensions = dims


### Write Data
message("Write Data")


setwd(path_)
csv_name = paste0("seed_",seed_,"-",design,"-dims_",dims,"-",fname,".csv")
write_csv(x = out,path = csv_name)
