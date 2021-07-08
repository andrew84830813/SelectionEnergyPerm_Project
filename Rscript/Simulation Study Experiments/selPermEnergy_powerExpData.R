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


## Start Local Cluster
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)



## Notes
## update path for your machine


######===========================================================-
### Read External Arguments
######===========================================================-
args = c(2,1,1,1,2)
args = commandArgs(trailingOnly = TRUE)
seed_ = as.numeric(args[1]) # random seed selection
sampleDesign = as.numeric(args[2]); design = if_else(sampleDesign==1,"Balanced","Unbalanced")
shiftParm = as.numeric(args[3])
scenario = as.numeric(args[4])
label_parm = 2; trueDistr = if_else(label_parm==1,"same","different")
######===========================================================-



###----------------------------------------*
## Power  Analysis Simulation Scenarios
###----------------------------------------*
sparsePercent = 0.85
if(sampleDesign==1){
  g1 = 40
  g2 = 40
}else{
  g1 = 20
  g2 = 60
}




switch(scenario,
       {
         ## Scenario-1 = Sim from 16S Data with mean shift
         mdwgs = readRDS("Output/16sModel.RDS");
         set.seed(seed_);
         dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
         dat = sample_n(dat,size = g1+g2,replace = F);
         labels = sample(c(rep("S1",g1),rep("S2",g2)));
         dat = data.frame(Status = labels,dat);
         fname = "16S_meanShift";
         path_ = "/nas/longleaf/home/andrew84/selEnergyPerm2/Results/powerExpData/s1/";
         #process shift
         procData = processCompData(dat,minPrevalence = sparsePercent);
         dat = procData$processedData;
         impFact = procData$impFactor;
         y = dat[,-1];
         bool = colSums(y)==0;
         y = y[,!bool];
         dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));

         if(label_parm==2){
           parm = seq(.5,.95,length.out = 4)[shiftParm];
           shiftPer = seq(1.25,1.25,length.out = 4)[shiftParm]
           dat = simFromExpData.largeMeanShft(raMatrix = dat[,-1],
                                              n1 = g1,n2 = g2,
                                              featureShiftPercent =  shiftPer,
                                              perFixedFeatures = parm)
         }

       },

       {
         ## Scenario-2 = Sim from 16S Data with cov shift
         mdwgs = readRDS("Output/16sModel.RDS");
         set.seed(seed_);
         dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
         dat = sample_n(dat,size = g1+g2,replace = F);
         labels = sample(c(rep("S1",g1),rep("S2",g2)));
         dat = data.frame(Status = labels,dat);
         fname = "16S_covShift";
         path_ = "/nas/longleaf/home/andrew84/selEnergyPerm2/Results/powerExpData/s2/";
         procData = processCompData(dat,minPrevalence = sparsePercent);
         dat = procData$processedData;
         impFact = procData$impFactor;
         y = dat[,-1];
         bool = colSums(y)==0;
         y = y[,!bool];
         dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));

         if(label_parm==2){

           cc = seq(0.1,4,length.out = 4)
           sp = seq(1.2,1.1,length.out = 4)
           dat = simFromExpData.covarianceShift(raMatrix = dat[,-1],
                                                n1 = g1,n2 = g2,
                                                maxCov =cc[shiftParm],
                                                shiftPercent = sp[shiftParm],
                                                perFixedFeatures = .9)


         }

       },

       {
         ## Scenario-3 = Sim from WGS Data with mean shift
         mdwgs = readRDS("Output/wgsModel.RDS");
         set.seed(seed_);
         dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
         dat = sample_n(dat,size = g1+g2);
         labels = sample(c(rep("S1",g1),rep("S2",g2)));
         dat = data.frame(Status = labels,dat);
         fname = "WGS_meanShift";
         path_ = "/nas/longleaf/home/andrew84/selEnergyPerm2/Results/powerExpData/s3_/";
         procData = processCompData(dat,minPrevalence = sparsePercent);
         dat = procData$processedData;
         impFact = procData$impFactor;
         y = dat[,-1];
         bool = colSums(y)==0;
         y = y[,!bool];
         dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));

         if(label_parm==2){
           parm = seq(.5,.95,length.out = 4)[shiftParm];
           shiftPer = seq(1.25,1.25,length.out = 4)[shiftParm]
           dat = simFromExpData.largeMeanShft(raMatrix = dat[,-1],
                                              n1 = g1,n2 = g2,
                                              featureShiftPercent =  shiftPer,
                                              perFixedFeatures = parm)
         }
       },


       {
         ## Scenario-4 = Sim from WGS Data with cov shift
         mdwgs = readRDS("Output/wgsModel.RDS");
         set.seed(seed_);
         dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
         dat = sample_n(dat,size = g1+g2);
         labels = sample(c(rep("S1",g1),rep("S2",g2)));
         dat = data.frame(Status = labels,dat);
         fname = "WGS_covShift";
         path_ = "/nas/longleaf/home/andrew84/selEnergyPerm2/Results/powerExpData/s4/";
         procData = processCompData(dat,minPrevalence = sparsePercent);
         dat = procData$processedData;
         impFact = procData$impFactor;
         y = dat[,-1];
         bool = colSums(y)==0;
         y = y[,!bool];
         dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));

         if(label_parm==2){


           cc = seq(0.1,12,length.out = 4)
           sp = seq(1.25,1.25,length.out = 4)
           dat = simFromExpData.covarianceShift(raMatrix = dat[,-1],
                                                n1 = g1,n2 = g2,
                                                maxCov =cc[shiftParm],
                                                shiftPercent = sp[shiftParm],
                                                perFixedFeatures = .9)
         }

       }

)


#visualizeData_pca(dat)

message("seed",seed_,"-",design,"-shiftParm",shiftParm,"-",fname,label_parm)

whichCoef_i = shiftParm
## Pre - Process , Compute Ratios, and Distance Matrix
{
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

}

message("Distr Difference = ",trueDistr)
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
  trueF = adonis2(d.compdat ~ Type,
                  data =  data.frame(Type = lrs[,1]),
                  permutations = 0)$F[1]
  null.Stats = foreach (x =1:permRep,.combine = c) %do%{
    message("PERMANOVA ",x," of ",permRep)
    f.PERMANOVA_aitch = adonis2(d.compdat ~ Type,
                                data =  data.frame(Type = labels[seeds[,x]]),
                                permutations = 0)$F[1]
  }
})


p_ = (sum(trueF < null.Stats ,na.rm = T) + 1) / (permRep+1)
permanova_aitch = data.frame(numRatios = ncol(lrs[,-1]),compTime = distMatTime[3]+tt[3] ,
                             Statistic = "permanovaF",
                             Value = trueF,
                             Seed = seed_,
                             permNumber = permRep,
                             type = "Empirical", test = "PERMANOVA",
                             p = p_)
message("permanova complete")


##-------------------------*
### Energy Test
##-------------------------*
#energy::eqdist.etest(x = d.compdat,sizes = sz$Freq,distance = T,R = 100)
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
### Selperm Energy
##-------------------------*
tt=system.time({

  ##-------------------------------------------*
  ## Run selection Energy Routine ####
  ##-------------------------------------------*
  sep.true = selectionEnergy.scaled(inputData = dat)
  ## Get final features ####
  #features.df = data.frame(Status = sep.true$finalSubset[,1],sep.true$mstSubset)
  features.df = sep.true$finalSubset
  ## Retiecve DCV ####
  lrs_dcv = sep.true$retainedDCV
  rxDCV = sep.true$rawDCV[,2:6]
  ## compute Metric on Final Subet ####
  spEPerf = featureSlectionPerformance(tbl = features.df,Method_Name = "spE",cvRepeats = 1,plot_ = F)
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
        sep.null = selectionEnergy.scaled(inputData = null_data,targetFeats = targetRatios,optimizationMetric = optMetric)

        ## Compute Null Metrics
        #nullSubset = data.frame(Status = sep.null$finalSubset[,1],sep.null$mstSubset)
        nullSubset = sep.null$finalSubset
        spEPerf.null = featureSlectionPerformance(tbl = nullSubset,Method_Name = "spE",cvRepeats = 1,plot_ = F)
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
out$shiftParm = shiftParm
out$trueDistr = trueDistr


### Write Data
message("Write Data")


setwd(path_)
csv_name = paste0("seed_",seed_,"-",design,"-dims_",shiftParm,"-",trueDistr,"-",fname,".csv")
write_csv(x = out,path = csv_name)
