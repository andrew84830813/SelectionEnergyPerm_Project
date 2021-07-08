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




out.df = foreach(seed_ = 1:100,.combine = rbind,.packages = c("Matrix","data.table","reshape2","vegan","tidyverse","compositions"))%dopar%{
  out = data.frame()

  for(dims in c(50,150,250)){
    for(sampleDesign in 1:2){
      for(scenario in 1:4){


        if(sampleDesign==1){
          g1 = 40
          g2 = 40
        }else{
          g1 = 20
          g2 = 60
        }

        design = if_else(sampleDesign==1,"Balanced","Unbalanced")

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

        permanova_aitch$scenario = fname
        permanova_aitch$sampleDesign = design
        permanova_aitch$dimensions = dims

        out = rbind(out,permanova_aitch)



        #write_csv(x = out,path = "Results/permanovaCorrection_null.csv")



      }
    }
  }

  out
}





### Write Data
message("Write Data")
write_csv(x = out.df,path = "Results/permanovaCorrection_null.csv")


