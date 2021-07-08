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
library(ggridges)


clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)


md = read_csv("Data/hmp2_metadata (5).csv")
md = data.frame(md)

## ihmp data
data = data.frame(read_tsv("Data/hmp2_mgx_taxonomy.tsv"))
data = separate(data,col = 1,into = c("K","P","C","O","F","Genus","Species"),sep = "\\|")
## remove ungroups in row 1
data = data[,-6:-1]
data = data %>%
  filter(!is.na(Species)) %>%
  group_by(Species) %>%
  summarise_all(.funs = sum)
## gather and spread by features samples
data = gather(data,key = 'External.ID',value = "Counts",2:ncol(data))
data = spread(data,"Species","Counts")
id = str_split(data$External.ID,pattern = "_profile",n = 2,simplify = T)[,1]
data$External.ID = id
colnames(data) = str_replace_all(colnames(data),pattern = "\\[",replacement = "")
colnames(data) = str_replace_all(colnames(data),pattern = "\\]",replacement = "")
colnames(data) = str_replace_all(colnames(data),pattern = "s__",replacement = "")




## IBD vs CD
samps = md %>%
  filter(!is.na(fecalcal),data_type=="metagenomics") %>%
  filter(diagnosis %in% c("UC","CD")) %>%
  dplyr::select(External.ID,Participant.ID,diagnosis,visit_num,fecalcal,SIBDQ.Score,SES.CD.Score,hbi,sccai,fecalcal_ng_ml)
labels = samps$diagnosis
dat = left_join(samps,data)
rownames(dat) = samps$External.ID
dat = data.frame(Status = labels,dat[,-ncol(samps):-1])



# ALL ----------------------------------------------------------------------
## metada for fecalcalprotecin
samps = md %>%
  filter(!is.na(fecalcal),data_type=="metagenomics") %>%
  #filter(diagnosis == "CD") %>%
  filter(fecalcal<50 | fecalcal>120) %>%
  dplyr::select(External.ID,Participant.ID,diagnosis,visit_num,fecalcal,SIBDQ.Score,SES.CD.Score,hbi,sccai,fecalcal_ng_ml) %>%
  mutate(level = if_else(fecalcal>120,"abnormal","normal"))
# mutate(hbi_decode = if_else(hbi>4,"active","remission")) %>%
# filter(!is.na(hbi_decode))
labels = samps$hbi_decode
labels = samps$level
dat = left_join(samps,data)
rownames(dat) = samps$External.ID
dat = data.frame(Status = labels,dat[,-ncol(samps):-1])
## PRe Process Data
dat = data.frame(dat)
table(dat[,1])
maxSparisty = .9
procData = processCompData(dat,minPrevalence = maxSparisty)
dat = procData$processedData
impFact = procData$impFactor
minorityClass = procData$minClss
majorityClass = procData$majClass
#impClass = "pos"
## Pre Process -  Removbe col with all 0's
y = dat[,-1]
bool = colSums(y)==0
y = y[,!bool]
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact))
samps = left_join(data.frame(External.ID=rownames(dat)),samps)
dx_visitnum = paste0(samps$diagnosis,samps$visit_num)
blocks = factor(dx_visitnum)
samps$perm_blocks = blocks
message("preProcess Complete")
write_csv(dat,"Output/fecalcalprotecin_ALL.csv")
write_csv(samps,"Output/fecalcalprotecin_ALL_metadata.csv")









# CD ----------------------------------------------------------------------
## metada for fecalcalprotecin
samps = md %>%
  filter(!is.na(fecalcal),data_type=="metagenomics") %>%
  filter(diagnosis == "CD") %>%
  filter(fecalcal<50 | fecalcal>120) %>%
  dplyr::select(External.ID,Participant.ID,diagnosis,visit_num,fecalcal,SIBDQ.Score,SES.CD.Score,hbi,sccai,fecalcal_ng_ml) %>%
  mutate(level = if_else(fecalcal>120,"abnormal","normal"))
  # mutate(hbi_decode = if_else(hbi>4,"active","remission")) %>%
  # filter(!is.na(hbi_decode))
labels = samps$hbi_decode
labels = samps$level
dat = left_join(samps,data)
rownames(dat) = samps$External.ID
dat = data.frame(Status = labels,dat[,-ncol(samps):-1])
## PRe Process Data
dat = data.frame(dat)
table(dat[,1])
maxSparisty = .9
procData = processCompData(dat,minPrevalence = maxSparisty)
dat = procData$processedData
impFact = procData$impFactor
minorityClass = procData$minClss
majorityClass = procData$majClass
#impClass = "pos"
## Pre Process -  Removbe col with all 0's
y = dat[,-1]
bool = colSums(y)==0
y = y[,!bool]
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact))
samps = left_join(data.frame(External.ID=rownames(dat)),samps)
message("preProcess Complete")
## Compute Ratios
dat.f = data.frame(Status = factor(dat[,1]), fastImputeZeroes( compositions::clo(dat[,-1]) , impFactor  = impFact ) )
lrs = calcLogRatio(dat.f)
message("logratios computed")
write_csv(dat,"Output/fecalcalprotecin_CD.csv")
write_csv(samps,"Output/fecalcalprotecin_CD_metadata.csv")






# UC ----------------------------------------------------------------------
## metada for fecalcalprotecin
samps = md %>%
  filter(!is.na(fecalcal),data_type=="metagenomics") %>%
  filter(diagnosis == "UC") %>%
  filter(fecalcal<50 | fecalcal>120) %>%
  dplyr::select(External.ID,Participant.ID,diagnosis,visit_num,fecalcal,SIBDQ.Score,SES.CD.Score,hbi,sccai,fecalcal_ng_ml) %>%
  mutate(level = if_else(fecalcal>120,"abnormal","normal"))
  # mutate(hbi_decode = if_else(sccai>2,"active","remission")) %>%
  # filter(!is.na(hbi_decode))

labels = samps$hbi_decode
labels = samps$level
dat = left_join(samps,data)
rownames(dat) = samps$External.ID
dat = data.frame(Status = labels,dat[,-ncol(samps):-1])
## PRe Process Data
dat = data.frame(dat)
table(dat[,1])
maxSparisty = .9
procData = processCompData(dat,minPrevalence = maxSparisty)
dat = procData$processedData
impFact = procData$impFactor
minorityClass = procData$minClss
majorityClass = procData$majClass
#impClass = "pos"
## Pre Process -  Removbe col with all 0's
y = dat[,-1]
bool = colSums(y)==0
y = y[,!bool]
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact))
samps = left_join(data.frame(External.ID=rownames(dat)),samps)
message("preProcess Complete")
## Compute Ratios
dat.f = data.frame(Status = factor(dat[,1]), fastImputeZeroes( compositions::clo(dat[,-1]) , impFactor  = impFact ) )
lrs = calcLogRatio(dat.f)
message("logratios computed")
write_csv(dat,"Output/fecalcalprotecin_UC.csv")
write_csv(samps,"Output/fecalcalprotecin_UC_metadata.csv")




##-------------------------*
### Selperm ENergy
##-------------------------*
dcv_folds = 5
px = 500
sep.true = selectionEnergy.scaled(inputData = dat,nreps_energy = 5e3,
                                  dcv_nfold = dcv_folds,patience = px,
                                  optimizationMetric = "combinedF")
x = sep.true$retainedDCV
## Get final features ####
feature.df = sep.true$finalSubset
## Retiecve DCV ####
lrs_dcv = sep.true$retainedDCV
rxDCV = sep.true$rawDCV[,2:6]
## compute Metric on Final Subet ####
spEPerf = selEnergyPermR::featureSlectionPerformance(tbl = feature.df,Method_Name = "spE",
                                                     plot_ = T,nreps_energy = 1e6)
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
## permutation design
h1 <- how(blocks = blocks)
## Set Common Permutations
permRep = 150
seeds = shuffleSet(nrow(dat),nset = permRep,control = h1)
seeds = data.frame(t(seeds))



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
      sep.null = selectionEnergy.scaled(inputData = null_data,dcv_nfold = dcv_folds,
                                        targetFeats = targetRatios,
                                        optimizationMetric = optMetric,
                                        nreps_energy = 5e3)

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

      null.testStat = rbind(null.testStat,data.frame(Sd = r,null_tstat,optMetric))


    }))

    message("selPermEnergy ",r," of ",permRep)

  }
})

