
library(tidyr)
library(dplyr)
library(stringr)
library(readr)

load("Data/DIABIMMUNE_Karelia_metadata.RData")
comb_metadata = read_csv("Data/diabimmune_comb_metadata.csv")



## Process WGS
df = read_tsv("Data/diabimmune_karelia_metaphlan_table.txt")
df = separate(df,col = 1,into = c("k","p","c","o","f","g","species","t"),sep = "\\|")
df = df %>%
  filter(!is.na(species)) %>%
  filter(is.na(t)) %>%
  dplyr::select(-k,-p,-c,-o,-f,-g,-t)
df$species = str_replace_all(df$species,"s__","")
df = df %>%
  gather("SampleID","counts",2:ncol(df)) %>%
  group_by(SampleID,species,) %>%
  summarise(counts = sum(counts)) %>%
  spread("species","counts",fill = 0)
rowSums(df[,-1])



# All data
mgx_samples = comb_metadata %>%
  filter(!is.na(gid_wgs))
samps = metadata %>%
  filter(SampleID %in% mgx_samples$sampleID) %>%
  select(subjectID,SampleID,age_at_collection,collection_month,country,allergy_milk,allergy_egg,allergy_peanut,gid_wgs) %>%
  filter(!is.na(gid_wgs)) %>%
  mutate(Status = if_else(allergy_egg+allergy_peanut+allergy_milk>0,"FA","none")) %>%
  filter(!is.na(Status)) %>%
  filter(collection_month<=24) %>%
  mutate(SampleID = gid_wgs) %>%
  select(SampleID,country,Status,collection_month,subjectID)

dat = left_join(samps,df)
label = factor(dat$Status)
dat = data.frame(Status = label,dat[,-ncol(samps):-1])
rs = rowSums(dat[,-1])
bool = !is.na(rs)
dat =dat[bool,]
samps = samps[bool,]
month_country = paste0(samps$country,samps$collection_month)

samps %>%
  group_by(Status,subjectID,country) %>%
  summarise(n = n()) %>%
  group_by(country) %>%
  summarise(n = n())


# First 6 months ----------------------------------------------------------
mgx_samples = comb_metadata %>%
  filter(!is.na(gid_wgs))
samps = metadata %>%
  filter(SampleID %in% mgx_samples$sampleID) %>%
  select(subjectID,SampleID,age_at_collection,collection_month,country,allergy_milk,allergy_egg,allergy_peanut,gid_wgs) %>%
  filter(!is.na(gid_wgs)) %>%
  mutate(Status = if_else(allergy_egg+allergy_peanut+allergy_milk>0,"FA","none")) %>%
  filter(!is.na(Status)) %>%
  filter(collection_month<=6) %>%
  mutate(SampleID = gid_wgs) %>%
  select(SampleID,country,Status,collection_month,subjectID)
dat = left_join(samps,df)
label = factor(dat$Status)
dat = data.frame(Status = label,dat[,-ncol(samps):-1])
rs = rowSums(dat[,-1])
bool = !is.na(rs)
dat =dat[bool,]
samps = samps[bool,]
month_country = paste0(samps$country,samps$collection_month)
blocks = factor(month_country)
n_distinct(blocks)
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

## our imputation
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact))

d1 = gather(dat,"taxa","count",2:ncol(dat))
d1$id =  rownames(dat)
samps$permutation_blocks = blocks
write_csv(dat,"Output/diabimmune_f6.csv")
write_csv(samps,"Output/diabimmune_f6_metadata.csv")





# 6-12 months -------------------------------------------------------------
mgx_samples = comb_metadata %>%
  filter(!is.na(gid_wgs))
samps = metadata %>%
  filter(SampleID %in% mgx_samples$sampleID) %>%
  select(subjectID,SampleID,age_at_collection,collection_month,country,allergy_milk,allergy_egg,allergy_peanut,gid_wgs) %>%
  filter(!is.na(gid_wgs)) %>%
  mutate(Status = if_else(allergy_egg+allergy_peanut+allergy_milk>0,"FA","none")) %>%
  filter(!is.na(Status)) %>%
  filter(collection_month>6 & collection_month<=12) %>%
  mutate(SampleID = gid_wgs) %>%
  select(SampleID,country,Status,collection_month,subjectID)
dat = left_join(samps,df)
rownames(dat) = samps$SampleID
label = factor(dat$Status)
dat = data.frame(Status = label,dat[,-ncol(samps):-1])
rs = rowSums(dat[,-1])
bool = !is.na(rs)
dat =dat[bool,]
samps = samps[bool,]
n_distinct(blocks)
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

zCompositions::multRepl()
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact),row.names = rownames(dat))
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact),row.names = rownames(dat))
samps = left_join(data.frame(SampleID = rownames(dat)),samps)
d2 = gather(dat,"taxa","count",2:ncol(dat))
d2$id =  rownames(dat)

# ## Down sample
# tt = table(samps$Status)
# m = dat %>%
#   filter(Status==majorityClass)
# mn = dat %>%
#   filter(Status==minorityClass)
# m = sample_n(m,size = tt[minorityClass],replace = F)
# dat = rbind(mn,m)
# samps = left_join(data.frame(SampleID = rownames(dat)),samps)
# month_country = paste0(samps$country,samps$collection_month)
# blocks = factor(month_country)
# samps$permutation_blocks = blocks


write_csv(dat,"Output/diabimmune_6-12.csv")
write_csv(samps,"Output/diabimmune_6-12_metadata.csv")



# 12-18 months ------------------------------------------------------------
mgx_samples = comb_metadata %>%
  filter(!is.na(gid_wgs))
samps = metadata %>%
  filter(SampleID %in% mgx_samples$sampleID) %>%
  select(subjectID,SampleID,age_at_collection,collection_month,country,allergy_milk,allergy_egg,allergy_peanut,gid_wgs) %>%
  filter(!is.na(gid_wgs)) %>%
  mutate(Status = if_else(allergy_egg+allergy_peanut+allergy_milk>0,"FA","none")) %>%
  filter(!is.na(Status)) %>%
  filter(collection_month>12 & collection_month<=18) %>%
  mutate(SampleID = gid_wgs) %>%
  select(SampleID,country,Status,collection_month,subjectID)
dat = left_join(samps,df)
rownames(dat) = samps$SampleID
label = factor(dat$Status)
dat = data.frame(Status = label,dat[,-ncol(samps):-1])
rs = rowSums(dat[,-1])
bool = !is.na(rs)
dat =dat[bool,]
samps = samps[bool,]
n_distinct(blocks)
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
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact),row.names = rownames(dat))
samps = left_join(data.frame(SampleID = rownames(dat)),samps)
d3 = gather(dat,"taxa","count",2:ncol(dat))
d3$id =  rownames(dat)


# ## Down sample
# tt = table(samps$Status)
# m = dat %>%
#   filter(Status==majorityClass)
# mn = dat %>%
#   filter(Status==minorityClass)
# m = sample_n(m,size = tt[minorityClass],replace = F)
# dat = rbind(mn,m)
# samps = left_join(data.frame(SampleID = rownames(dat)),samps)
# month_country = paste0(samps$country,samps$collection_month)
# blocks = factor(month_country)
# samps$permutation_blocks = blocks

write_csv(dat,"Output/diabimmune_12-18.csv")
write_csv(samps,"Output/diabimmune_12-18_metadata.csv")



# 18-24 months ------------------------------------------------------------
mgx_samples = comb_metadata %>%
  filter(!is.na(gid_wgs))
samps = metadata %>%
  filter(SampleID %in% mgx_samples$sampleID) %>%
  select(subjectID,SampleID,age_at_collection,collection_month,country,allergy_milk,allergy_egg,allergy_peanut,gid_wgs) %>%
  filter(!is.na(gid_wgs)) %>%
  mutate(Status = if_else(allergy_egg+allergy_peanut+allergy_milk>0,"FA","none")) %>%
  filter(!is.na(Status)) %>%
  filter(collection_month>18 & collection_month<=24) %>%
  mutate(SampleID = gid_wgs) %>%
  select(SampleID,country,Status,collection_month,subjectID)
dat = left_join(samps,df)
rownames(dat) = samps$SampleID
label = factor(dat$Status)
dat = data.frame(Status = label,dat[,-ncol(samps):-1])
rs = rowSums(dat[,-1])
bool = !is.na(rs)
dat =dat[bool,]
samps = samps[bool,]
n_distinct(blocks)
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
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact),row.names = rownames(dat))
samps = left_join(data.frame(SampleID = rownames(dat)),samps)
d4 = gather(dat,"taxa","count",2:ncol(dat))
d4$id =  rownames(dat)
# ## Down sample
# tt = table(samps$Status)
# m = dat %>%
#   filter(Status==majorityClass)
# mn = dat %>%
#   filter(Status==minorityClass)
# m = sample_n(m,size = tt[minorityClass],replace = F)
# dat = rbind(mn,m)
# samps = left_join(data.frame(SampleID = rownames(dat)),samps)
# month_country = paste0(samps$country,samps$collection_month)
# blocks = factor(month_country)
# samps$permutation_blocks = blocks
#


write_csv(dat,"Output/diabimmune_18-24.csv")
write_csv(samps,"Output/diabimmune_18-24_metadata.csv")






##
all = rbind(d1,d2) %>%
  rbind(d3) %>%
  rbind(d4)

all = spread(all,"taxa","count")




