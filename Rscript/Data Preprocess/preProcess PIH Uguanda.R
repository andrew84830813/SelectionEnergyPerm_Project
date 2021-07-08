
## Required Packages
library(dplyr)
library(readr)
library(tidyr)
library(stringr)




#### Ugaunda
metaData = data.frame(read_csv(file = "Data/SamplesByMetadata_otuDADA2_PIH_Uganda_RSRC_Characteristics.csv"))
otuData = read_csv(file = "Data/SamplesByMetadata_otuDADA2_PIH_Uganda_RSRC_TaxaRelativeAbundance.csv")[,-10]

#Merge Taxa ID's
otuData = otuData %>%
  unite("taxa", Phylum:Genus, remove = T,sep = "|")
otuData = data.frame(otuData)
otuData = otuData[,-c(2,3,5)]
otuData = otuData %>%
  group_by(Sample.ID,taxa) %>%
  summarise(Abundance = sum(Absolute.Abundance))
#final ASV Table
otuData = spread(otuData,key = "taxa",value = "Abundance",fill = 0)


#Metadata
metaData = metaData[,-c(2,5,6)]
metaData = spread(metaData,key = "Property",value = "Value")
sampleLabels = metaData[,c(1,10)]


#Merge Taxa and Labels
df = left_join(sampleLabels,otuData)
df = df[,-1]
colnames(df)[1] = "Status"
df$Status = if_else(df$Status=="Postinfectious hydrocephalus (PIH)","pos","neg")
write_csv(df,path = "Data/Uguanda_PIH.csv")



