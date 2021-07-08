library(GUniFrac)
library(biomformat)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)



## read metadata
sampleInfo = read_tsv(file = "Data/10249_20200224-135750.txt")

##Read biom files
#nn = read_biom(biom_file = "219_otu_table.biom")
nn = read_biom(biom_file = "Data/otu_table.biom")

fnames = dir("Data/delivery_modeStudy/")

otu.df = data.frame()
si = data.frame()
for(f in fnames){
  nn = read_biom(biom_file = paste0("Data/delivery_modeStudy/",f))

  ##OTU Designation
  md = observation_metadata(nn)
  md = data.frame(otu_id = rownames(md),md)

  ##get otu datat
  otu_table = biom_data(nn)
  otu_table = as.matrix(otu_table)
  otu_table = data.frame(otu_id = rownames(otu_table),otu_table )
  otu_table = gather(otu_table,"sample_name","counts",2:ncol(otu_table))
  otu_table = left_join(otu_table,md)
  ##filter by family|genus
  otu_table = otu_table %>%
    filter(taxonomy5!="f__") %>%
    mutate(fam_genus = paste0(taxonomy5,"|",taxonomy6)) %>%
    select(sample_name,fam_genus,counts) %>%
    group_by(sample_name,fam_genus) %>%
    summarise(counts = sum(counts))

  ## merge metadata with count table
  sampleInfo1 = sampleInfo
  sampleInfo1$sample_name = paste0("X",sampleInfo1$sample_name)
  sampleInfo1 = sampleInfo1 %>%
    filter(sample_name%in%otu_table$sample_name)
  si = rbind(si,sampleInfo1)
  otu_table = left_join(sampleInfo1,otu_table)

  ## select only children with known delivery mode
  otu_table =otu_table %>%
    filter(mom_child=="C") %>%
    filter(delivery!="na")

  otu.df = rbind(otu.df,otu_table)

}

## final process otu.df
otu.df1 = otu.df %>%
  select(sample_name,delivery,host_subject_id,month,fam_genus,counts)%>%
  spread("fam_genus","counts",fill = 0)

##get first 6 months
for(m in 1:6){
  ph = otu.df1 %>%
    filter(month==m)
  write_csv(ph,file = paste0("Output/deliveryMode_month",m,".csv"))
}










## Preprocess birth data
fnames = dir("Output/")
bool = str_detect(fnames,pattern = "deliveryMode")
fnames = fnames[bool]

for(f in fnames){
  dat = read_csv(file = paste0("Output/",f))
  md = dat[,1:4]
  dat = data.frame(Status= dat$delivery,dat[,-4:-1])
  f = paste0("Output/",f)
  ## fix column names
  colnames(dat) = str_replace(colnames(dat),pattern  = "__\\.",replacement = "__");
  colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = "#");
  colnames(dat) = str_replace(colnames(dat),pattern  = "f__",replacement = "");
  colnames(dat) = str_replace(colnames(dat),pattern  = "\\.g__",replacement = ".");
  colnames(dat) = str_replace(colnames(dat),pattern  = "#g__",replacement = ".");
  colnames(dat) = str_replace(colnames(dat),pattern  = "g__",replacement = "__");
  colnames(dat) = str_replace(colnames(dat),pattern  = "\\.\\.",replacement = ".")
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
  dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact))
  write_csv(dat,file = f)
  mm = str_replace(f,"\\.csv","_metadata.csv")
  write_csv(md,file = mm)

}
