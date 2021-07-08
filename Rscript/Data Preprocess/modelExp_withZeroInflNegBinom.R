library(tidyr)
library(dplyr)
library(stringr)
library(readr)



## Start Cluster
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)



dat.16s = read_csv("Data/16S_microbiomeHD_goodrich_Healthy.csv")
rs = rowSums(dat.16s[,-1])
dat.wgs = read_csv("Data/WGS_curatedMetagenome_ZeeviD2015_controls.csv")
rs = rowSums(dat.wgs[,-1])

## Pre Process -  Removbe col with all 0's
procData = processCompData(dat.16s,minPrevalence = .97)
dat.16s = procData$processedData
y = dat.16s[,-1]
bool = colSums(y)==0
y = y[,!bool]
y.16s = y

## Pre Process -  Removbe col with all 0's
procData = processCompData(dat.wgs,minPrevalence = .97)
dat.wgs = procData$processedData
y = dat.wgs[,-1]
bool = colSums(y)==0
y = y[,!bool]
y.wgs = y

## Fit Model
ft.wgs =zinbwave::zinbFit(t(y.wgs),zeroinflation = T)
saveRDS(ft.wgs,file = "Output/wgsModel.RDS")
mdwgs = readRDS("Output/wgsModel.RDS")


ft.16s =zinbwave::zinbFit(t(y.16s),zeroinflation = T)
md16s = readRDS("Output/16sModel.RDS")
saveRDS(ft.16s,file = "Output/16sModel.RDS")



## Example Simulation

## Simulate 16S
mdwgs = readRDS("Output/16sModel.RDS");
set.seed(08272008);
dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts))
colnames(dat) = colnames(dat.16s[,-1])
dat = sample_n( data.frame(Data = "Simulated",dat), size = 428,replace = F)
dat.16s = sample_n ( data.frame(Data = "True",dat.16s[,-1]), size = 428,replace = F)
cb = rbind(dat,dat.16s)
colnames(cb) = c("Status",paste0("V",2:ncol(cb)))
cbx = fastImputeZeroes(cb[,-1])
#lrs = data.frame(calcLogRatio( fastImputeZeroes(cb)))
pc = Rtsne::Rtsne( clo(cbx) )
pc.df = data.frame(Data = cb[,1],pc$Y)

tiff(filename = "Figures/dataSimulation_tsne_16s.tiff",width = 4,height = 3,units = "in",res = 300)
ggplot(pc.df,aes(X1,X2))+
  geom_point(aes(col = Data),alpha = .5  )+
  ggthemes::geom_rangeframe()+
  ggthemes::theme_gdocs()+
  theme(legend.position = c(.85,.9))+
  theme(
    axis.line = element_line(),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0,face = "bold",size = 6,hjust = 0),
    axis.ticks = element_line(colour = "black",size = 1),
    #axis.title.y = element_blank(),
    axis.title = element_text(size = 8,face = "bold"),
   # legend.position = "top",
    legend.text = element_text(size = 8),
    legend.margin = margin(0,0,0,0,unit = "cm"),
    legend.box.spacing = unit(0,units = "in"),
    legend.box.margin = margin(0,0,0.1,0,unit = "cm"),
    legend.key.height = unit(.1,units = "in"),
    axis.text =  element_text(size = 8),
    #axis.text.y = element_blank(),
    legend.title =element_text(size = 8),
    plot.caption = element_text(size = 8,face = "italic"))
dev.off()

ct.sim = colMeans(clo(dat[,-1]))
sim.df = data.frame(Data = "Sim",Ra = ct.sim)
ct.true = colMeans(clo(dat.16s[,-1]))
true.df = data.frame(Data = "True",Ra = ct.true)
al = rbind(sim.df,true.df)

tiff(filename = "Figures/dataSimulation_histo_16s.tiff",width = 4,height = 3,units = "in",res = 300)
ggplot(al,aes(Ra))+
  geom_histogram(col = "black",fill = "skyblue")+
  ggthemes::theme_clean()+
  xlab("Mean Relative Abundance")+
  facet_wrap(Data~.,nrow = 2)+
  theme(
    axis.line = element_line(),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0,face = "bold",size = 6,hjust = 0),
    axis.ticks = element_line(colour = "black",size = 1),
    #axis.title.y = element_blank(),
    axis.title = element_text(size = 8,face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 6),
    legend.margin = margin(0,0,0,0,unit = "cm"),
    legend.box.spacing = unit(0,units = "in"),
    legend.box.margin = margin(0,0,0.1,0,unit = "cm"),
    legend.key.height = unit(.1,units = "in"),
    axis.text =  element_text(size = 8),
    #axis.text.y = element_blank(),
    legend.title =element_blank(),
    plot.caption = element_text(size = 8,face = "italic"))
dev.off()


ct.sim = colMeans(clo(dat[,-1]))
ct.true = colMeans(clo(dat.16s[,-1]))
hist(ct.sim)
hist(ct.true)
ks.test(x = ct.sim,y = ct.true)
cor.test(ct.sim,ct.true)
plot(ct.sim,ct.true)





## Simulate WGS
mdwgs = readRDS("Output/wgsModel.RDS");
set.seed(08272008);
dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts))
colnames(dat) = colnames(dat.wgs[,-1])
dat = sample_n( data.frame(Data = "Simulated",dat), size = 900,replace = F)
dat.wgs = sample_n ( data.frame(Data = "True",dat.wgs[,-1]), size = 900,replace = F)
cb = rbind(dat,dat.wgs)
cbx = fastImputeZeroes(cb[,-1])
#lrs = data.frame(calcLogRatio( fastImputeZeroes(cb)))
pc = Rtsne::Rtsne( clo(cbx) )
pc.df = data.frame(Data = cb[,1],pc$Y)
setwd("C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/selectionEnergyPermutation/")

tiff(filename = "Figures/dataSimulation_tsne_wgs.tiff",width = 4,height = 3,units = "in",res = 300)
ggplot(pc.df,aes(X1,X2))+
  geom_point(aes(col = Data),alpha = .5)+
  ggthemes::geom_rangeframe()+
  ggthemes::theme_gdocs()+
  theme(legend.position = c(.15,.8))+
  theme(
    axis.line = element_line(),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0,face = "bold",size = 6,hjust = 0),
    axis.ticks = element_line(colour = "black",size = 1),
    #axis.title.y = element_blank(),
    axis.title = element_text(size = 8,face = "bold"),
    legend.text = element_text(size = 8),
    legend.margin = margin(0,0,0,0,unit = "cm"),
    legend.box.spacing = unit(0,units = "in"),
    legend.box.margin = margin(0,0,0.1,0,unit = "cm"),
    legend.key.height = unit(.1,units = "in"),
    axis.text =  element_text(size = 8),
    #axis.text.y = element_blank(),
    legend.title =element_text(size = 8),
    plot.caption = element_text(size = 8,face = "italic"))
dev.off()

ct.sim = colMeans(clo(dat[,-1]))
sim.df = data.frame(Data = "Sim",Ra = ct.sim)
ct.true = colMeans(clo(dat.wgs[,-1]))
true.df = data.frame(Data = "True",Ra = ct.true)
al = rbind(sim.df,true.df)

tiff(filename = "Figures/dataSimulation_histo_wgs.tiff",width = 4,height = 3,units = "in",res = 300)
ggplot(al,aes(Ra))+
  geom_histogram(col = "black",fill = "skyblue")+
  ggthemes::theme_clean()+
  xlab("Mean Relative Abundance")+
  facet_wrap(Data~.,nrow = 2)+
  theme(
    axis.line = element_line(),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0,face = "bold",size = 6,hjust = 0),
    axis.ticks = element_line(colour = "black",size = 1),
    #axis.title.y = element_blank(),
    axis.title = element_text(size = 8,face = "bold"),
    legend.text = element_text(size = 8),
    legend.margin = margin(0,0,0,0,unit = "cm"),
    legend.box.spacing = unit(0,units = "in"),
    legend.box.margin = margin(0,0,0.1,0,unit = "cm"),
    legend.key.height = unit(.1,units = "in"),
    axis.text =  element_text(size = 8),
    #axis.text.y = element_blank(),
    legend.title =element_text(size = 8),
    plot.caption = element_text(size = 8,face = "italic"))
dev.off()


hist(ct.sim,breaks = 10)
hist(ct.true,breaks = 10)
ks.test(x = ct.sim,y = ct.true)
cor.test(ct.sim,ct.true)
