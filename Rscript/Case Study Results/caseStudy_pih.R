rm(list = ls())
gc()



library(compositions)
library(data.table)
library(reshape2)
library(doParallel)
library(igraph)
library(caret)
library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(glmnet)
library(vegan)
library(PRROC)
library(diffCompVarRcpp)
library(selEnergyPermR)
library(simplexDataAugmentation)
library(rstatix)




source("Helper Functions/functions.R")


# Permutation Result ------------------------------------------------------

### read SelEnergyPErm Resutls
fnames = dir("Results/Case_Study/")
bool = str_detect(fnames,pattern = "pih")
fnames = fnames[bool]
results.df = data.frame()
for(f in fnames){
  ph = read_csv(file = paste0("Results/Case_Study/",f))
  results.df = rbind(results.df,ph)
}
t1.null = results.df %>% 
  filter(Type=="null")
t1.empi = results.df %>% 
  filter(Type!="null")
ph = rbind(t1.null,t1.empi[1,])
testStat =  unique(t1.empi$tstat)
dt = t1.null$tstat
1-ecdf(dt)(testStat)
hist(dt)
mm = mean(dt)
ss = sd(dt)
pval = 1-pnorm(testStat,mean = mm,sd = ss)
# sep chat=rt
x = dt
permRep = length(x)
x = rnorm(1000,mean = mean(x),sd(x))
t_stat =  unique(t1.empi$tstat)
tstat_p =   (sum(t1.null$tstat>t_stat)+1)/ (nrow(t1.null)+1)
xden = density(x)$x
yden = density(x)$y
xx.df = data.frame(X = xden)
x.df = data.frame(X = xden,Y = yden)
jitter.df = data.frame(X = x,jitter = rnorm(length(x),quantile(yden,probs = .15),sd = 0.001))

pdf(file = "Figures/caseStudy_pih_sepResults.pdf",width = 3 ,height = 2.5)
ggplot(x.df,aes(X,Y))+
  geom_histogram(data = jitter.df,aes(y = ..density..),colour = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(x), sd = sd(x)),col = "blue")+
  geom_point(data = jitter.df,aes(X,jitter),col = "royalblue3",alpha = .2)+
  ggthemes::theme_tufte()+
  ggthemes::geom_rangeframe() +
  geom_vline(xintercept = t_stat,lty = "dashed",col = "red",size = .75)+
  xlab("Combined F-Statistic")+
  ylab("Density")+
  labs(caption = paste0("cF = ",round(t_stat,digits = 4)," , pval = ",
                        round(tstat_p,digits = 6), " (permutations = ",permRep,")"))+
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
    plot.caption = element_text(size = 8))
dev.off()








##-----------------------------------------*
### Read Data ####
##-----------------------------------------*
## Load Data
## Scenario-3 = Sim from additive log normal
fname = "16S_Uguanda_PIH";
dat = read_csv("Output/Uguanda_PIH.csv")



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


## Aggregate to Genus
dat = data.frame(id = 1:nrow(dat),dat)
dat = gather(dat,"taxa","count", 3:ncol(dat))
dat = separate(dat,col = 3,into = c("p","c","o","family","genus","s"))
dat[dat=="null"] = NA
dat = dat %>%
  filter(!is.na(family)) %>%
  mutate(genus = if_else(is.na(genus),"g.",genus)) %>%
  mutate(family_genus = paste0(family,"_",genus)) %>%
  group_by(id,Status,family_genus) %>%
  summarise(count = sum(count)) %>%
  spread("family_genus","count") %>%
  ungroup() %>%
  select(-id)
dat = data.frame(dat)


## SelENergyPerm Impute
message("preProcess Complete")

 

## Compute Ratios
dat.f = data.frame(Status = factor(dat[,1]), fastImputeZeroes( compositions::clo(dat[,-1]) , impFactor  = impFact ) )
lrs = calcLogRatio(dat.f)
message("logratios computed")

## Get SelEnergy Signal
sep.true = selectionEnergy.scaled(inputData = dat,nreps_energy = 1000,
                                  dcv_nfold = 1,patience = 500)

x = sep.true$retainedDCV
feature.df = sep.true$finalSubset
lrs_dcv = sep.true$retainedDCV
lrs_dcv = lrs_dcv %>%
  filter(Ratio%in% colnames(feature.df[,-1]))
xx = getLogratioFromList(Ratio = colnames(feature.df[,-1]),raMatrix = dat,Class = "test")
feature.df = data.frame(Status = dat[,1],xx)
write_csv(lrs_dcv,"Output/pih_signature_dcv.csv")
write_csv(feature.df,"Output/pih_signature.csv")



##read signature
feature.df = read.csv("Output/pih_signature.csv")

# ROC Anlalysis -----------------------------------------------------------

tc1  <- trainControl(method="repeatedcv",
                     repeats = 50,
                     number = 10,
                     #search = "random",
                     #sampling = "rose",
                     #number=100,
                     seeds = NULL,
                     classProbs = TRUE,
                     savePredictions = T,
                     allowParallel = TRUE,
                     summaryFunction = caret::multiClassSummary
                     #summaryFunction = prSummary
)


pb = subset(dat,select =  "Paenibacillaceae_Paenibacillus")
glm.mdl1=train(x = data.frame(pb),
               y = dat[,1],
               probMethod = "softmax",
               method = "rf",
               #tuneGrid = expand.grid(ncomp = 1),
               metric = "AUC",
               trControl = tc1
)

pb =  glm.mdl1$results
pb$type =  "Paenibacillus alone"
# pb = glm.mdl1$resample
# pb = separate(pb,col = 15,into = c("fold","rep"),remove = F)
# pb$type =  "Paenibacillus alone"
# pb = pb %>% 
#   group_by(rep,type) %>% 
#   summarise(auc = mean(AUC))




glm.mdl1=train(x = data.frame(lrs[,-1]),
               y = lrs[,1],
               #probMethod = "softmax",
               method = "rf",
               tuneGrid = expand.grid(mtry = round(sqrt(ncol(lrs[,-1])))),
               metric = "AUC",
               trControl = tc1
)
all = glm.mdl1$results
all$type = "all logratios"


glm.mdl1 = train(x = feature.df[,-1],
               y = feature.df[,1],
               probMethod = "softmax",
               method = "rf",proximity = T,importance=T,
               tuneGrid = expand.grid(mtry = round(sqrt(ncol(feature.df[,-1])))),
               metric = "AUC",
               trControl = tc1
)
vi = varImp(glm.mdl1,scale = F)
vi = data.frame(taxa  =rownames(vi$importance),Importance = vi$importance[,1])
vi$taxa = str_replace(vi$taxa,pattern = "___"," / ")

## roc
preds = glm.mdl1$pred
signature.roc = glm.mdl1$results
signature.roc$type = "signature"

cc = rbind(all,signature.roc)
cc = rbind(cc,pb)
cc$ymin = cc$AUC - (1.96*cc$AUCSD)/sqrt(50)
cc$ymax = cc$AUC + (1.96*cc$AUCSD)/sqrt(50)


## Figure
## Read Results
pdf(file = "Figures/caseStudy_pih_rocComp.pdf",width = 2.5 ,height = 1.7)
ggplot(cc,aes(type,AUC,col = type))+
  geom_errorbar(aes(ymin = ymin,ymax = ymax),position = "dodge")+
  geom_point(position = position_dodge2(width = .9))+
  theme_bw()+
  ylab("AUC")+
  theme(axis.text = element_text(size = 7),
        plot.title = element_text(size = 7,hjust = .5),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold",size = 7),
        #axis.title.x = element_text(size = 8,face = "bold"),
        axis.title = element_text(size = 7,face = "bold"),
        axis.title.x = element_blank(),
        #axis.title.x =  = element_blank(),
        legend.text = element_text(size = 7),
        #legend.direction = "horizontal",
        legend.margin = margin(0,0,0,0,unit = "cm"),legend.box.spacing = unit(0,units = "in"), 
        legend.title = element_blank(),
        #legend.box.background = element_rect(colour = "black"),legend.key = element_rect(colour = "black"),
        legend.position = "top",
        legend.key.size = unit(.15,units = "in"),
        #legend.background = element_rect(colour = "black"),
        axis.text.x = element_blank(),
        plot.caption = element_text(size = 7))

dev.off()





# PCA Analysis ------------------------------------------------------------

## PCA
feature.df = data.frame(feature.df)
pc = prcomp(feature.df[,-1])
pc.df = data.frame(labels = feature.df[,1],pc$x)
plloading = data.frame(taxa = rownames(pc$rotation),pc$rotation)
plloading$taxa = str_replace(plloading$taxa,pattern = "___"," / ")

summary(pc)


pdf(file = "Figures/caseStudy_pih_pca.pdf",width = 2.5 ,height = 1.7)
ggplot(pc.df,aes(PC1,PC2,col = labels))+
  geom_point(size = 3,alpha = .75)+
  theme_bw()+
  xlab(paste0("PC1 (78.48% of Variation)"))+
  ylab(paste0("PC2 (3.88% of Variation)"))+
  theme(legend.position = "top",legend.title = element_blank())+
  ggsci::scale_color_lancet()+
  theme(axis.text = element_text(size = 7),panel.grid = element_blank(),
        axis.title.x = element_text(size = 7,face = "bold"),
        axis.title.y = element_text(size = 7,face = "bold"),
        legend.text = element_text(size = 7),
        #legend.direction = "horizontal",
        legend.margin = margin(0,0,.1,0,unit = "cm"),
        legend.box.spacing = unit(0,units = "in"), 
        legend.title = element_blank(),
        #legend.box.background = element_rect(colour = "black"),legend.key = element_rect(colour = "black"),
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        #legend.background = element_rect(colour = "black"),
        axis.text.x = element_text(size = 7),plot.caption = element_text(size = 7))
dev.off()






# Log ratio means ---------------------------------------------------------
rp = data.frame(id = 1:nrow(feature.df),Status = feature.df[,1],feature.df[,-1])
rp = rp %>%
  gather("taxa","log_ratio",3:ncol(rp))
rp$taxa = str_replace(rp$taxa,pattern = "___"," / ")
dd = rp %>%
  group_by(taxa) %>%
  wilcox_test(data =., log_ratio ~ Status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  dplyr::arrange(dplyr::desc(-p))


rp = data.frame(id = 1:nrow(feature.df),Status = feature.df[,1],feature.df[,-1])
rp = rp %>% 
  gather("taxa","log_ratio",3:ncol(rp)) %>% 
  group_by(Status,taxa) %>% 
  summarise(n = n(),
            lb =  mean(log_ratio) - ((1.96*mean(log_ratio))/sqrt(n)),
            ub =  mean(log_ratio) + ((1.96*mean(log_ratio))/sqrt(n)),
            log_ratio = mean(log_ratio),
  )
rp$taxa = str_replace(rp$taxa,pattern = "___"," / ")

rp = left_join(rp,dd)
rp = left_join(rp,plloading[,1:2])
rp = rp %>% 
  arrange(desc(PC1))
rp$taxa = factor(rp$taxa,levels = unique(rp$taxa))

pdf(file = "Figures/caseStudy_pih_lrMeans.pdf",width = 5.5 ,height = 4.5)
ggplot(rp,aes(taxa,log_ratio,fill  = Status,label = p.adj.signif))+
  geom_col(position = position_identity(),alpha = .5,col = "black",width = .5)+
  geom_errorbar(aes(ymin = lb,ymax = ub),width = .15,size = .5)+
  geom_point(aes(fill = Status),pch = 21,size = 2)+
  theme_bw()+
  coord_flip()+
  ggtitle("Log Ratio Means")+
  geom_text(nudge_y = 1,size = 3,fontface = "bold",y = -8,col = "red" )+
  scale_y_continuous(limits = c(-8.5,7))+
  ylab("Denominator enriched <-- Log Ratio --> Numerator Enriched")+
  #geom_hline(yintercept = 0)+
  ggsci::scale_fill_lancet()+
  ggsci::scale_color_lancet()+
  theme(legend.position = "top",plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 7),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 7),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        #legend.background = element_rect(colour = "black")
  )
dev.off()



## loadings
rp1=  rp %>% 
  group_by(taxa) %>% 
  summarise(loading = mean(PC1))
rp1 = left_join(rp1,vi)
rp1$taxa = factor(rp1$taxa,levels = levels(rp$taxa))


pdf(file = "Figures/caseStudy_pih_loadings.pdf",width = 2.5 ,height = 4.5)
ggplot(rp1,aes(taxa,loading,fill = Importance))+
  geom_col(position = position_identity(),alpha = 1,col = "black",width = .5)+
  coord_flip()+
  ggtitle("Log Ratio Means")+
  #geom_text(nudge_y = 1,size = 3,fontface = "bold",y = -8 )+
  #scale_y_continuous(limits = c(-8,7))+
  ylab("PC1 Loading")+
  scale_fill_distiller(palette = "Purples",direction = 1)+
  theme_bw()+
  theme(legend.position = c(.2,.5),
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.y =element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 7),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        #legend.background = element_rect(colour = "black")
  )
dev.off()



##-----------------------------------------*
## LR Network ####
##-----------------------------------------*

imp.df = data.frame(Ratio = vi$Ratio,Imp = vi$Importance)
keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g)$weight = if_else((imp.df$Imp)<0,0,(imp.df$Imp))
#g = minimum.spanning.tree(graph = g,weights = -E(g)$weight)
w_degree  = strength(g,mode = "total")
plot(g,layout = layout_with_fr)
toGephi(g,"Output/caseStudy_pih")


imp = glm.mdl1$finalModel
imp = randomForest::importance(imp)


















