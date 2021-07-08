rm(list = ls())
gc()

## Load/Install Required Packages
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
library(rstatix)



## Load Helper Functions
source("Helper Functions/functions.R")


## Start Cluster
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)




# Read Data and Compute Signatures ---------------------------------------------------------------

dat =  read_csv("Output/diabimmune_f6.csv")
md = read_csv("Output/diabimmune_f6_metadata.csv");md = data.frame(md)
sep.true = selectionEnergy.scaled(inputData = dat,nreps_energy = 1000,
                                  dcv_nfold = 1,patience = 500)

x = sep.true$retainedDCV
feature.df = sep.true$finalSubset
lrs_dcv = sep.true$retainedDCV
lrs_dcv = lrs_dcv %>%
  filter(Ratio%in% colnames(feature.df[,-1]))
xx = getLogratioFromList(Ratio = colnames(feature.df[,-1]),raMatrix = dat,Class = "test")
feature.df = data.frame(Status = dat[,1],xx)
write_csv(lrs_dcv,"Output/diabimmune_f6_signature_dcv.csv")
write_csv(feature.df,"Output/diabimmune_f6_signature.csv")


dat =  read_csv("Output/diabimmune_6-12.csv")
md = read_csv("Output/diabimmune_6-12_metadata.csv");md = data.frame(md)
sep.true = selectionEnergy.scaled(inputData = dat,nreps_energy = 1000,
                                  dcv_nfold = 1,patience = 500)
x = sep.true$retainedDCV
feature.df = sep.true$finalSubset
lrs_dcv = sep.true$retainedDCV
lrs_dcv = lrs_dcv %>%
  filter(Ratio%in% colnames(feature.df[,-1]))
xx = getLogratioFromList(Ratio = colnames(feature.df[,-1]),raMatrix = dat,Class = "test")
feature.df = data.frame(Status = dat[,1],xx)
write_csv(lrs_dcv,"Output/diabimmune_6-12_signature_dcv.csv")
write_csv(feature.df,"Output/diabimmune_6-12_signature.csv")




# Permutation Result ------------------------------------------------------

### read SelEnergyPErm Resutls
fnames = dir("Results/Case_Study/")
bool = str_detect(fnames,pattern = "diab")
fnames = fnames[bool]
results.df = data.frame()
for(f in fnames){
  ph = read_csv(file = paste0("Results/Case_Study/",f))
  results.df = rbind(results.df,ph)
}

time = unique(results.df$fname)
selEnergyPerm.results = data.frame()
scaledResults = data.frame()
for(i in time){
  t1.null = results.df %>%
    filter(fname==i,Type=="null")
  t1.empi = results.df %>%
    filter(fname==i,Type!="null")
  ph = rbind(t1.null,t1.empi[1,])
  ph$tstat  = scale(ph$tstat)

  scaledResults =rbind(scaledResults,ph )
  hist(t1.null$tstat)
  tstat = t1.empi$tstat[1]
  p = (sum(t1.null$tstat>tstat)+1)/ (nrow(t1.null)+1)
  selEnergyPerm.results = rbind(selEnergyPerm.results,data.frame(Time = i,tstat,p))
}

selEnergyPerm.results$p.adj = p.adjust(selEnergyPerm.results$p,method = "BH")
selEnergyPerm.results$logp = -1*log(selEnergyPerm.results$p.adj)
selEnergyPerm.results$Time  = str_split(selEnergyPerm.results$Time,pattern = "_",simplify = T)[,2]
selEnergyPerm.results$col = if_else(selEnergyPerm.results$p.adj<0.05,"darkgreen","grey")
ggplot(selEnergyPerm.results,aes(Time,logp))+
  geom_col(fill = selEnergyPerm.results$col,col ="black")+
  ylab("-log(p.adj)")+
  geom_hline(yintercept = -log(0.05),lty = "dashed",col = "red")+
  theme_classic()


scaledResults$fname  = str_split(scaledResults$fname,pattern = "_",simplify = T)[,2]
unique(scaledResults$fname)
scaledResults$fname = factor(scaledResults$fname,levels = c("f6","6-12","12-18","18-24"),
                             labels = c("(0,6]","(6-12]","(12-18]","(18-24]") )
scaledResults$label = "null"
scaled.null = scaledResults %>%
  filter(Type=="null")
scaled.emp = scaledResults %>%
  filter(Type!="null")
scaled.emp$label = paste0("p = ",round(selEnergyPerm.results$p.adj,digits = 3))
scaled.emp$col = if_else(selEnergyPerm.results$p.adj<0.05,"red","black")

tiff(filename = "Figures/selEneryResults_diabimmune.tiff",width = 3.5,height = 2.2,units = "in",res = 300)

pdf(file = "Figures/caseStudy_diabimmune_sepResults.pdf",width = 3 ,height = 2.5)
ggplot(scaled.null,aes(fname,tstat))+
  geom_violin(trim = T,width = .75)+
  geom_jitter(alpha = .1,width = .1,col = "grey40")+
  theme_bw()+
  xlab("Collection Month")+
  ylab("SelEnergyPerm Test Statistic")+
  scale_y_continuous(limits = c(-2.5,6))+
  geom_point(data = scaled.emp,aes(fname,tstat),col = "black",size = 2,pch =21,fill = scaled.emp$col)+
  #geom_line(data = scaled.emp,aes(fname,tstat,group =1),col = "red",size = .75)+
  geom_text(data =scaled.emp,aes(label = label),col = scaled.emp$col,nudge_y = 1,size = 3,fontface = "bold",y = 6 )+
  theme(axis.text = element_text(size = 7),
        panel.grid = element_blank(),
        #axis.title.x = element_text(size = 8,face = "bold"),
        axis.title = element_text(size = 7,face = "bold"),
        legend.text = element_text(size = 7),
        #legend.direction = "horizontal",
        legend.margin = margin(0,0,0,0,unit = "cm"),legend.box.spacing = unit(0,units = "in"),
        legend.title = element_blank(),
        #legend.box.background = element_rect(colour = "black"),legend.key = element_rect(colour = "black"),
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        #legend.background = element_rect(colour = "black"),
        axis.text.x = element_text(size = 7),
        plot.caption = element_text(size = 7))
dev.off()






# AUC Comparision ---------------------------------------------------------

## first 6
dat = read_csv("Output/diabimmune_f6.csv")
lrs = calcLogRatio(dat)
feature.df = read_csv("Output/diabimmune_f6_signature.csv");feature.df = data.frame(feature.df)
md = read_csv("Output/diabimmune_f6_metadata.csv");md =data.frame(md)

classes = unique(feature.df$Status)
## Train Control
tc1  <- trainControl(method="none",
                     repeats = 1,
                     number = 5,
                     seeds = NULL,
                     classProbs = TRUE,
                     savePredictions = T,
                     allowParallel = TRUE,
                     summaryFunction = caret::multiClassSummary
)

## Stratified K-Fold Leave-one-sample-out
## Which samples to keep
classes = unique(feature.df$Status)
nfolds = 10
md1 = md %>%
  group_by(subjectID,Status) %>%
  summarise(n = n())
folds_seeds = list()
for(s in 1:20){
  f = caret::createFolds(md1$Status,k = nfolds,list = F)
  md1$folds = f
  md2 = left_join(md,md1)
  folds_seeds[[s]]  = data.frame(seed = s,folds = md2$folds)
}

## All logratios
pred = data.frame()
preds = data.frame()
nreps = 20
for(sd in 1:nreps){
  flds = folds_seeds[[sd]]

  tp = foreach::foreach(f = 1:max(flds$folds),.combine = rbind,.packages = "caret")%dopar%{

    ## partition data
    outer_indices = f==flds$folds
    trainx = lrs[!outer_indices,-1]
    ytrain = lrs[!outer_indices,1]
    testx = lrs[outer_indices,-1]
    ytest = lrs[outer_indices,1]

    glm.mdl1 = caret::train(x = trainx,
                            y = ytrain,
                            probMethod = "softmax",
                            method = "pls",
                            tuneGrid = expand.grid(ncomp = 2),
                            metric = "AUC",
                            trControl = tc1
    )
    p = data.frame(Seed = sd,
                   Fold = f,
                   md[outer_indices,],
                   predict.train(glm.mdl1,newdata = testx,type = "prob")
    )
    rc = pROC::auc(ytest,p[,as.character(classes[1])])
    p$auc = as.numeric(rc)
    p
  }
  preds = rbind(preds,tp)

  message(sd)


}
auc.ovr = preds %>%
  group_by(Seed) %>%
  summarise(auc= mean(auc))
auc.ovr$type  = "all_logratios"
auc.ovr$time = "[0,6)"

## Signature
pred = data.frame()
preds = data.frame()
nreps = 20
for(sd in 1:nreps){
  flds = folds_seeds[[sd]]

  tp = foreach::foreach(f = 1:max(flds$folds),.combine = rbind,.packages = "caret")%dopar%{

    ## partition data
    outer_indices = f==flds$folds
    trainx = feature.df[!outer_indices,-1]
    ytrain = feature.df[!outer_indices,1]
    testx = feature.df[outer_indices,-1]
    ytest = feature.df[outer_indices,1]

    glm.mdl1 = caret::train(x = trainx,
                            y = ytrain,
                            probMethod = "softmax",
                            method = "pls",
                            tuneGrid = expand.grid(ncomp = 2),
                            metric = "AUC",
                            trControl = tc1
    )
    p = data.frame(Seed = sd,
                   Fold = f,
                   md[outer_indices,],
                   predict.train(glm.mdl1,newdata = testx,type = "prob")
    )
    rc = pROC::auc(ytest,p[,as.character(classes[1])])
    p$auc = as.numeric(rc)
    p
  }
  preds = rbind(preds,tp)

  message(sd)


}
auc.signature = preds %>%
  group_by(Seed) %>%
  summarise(auc= mean(auc))
auc.signature$type  = "signature"
auc.signature$time = "[0,6)"

f6 = rbind(auc.ovr,auc.signature)

## 6-12 months
dat = read_csv("Output/diabimmune_6-12.csv")
lrs = calcLogRatio(dat)
feature.df = read_csv("Output/diabimmune_6-12_signature.csv");feature.df = data.frame(feature.df)
md = read_csv("Output/diabimmune_6-12_metadata.csv");md =data.frame(md)

classes = unique(feature.df$Status)
## Train Control
tc1  <- trainControl(method="none",
                     repeats = 1,
                     number = 5,
                     seeds = NULL,
                     classProbs = TRUE,
                     savePredictions = T,
                     allowParallel = TRUE,
                     summaryFunction = caret::multiClassSummary
)

## Stratified K-Fold Leave-one-sample-out
## Which samples to keep
classes = unique(feature.df$Status)
nfolds = 10
md1 = md %>%
  group_by(subjectID,Status) %>%
  summarise(n = n())
folds_seeds = list()
for(s in 1:20){
  f = caret::createFolds(md1$Status,k = nfolds,list = F)
  md1$folds = f
  md2 = left_join(md,md1)
  folds_seeds[[s]]  = data.frame(seed = s,folds = md2$folds)
}

## All logratios
pred = data.frame()
preds = data.frame()
nreps = 20
for(sd in 1:nreps){
  flds = folds_seeds[[sd]]

  tp = foreach::foreach(f = 1:max(flds$folds),.combine = rbind,.packages = "caret")%dopar%{

    ## partition data
    outer_indices = f==flds$folds
    trainx = lrs[!outer_indices,-1]
    ytrain = lrs[!outer_indices,1]
    testx = lrs[outer_indices,-1]
    ytest = lrs[outer_indices,1]

    glm.mdl1 = caret::train(x = trainx,
                            y = ytrain,
                            probMethod = "softmax",
                            method = "pls",
                            tuneGrid = expand.grid(ncomp = 2),
                            metric = "AUC",
                            trControl = tc1
    )
    p = data.frame(Seed = sd,
                   Fold = f,
                   md[outer_indices,],
                   predict.train(glm.mdl1,newdata = testx,type = "prob")
    )
    rc = pROC::auc(ytest,p[,as.character(classes[1])])
    p$auc = as.numeric(rc)
    p
  }
  preds = rbind(preds,tp)

  message(sd)


}
auc.ovr = preds %>%
  group_by(Seed) %>%
  summarise(auc= mean(auc))
auc.ovr$type  = "all_logratios"
auc.ovr$time = "[6,12)"

## Signature
pred = data.frame()
preds = data.frame()
nreps = 20
for(sd in 1:nreps){
  flds = folds_seeds[[sd]]

  tp = foreach::foreach(f = 1:max(flds$folds),.combine = rbind,.packages = "caret")%dopar%{

    ## partition data
    outer_indices = f==flds$folds
    trainx = feature.df[!outer_indices,-1]
    ytrain = feature.df[!outer_indices,1]
    testx = feature.df[outer_indices,-1]
    ytest = feature.df[outer_indices,1]

    glm.mdl1 = caret::train(x = trainx,
                            y = ytrain,
                            probMethod = "softmax",
                            method = "pls",
                            tuneGrid = expand.grid(ncomp = 2),
                            metric = "AUC",
                            trControl = tc1
    )
    p = data.frame(Seed = sd,
                   Fold = f,
                   md[outer_indices,],
                   predict.train(glm.mdl1,newdata = testx,type = "prob")
    )
    rc = pROC::auc(ytest,p[,as.character(classes[1])])
    p$auc = as.numeric(rc)
    p
  }
  preds = rbind(preds,tp)

  message(sd)


}
auc.signature = preds %>%
  group_by(Seed) %>%
  summarise(auc= mean(auc))
auc.signature$type  = "signature"
auc.signature$time = "[6,12)"

s612 =  rbind(auc.ovr,auc.signature)



## COmbine and plot results
all.df = rbind(f6,s612) %>%
  group_by(type,time) %>%
  summarise(mean_cl_normal(auc) )

pdf(file = "Figures/caseStudy_diabimmune_rocComparision.pdf",width = 2.5 ,height = 2.5)
ggplot(all.df,aes(type,y,col = type))+
  geom_errorbar(aes(ymin = ymin,ymax = ymax),position = "dodge")+
  geom_point(position = position_dodge2(width = .9))+
  theme_bw()+
  facet_grid(.~time)+
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




# DCV Strength Analysis  --------------------------------------------------

## Month 0
dcv = read_csv(file = "Output/diabimmune_f6_signature_dcv.csv")
trainx = read.csv(file = "Output/diabimmune_f6_signature.csv")
glm.mdl1 = caret::train(x = trainx[,-1],
                        y = trainx[,1],
                        probMethod = "softmax",
                        method = "pls",
                        tuneGrid = expand.grid(ncomp = 2),
                        metric = "AUC",
                        trControl = tc1
)
vi = varImp(glm.mdl1)
vi = data.frame(Ratio = rownames(vi$importance),rowmean = as.numeric(vi$importance$Overall))
str_month0 = dcvStrength(dcv_mat = list(dcv = vi))
str_month0 = str_month0 %>%
  top_n(5,Str)
str_month0$month = "[0,6)"
## Month 1
dcv = read_csv(file = "Output/diabimmune_6-12_signature_dcv.csv")
trainx = read.csv(file = "Output/diabimmune_6-12_signature.csv")
glm.mdl1 = caret::train(x = trainx[,-1],
                        y = trainx[,1],
                        probMethod = "softmax",
                        method = "pls",
                        tuneGrid = expand.grid(ncomp = 2),
                        metric = "AUC",
                        trControl = tc1
)
vi = varImp(glm.mdl1)
vi = data.frame(Ratio = rownames(vi$importance),rowmean = as.numeric(vi$importance$Overall))
str_month1 = dcvStrength(dcv_mat = list(dcv = vi))
str_month1 = str_month1 %>%
  top_n(5,Str)
str_month1$month = "[6,12)"
str_all = rbind(str_month0,str_month1)
str_all = str_all %>%
  group_by(Node,month) %>%
  summarise(Strength = sum(Str)) %>%
  spread("Node","Strength",fill = 0)
str_all1 = data.frame(month = str_all$month,clo(str_all[,-1]))
str_all = str_all1 %>%
  gather("Node","Strength",2:ncol(str_all))
ord = str_all %>%
  group_by(Node) %>%
  summarise(str = mean(Strength)) %>%
  arrange(desc(str))
str_all$Node = factor(str_all$Node,levels = ord$Node)

pdf(file = "Figures/caseStudy_diabimmune_dcvImp.pdf",width = 3 ,height = 2.5)
ggplot(str_all,aes(month,Strength,fill = Node))+
  geom_col(col = "black")+
  #coord_flip()+
  theme_minimal()+
  ylab("Relative PLS-DA Var. Importance Strength")+
  ggsci::scale_fill_simpsons()+
  theme(legend.position = "right",
        #plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 7,face = "bold"),
        axis.title.x = element_blank(),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-1,-10,-10),
        axis.text = element_text(size = 7),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")
  )
dev.off()






# LR Means Network --------------------------------------------------------

## Month 1
ra.dat = read_csv("Output/diabimmune_f6.csv")
feature.df = read_csv("Output/diabimmune_f6_signature.csv");feature.df = data.frame(feature.df)
md = read_csv("Output/diabimmune_f6_metadata.csv");md = data.frame(md)
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

classes = unique(rp$Status)
## CLass 1
rpc1 = rp %>%
  filter(Status==classes[1]) %>%
  arrange(desc(log_ratio))
rpc1$edgeweight = abs(rpc1$log_ratio)
rpc1$edgeweightCol = if_else(rpc1$log_ratio>0,"pos","neg")
imp.df = data.frame(Ratio = rpc1$taxa,Imp = rpc1$edgeweight)
keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = " / ",remove = F)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g)$weight = imp.df$Imp
#g = minimum.spanning.tree(graph = g,weights = -E(g)$weight)
w_degree  = strength(g,mode = "total")
w_degree = data.frame(Taxa = names(w_degree),Str = as.numeric(w_degree))

ra = ra.dat %>%
  dplyr::select(Status,one_of(w_degree$Taxa)) %>%
  filter(Status==classes[1])
ra = data.frame(clo(ra[,-1]))
cmeans =compositions::mean.acomp(acomp(ra))
ra = data.frame(ID = names(cmeans),meanAbund = (as.numeric(cmeans)))

w_degree  = strength(g,mode = "total")
w_degree = data.frame(ID = names(w_degree),
                      Label = names(w_degree),
                      Str = as.numeric(w_degree))
w_degree1 = w_degree %>%
  arrange(desc(Str)) %>%
  top_n(n = 5,wt = Str)
w_degree = w_degree %>%
 mutate(Label = if_else(ID %in% w_degree1$Label,Label,""))

plot(g,layout = layout_with_fr)
toGephi(g,"Output/diabimmune_f6_FA")
## Update Attributes
attbs = read_csv("Output/diabimmune_f6_FA_attributes_.csv")
attbs = left_join(attbs[,-2],w_degree,by = "ID")
attbs = left_join(attbs,ra)
write_csv(attbs,"Output/diabimmune_f6_FA_attributes_.csv")
ew = read_csv("Output/diabimmune_f6_FA_edgeWeights_.csv")
ew$taxa = paste0(ew$Source," / ",ew$Target)
ew = left_join(ew,rpc1,by = "taxa")
write_csv(ew,"Output/diabimmune_f6_FA_edgeWeights_.csv")

library(rstatix)
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
rp$taxa = factor(rp$taxa,levels = rpc1$taxa)
rp = left_join(rp,dd)
tiff(filename = "Figures/Manuscript Figures/Additional File X3.tiff",width = 6,height = 8,units = "in",res = 300)
ggplot(rp,aes(reorder(taxa,log_ratio),log_ratio,fill  =Status,label = p.adj.signif ))+
  geom_col(position = position_identity(),alpha = .5,col = "black",width = .5)+
  geom_errorbar(aes(ymin = lb,ymax = ub),width = .15,size = .5)+
  #geom_point(aes(fill = Status),pch = 21,size = 2)+
  theme_bw()+
  coord_flip()+
  scale_y_continuous(limits = c(-15,10))+
  geom_text(nudge_y = 1,size = 3,fontface = "bold",y = -15 )+
  ggtitle("(0,6] Months")+
  #scale_x_discrete(labels = parse(text = loading_$ratio_name))+
  ylab("Denominator enriched <-- Log Ratio --> Numerator Enriched")+
  #geom_hline(yintercept = 0)+
  ggsci::scale_fill_jco()+
  ggsci::scale_color_jco()+
  theme(legend.position = "top",
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 7),
        axis.title.y = element_blank(),
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


## CLass 2
classes[2]
rpc1 = rp %>%
  filter(Status==classes[2]) %>%
  arrange(desc(log_ratio))
rpc1$edgeweight = abs(rpc1$log_ratio)
rpc1$edgeweightCol = if_else(rpc1$log_ratio>0,"pos","neg")
imp.df = data.frame(Ratio = rpc1$taxa,Imp = rpc1$edgeweight)
keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = " / ",remove = F)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g)$weight = imp.df$Imp
#g = minimum.spanning.tree(graph = g,weights = -E(g)$weight)
w_degree  = strength(g,mode = "total")
w_degree = data.frame(Taxa = names(w_degree),Str = as.numeric(w_degree))

ra = ra.dat %>%
  dplyr::select(Status,one_of(w_degree$Taxa)) %>%
  filter(Status==classes[1])
ra = data.frame(clo(ra[,-1]))
cmeans =compositions::mean.acomp(acomp(ra))
ra = data.frame(ID = names(cmeans),meanAbund = (as.numeric(cmeans)))

w_degree  = strength(g,mode = "total")
w_degree = data.frame(ID = names(w_degree),
                      Label = names(w_degree),
                      Str = as.numeric(w_degree))
w_degree1 = w_degree %>%
  arrange(desc(Str)) %>%
  top_n(n = 5,wt = Str)
w_degree = w_degree %>%
  mutate(Label = if_else(ID %in% w_degree1$Label,Label,""))

plot(g,layout = layout_with_fr)
toGephi(g,"Output/diabimmune_f6_none")
## Update Attributes
attbs = read_csv("Output/diabimmune_f6_none_attributes_.csv")
attbs = left_join(attbs[,-2],w_degree,by = "ID")
attbs = left_join(attbs,ra)
write_csv(attbs,"Output/diabimmune_f6_none_attributes_.csv")
ew = read_csv("Output/diabimmune_f6_none_edgeWeights_.csv")
ew$taxa = paste0(ew$Source," / ",ew$Target)
ew = left_join(ew,rpc1,by = "taxa")
write_csv(ew,"Output/diabimmune_f6_none_edgeWeights_.csv")


rpc2 = rp %>%
  filter(Status==classes[2]) %>%
  arrange(desc(log_ratio))
rp$taxa = factor(rp$taxa,levels = rpc1$taxa)

#tiff(filename = "Figures/lrmeans_foodAllergy.tiff",width = 3,height = 4,units = "in",res = 300)
ggplot(rp,aes(taxa,log_ratio,fill  =Status))+
  geom_col(position = position_identity(),alpha = .5,col = "black",width = .5)+
  geom_errorbar(aes(ymin = lb,ymax = ub),width = .15,size = .5)+
  geom_point(aes(fill = Status),pch = 21,size = 2)+
  theme_bw()+
  coord_flip()+
  #ggtitle("Log Ratio Means")+
  #scale_x_discrete(labels = parse(text = loading_$ratio_name))+
  ylab("Denominator enriched <-- Log Ratio --> Numerator Enriched")+
  #geom_hline(yintercept = 0)+
  ggsci::scale_fill_jco()+
  ggsci::scale_color_jco()+
  theme(legend.position = "top",plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 7),
        axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 7),
        #panel.grid = element_blank(),
        #legend.key.size = unit(.05,units = "in"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        #legend.background = element_rect(colour = "black")
  )
dev.off()



## Month 2
ra.dat = read_csv("Output/diabimmune_6-12.csv")
feature.df = read_csv("Output/diabimmune_6-12_signature.csv");feature.df = data.frame(feature.df)
md = read_csv("Output/diabimmune_6-12_metadata.csv");md = data.frame(md)
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

classes = unique(rp$Status)
## CLass 1
classes[1]
rpc1 = rp %>%
  filter(Status==classes[1]) %>%
  arrange(desc(log_ratio))
rpc1$edgeweight = abs(rpc1$log_ratio)
rpc1$edgeweightCol = if_else(rpc1$log_ratio>0,"pos","neg")
imp.df = data.frame(Ratio = rpc1$taxa,Imp = rpc1$edgeweight)
keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = " / ",remove = F)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g)$weight = imp.df$Imp
#g = minimum.spanning.tree(graph = g,weights = -E(g)$weight)
w_degree  = strength(g,mode = "total")
w_degree = data.frame(Taxa = names(w_degree),Str = as.numeric(w_degree))

ra = ra.dat %>%
  dplyr::select(Status,one_of(w_degree$Taxa)) %>%
  filter(Status==classes[1])
ra = data.frame(clo(ra[,-1]))
cmeans =compositions::mean.acomp(acomp(ra))
ra = data.frame(ID = names(cmeans),meanAbund = (as.numeric(cmeans)))

w_degree  = strength(g,mode = "total")
w_degree = data.frame(ID = names(w_degree),
                      Label = names(w_degree),
                      Str = as.numeric(w_degree))
w_degree1 = w_degree %>%
  arrange(desc(Str)) %>%
  top_n(n = 5,wt = Str)
w_degree = w_degree %>%
  mutate(Label = if_else(ID %in% w_degree1$Label,Label,""))

plot(g,layout = layout_with_fr)
toGephi(g,"Output/diabimmune_612_FA")
## Update Attributes
attbs = read_csv("Output/diabimmune_612_FA_attributes_.csv")
attbs = left_join(attbs[,-2],w_degree,by = "ID")
attbs = left_join(attbs,ra)
write_csv(attbs,"Output/diabimmune_612_FA_attributes_.csv")
ew = read_csv("Output/diabimmune_612_FA_edgeWeights_.csv")
ew$taxa = paste0(ew$Source," / ",ew$Target)
ew = left_join(ew,rpc1,by = "taxa")
write_csv(ew,"Output/diabimmune_612_FA_edgeWeights_.csv")




## CLass 2
classes[2]
rpc1 = rp %>%
  filter(Status==classes[2]) %>%
  arrange(desc(log_ratio))
rpc1$edgeweight = abs(rpc1$log_ratio)
rpc1$edgeweightCol = if_else(rpc1$log_ratio>0,"pos","neg")
imp.df = data.frame(Ratio = rpc1$taxa,Imp = rpc1$edgeweight)
keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = " / ",remove = F)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g)$weight = imp.df$Imp
#g = minimum.spanning.tree(graph = g,weights = -E(g)$weight)
w_degree  = strength(g,mode = "total")
w_degree = data.frame(Taxa = names(w_degree),Str = as.numeric(w_degree))

ra = ra.dat %>%
  dplyr::select(Status,one_of(w_degree$Taxa)) %>%
  filter(Status==classes[1])
ra = data.frame(clo(ra[,-1]))
cmeans =compositions::mean.acomp(acomp(ra))
ra = data.frame(ID = names(cmeans),meanAbund = (as.numeric(cmeans)))

w_degree  = strength(g,mode = "total")
w_degree = data.frame(ID = names(w_degree),
                      Label = names(w_degree),
                      Str = as.numeric(w_degree))
w_degree1 = w_degree %>%
  arrange(desc(Str)) %>%
  top_n(n = 5,wt = Str)
w_degree = w_degree %>%
  mutate(Label = if_else(ID %in% w_degree1$Label,Label,""))

plot(g,layout = layout_with_fr)
toGephi(g,"Output/diabimmune_612_none")
## Update Attributes
attbs = read_csv("Output/diabimmune_612_none_attributes_.csv")
attbs = left_join(attbs[,-2],w_degree,by = "ID")
attbs = left_join(attbs,ra)
write_csv(attbs,"Output/diabimmune_612_none_attributes_.csv")
ew = read_csv("Output/diabimmune_612_none_edgeWeights_.csv")
ew$taxa = paste0(ew$Source," / ",ew$Target)
ew = left_join(ew,rpc1,by = "taxa")
write_csv(ew,"Output/diabimmune_612_none_edgeWeights_.csv")



library(rstatix)
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
rp$taxa = factor(rp$taxa,levels = rpc1$taxa)
rp = left_join(rp,dd)
tiff(filename = "Figures/Manuscript Figures/Additional File X2.tiff",width = 6,height = 8,units = "in",res = 300)
ggplot(rp,aes(reorder(taxa,log_ratio),log_ratio,fill  =Status,label = p.adj.signif ))+
  geom_col(position = position_identity(),alpha = .5,col = "black",width = .5)+
  geom_errorbar(aes(ymin = lb,ymax = ub),width = .15,size = .5)+
  #geom_point(aes(fill = Status),pch = 21,size = 2)+
  theme_bw()+
  coord_flip()+
  scale_y_continuous(limits = c(-6,7))+
  geom_text(nudge_y = 1,size = 3,fontface = "bold",y = -5.5 )+
  ggtitle("(6,12] Months")+
  #scale_x_discrete(labels = parse(text = loading_$ratio_name))+
  ylab("Denominator enriched <-- Log Ratio --> Numerator Enriched")+
  #geom_hline(yintercept = 0)+
  ggsci::scale_fill_jco()+
  ggsci::scale_color_jco()+
  theme(legend.position = "top",
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 7),
        axis.title.y = element_blank(),
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







# Traditional Methods Analysis-----------------------------------------------------

pno = c()
p = c()



dat =  read_csv("Output/diabimmune_f6.csv");dat = data.frame(dat)
md =read_csv("Output/diabimmune_f6_metadata.csv");md = data.frame(md)
blocks = md$permutation_blocks

dat =  read_csv("Output/diabimmune_6-12.csv");dat = data.frame(dat)
md =read_csv("Output/diabimmune_6-12_metadata.csv");md = data.frame(md)
blocks = md$permutation_blocks


dat =  read_csv("Output/diabimmune_18-24.csv");dat = data.frame(dat)
md =read_csv("Output/diabimmune_18-24_metadata.csv");md = data.frame(md)
blocks = md$permutation_blocks


dat =  read_csv("Output/diabimmune_18-24.csv");dat = data.frame(dat)
md =read_csv("Output/diabimmune_18-24_metadata.csv");md = data.frame(md)
blocks = md$permutation_blocks





dat.f = data.frame(Status = factor(dat[,1]), fastImputeZeroes( compositions::clo(dat[,-1])  ) )
lrs = calcLogRatio(dat.f)
message("logratios computed")
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
  #d.compdat = parallelDist::parDist(as.matrix(dat[,-1]),method = "bray")

})
message("distance matrix computed")


## anoim
library(permute)
h1 <- how(nperm = 1000,blocks = factor(blocks))
ano = vegan::anosim(x = d.compdat,grouping = lrs[,1],permutations = h1)
pno = c(pno,ano$signif)

## PERMANOVA
library(vegan)
ttt = vegan::adonis2(d.compdat ~ Type,
                     data =  data.frame(Type = lrs[,1]),
                     permutations = h1)
p = c(p,ttt$`Pr(>F)`[1])
p.adjust(p,method = "BH")


# ## ENrgy missing block design
# energy_aitch =  energy::eqdist.etest(x = d.compdat,sizes = sz$Freq,distance = T,R = 1000)
#
# ## PERMDIS
# mod = vegan::betadisper(d.compdat,group = lrs[,1])
# md1 = (vegan::permutest(mod,permutations = h1))

