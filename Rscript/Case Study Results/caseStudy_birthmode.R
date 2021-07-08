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

## Start Cluster
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)




## Overview Metadata
md.df = data.frame()
otu.df = data.frame()
for(month_ in c("month0","month1","month2","month3")){

  file_ = paste0("Output/deliveryMode_",month_,".csv")
  dat = readr::read_csv(file_);dat = data.frame(ID = 1:nrow(dat),dat)
  dat = tidyr::gather(dat,"Taxa","Count",3:ncol(dat))

  dat$month = month_
  otu.df = rbind(otu.df,dat)
  path_ = paste0("Output/deliveryMode_",month_,"_metadata.csv")
  md = readr::read_csv(path_);md =data.frame(md)
  md.df = rbind(md.df,md)
}

md.df %>%
  group_by(month,host_subject_id) %>%
  summarise(n = n()) %>%
  group_by(month) %>%
  summarise(n = n())



# Permutation Results -----------------------------------------------------

fnames = dir("Results/Case_Study/")
bool  = str_detect(fnames,"birth")
fnames = fnames[bool]
results.df = data.frame()
for(f in fnames){
  ph = read_csv(file = paste0("Results/Case_Study/",f))
  results.df = rbind(results.df,ph)
}
birth.results = data.frame(results.df)
time = unique(birth.results$fname)
selEnergyPerm.results = data.frame()
scaledResults = data.frame()
for(i in time[1:4]){
  t1.null = birth.results %>%
    filter(fname==i,Type=="null")
  t1.empi = birth.results %>%
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
scaledResults$fname = factor(scaledResults$fname,levels = unique(scaledResults$fname),
                             labels = c("month0 \n (h=40) \n (ces.=19 ,vag.=55)",
                                        "month1\n (h=42) \n (ces.=31,vag.=31)" ,
                                        "month2\n (h=38) \n (ces.=23,vag.=28)" ,
                                        "month3\n (h=37) \n (ces.=24,vag.=19)"))
scaledResults$label = "null"
scaled.null = scaledResults %>%
  filter(Type=="null")
scaled.emp = scaledResults %>%
  filter(Type!="null")
scaled.emp$label = paste0("p = ",round(selEnergyPerm.results$p.adj,digits = 3))
scaled.emp$col = if_else(selEnergyPerm.results$p.adj<0.05,"red","black")

pdf(file = "Figures/caseStudy_birth_sepResults.pdf",width = 3.25 ,height = 2.7)
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









# Compute and Store SelEnergyPerm Features --------------------------------

for(i in 0:3){
  f = paste0("Output/deliveryMode_month",i,".csv")
  dat = read_csv(file = f)
  fname = paste0("birthDelivery_month",i)
  set.seed(08272008)
  px = 500
  sep.true = selEnergyPermR::selectionEnergy.scaled(inputData = dat,
                                                    nreps_energy = 1000,
                                                    dcv_nfold = 1,
                                                    optimizationMetric = NULL,
                                                    patience = px)
  x = sep.true$retainedDCV
  ## Get final features ####
  feature.df = sep.true$finalSubset
  ## Retiecve DCV ####
  lrs_dcv = sep.true$retainedDCV
  rxDCV = sep.true$rawDCV[,2:6]
  lrs_dcv = lrs_dcv %>%
    filter(Ratio %in% colnames(feature.df[,-1]))
  write_csv(feature.df,file = paste0("Output/",fname,"_features.csv"))
  write_csv(lrs_dcv,file = paste0("Output/",fname,"_dcv.csv"))
}










# Visualize log ratio means -----------------------------------------------

i = 0
for(i in 0:2){
  ## Month 1
  f = paste0("Output/birthDelivery_month",i,"_features.csv")
  f1 = paste0("Output/deliveryMode_month",i,".csv")

  ra.dat = read_csv(f1)
  feature.df = read_csv(file = f)
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
  for(c in 1:2){
    rpc1 = rp %>%
      filter(Status==classes[c]) %>%
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
      filter(Status==classes[c])
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
    fname = paste0("birth_",i,"_",classes[c])

    ## Export data in gephi format
    toGephi(g,paste0("Output/",fname))
    ## Node Attributes
    fn = paste0("Output/",fname,"_attributes_.csv")
    attbs = read_csv(fn)
    attbs = left_join(attbs[,-2],w_degree,by = "ID")
    attbs = left_join(attbs,ra)
    write_csv(attbs,fn)
    ## Edge Attributes
    fn = paste0("Output/",fname,"_edgeWeights_.csv")
    ew = read_csv(fn)
    ew$taxa = paste0(ew$Source," / ",ew$Target)
    ew = left_join(ew,rpc1,by = "taxa")
    write_csv(ew,fn)

  }
}




# Log Ratio Mean Group Comparisons ----------------------------------------


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
rp = rp %>%
  arrange(desc(log_ratio))
rp$taxa = factor(rp$taxa,levels = unique(rp$taxa))
rp$col = if_else(rp$p.adj.signif=="ns","black","red")

ggplot(rp,aes(taxa,log_ratio,fill  =Status,label = p.adj.signif ))+
  geom_col(position = position_identity(),alpha = .5,col = "black",width = .5)+
  geom_errorbar(aes(ymin = lb,ymax = ub),width = .15,size = .5)+
  #geom_point(aes(fill = Status),pch = 21,size = 2)+
  theme_bw()+
  coord_flip()+
  scale_y_continuous(limits = c(-2.5,2.5))+
  geom_text(nudge_y = 1,size = 3,fontface = "bold",y = -2.5 ,col = rp$col )+
  ggtitle("(0,6] Months")+
  ylab("Denominator enriched <-- Log Ratio --> Numerator Enriched")+
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




# PLS Importance Strength Analysis ---------------------------------------------------

## Month 0
dcv = read_csv(file = "Output/birthDelivery_month0_dcv.csv")
trainx = read.csv(file = "Output/birthDelivery_month0_features.csv")
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
str_month0$month = "month0"
## Month 1
dcv = read_csv(file = "Output/birthDelivery_month1_dcv.csv")
trainx = read.csv(file = "Output/birthDelivery_month1_features.csv")
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
str_month1 =dcvStrength(dcv_mat = list(dcv = vi))
str_month1 = str_month1 %>%
  top_n(5,Str)
str_month1$month = "month1"
## Month 2
dcv = read_csv(file = "Output/birthDelivery_month2_dcv.csv")
trainx = read.csv(file = "Output/birthDelivery_month2_features.csv")
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
str_month2 = dcvStrength(dcv_mat = list(dcv = vi))
str_month2 = str_month2 %>%
  top_n(5,Str)
str_month2$month = "month2"
## Month 3
dcv = read_csv(file = "Output/birthDelivery_month3_dcv.csv")
trainx = read.csv(file = "Output/birthDelivery_month3_features.csv")
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
str_month3 = dcvStrength(dcv_mat = list(dcv = vi))
str_month3 = str_month3 %>%
  top_n(5,Str)
str_month3$month = "month3"
##combined
str_all = rbind(str_month0,str_month1) %>%
  rbind(str_month2) %>%
  rbind(str_month3) %>%
  separate(col = 1,into =c("Family","Genus"),remove = F,sep = "\\." )
str_all$Family[str_all$Family=="S24"] = "S24-7"
str_all = str_all %>%
  group_by(Family,month) %>%
  summarise(Strength = sum(Str)) %>%
  spread("Family","Strength",fill = 0)
str_all1 = data.frame(month = str_all$month,clo(str_all[,-1]))
str_all = str_all1 %>%
  gather("Family","Strength",2:ncol(str_all))
ord = str_all %>%
  group_by(Family) %>%
  summarise(str = mean(Strength)) %>%
  arrange(desc(str))

str_all$Family = factor(str_all$Family,levels = ord$Family)

pdf(file = "Figures/caseStudy_birth_viStrengthAnlaysis.pdf",width = 3.25 ,height = 2.7)
ggplot(str_all,aes(month,Strength,fill = Family))+
  geom_col(col = "black")+
  #coord_flip()+
  theme_minimal()+
  ylab("Relative PLS-DA Var. Importance Strength")+
  ggsci::scale_fill_simpsons()+
  theme(legend.position = "right",
        #plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 7),
        axis.title.x = element_blank(),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-5,-5,-5),
        axis.text = element_text(size = 7),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        #legend.background = element_rect(colour = "black")
  )
dev.off()








# AUC Comparison Stratified by delivery mode and host----------------------------------------------------------

## Train Parms
tc1  <- trainControl(method="cv",
                     repeats = 1,
                     number = 5,
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
nfolds = 5
nreps = 20
num_seeds = 20
combined.auc = data.frame()


## Compute AUC with Cross validation
md.df = data.frame()
otu.df = data.frame()
for(month_ in c("month0","month1","month2","month3")){

  file_ = paste0("Output/deliveryMode_",month_,".csv")
  dat = readr::read_csv(file_);dat = data.frame(ID = 1:nrow(dat),dat)
  dat = tidyr::gather(dat,"Taxa","Count",3:ncol(dat))

  dat$month = month_
  otu.df = rbind(otu.df,dat)
  path_ = paste0("Output/deliveryMode_",month_,"_metadata.csv")
  md = readr::read_csv(path_);md =data.frame(md)
  md.df = rbind(md.df,md)
}

otu.df = tidyr::spread(otu.df,"Taxa","Count")
md = md.df %>%
  dplyr::group_by(host_subject_id,delivery) %>%
  dplyr::summarise(n= n()) %>%
  group_by(delivery) %>%
  summarise(n = n())

for(month_ in c("month0","month1","month2")){
  ## first 6
  file_ = paste0("Output/deliveryMode_",month_,".csv")
  dat = read_csv(file_)
  lrs = calcLogRatio(dat)
  path_ = paste0("Output/birthDelivery_",month_,"_features.csv")
  feature.df = read_csv(path_);feature.df = data.frame(feature.df)
  xx = DiCoVarML::getLogratioFromList(Ratio = colnames(feature.df[,-1]),raMatrix = dat,Class = "test")
  feature.df = data.frame(Status = dat[,1],xx)
  path_ = paste0("Output/deliveryMode_",month_,"_metadata.csv")
  md = read_csv(path_);md =data.frame(md)
  table(dat[,1])
  pppp = md %>%
    group_by(host_subject_id,delivery) %>%
    summarise(n = n())
  classes = unique(feature.df$Status)
  ## Train Control
  tc1  <- trainControl(method="none",
                       repeats = 1,
                       number = 5,
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

  ## Stratified K-Fold Leave-one-sample-out
  ## Which samples to keep
  classes = unique(feature.df$Status)
  md1 = md %>%
    group_by(host_subject_id,delivery) %>%
    summarise(n = n())
  folds_seeds = list()
  for(s in 1:num_seeds){
    f = caret::createFolds(md1$delivery,k = nfolds,list = F)
    md1$folds = f
    md2 = left_join(md,md1)
    folds_seeds[[s]]  = data.frame(seed = s,folds = md2$folds)
  }

  ## All logratios
  pred = data.frame()
  preds = data.frame()
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
  auc.ovr$time = month_

  ## Signature
  pred = data.frame()
  preds = data.frame()
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
  auc.signature$time = month_
  f6 = rbind(auc.ovr,auc.signature)
  combined.auc = rbind(combined.auc,f6)
}
write_csv(combined.auc,"Results/birthStudy_aucComparision.csv")

## Read Results
combined.auc = read_csv("Results/birthStudy_aucComparision.csv")
combined.auc = combined.auc %>%
  group_by(type,time) %>%
  summarise(mean_cl_normal(auc) )

pdf(file = "Figures/caseStudy_birth_rocCom.pdf",width = 1.8 ,height = 2.23)
ggplot(combined.auc,aes(type,y,col = type))+
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












# Traditional Methods -----------------------------------------------------

month_ = "month3"
file_ = paste0("Data/deliveryMode_",month_,".csv")
dat = readr::read_csv(file_);dat = data.frame(dat)
path_ = paste0("Data/deliveryMode_",month_,"_metadata.csv")
md = readr::read_csv(path_);md =data.frame(md)
blocks = factor(md$host_subject_id)


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

## PERMANOVA
library(vegan)
ttt = vegan::adonis2(d.compdat ~ Type,
                     data =  data.frame(Type = lrs[,1]),
                     permutations = h1)
p = c(p,ttt$`Pr(>F)`[1])
p.adjust(p,method = "BH")



## ENrgy
energy_aitch =  energy::eqdist.etest(x = d.compdat,sizes = sz$Freq,distance = T,R = 1000)

## PERMDIS
mod = vegan::betadisper(d.compdat,group = lrs[,1])
md1 = (vegan::permutest(mod,permutations = h1))
