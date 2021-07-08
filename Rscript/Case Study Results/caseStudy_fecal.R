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


# Read Data ---------------------------------------------------------------
dat =  read_csv("Output/fecalcalprotecin_ALL.csv");dat = data.frame(dat)
md = read_csv("Output/fecalcalprotecin_ALL_metadata.csv");md = data.frame(md)
blocks = md$perm_blocks
sep.true = selectionEnergy.scaled(inputData = dat,nreps_energy = 1000,
                                  dcv_nfold = 1,patience = 500,
                                  optimizationMetric = "combinedF")
x = sep.true$retainedDCV
feature.df = sep.true$finalSubset
lrs_dcv = sep.true$retainedDCV
lrs = calcLogRatio(dat)
xx = DiCoVarML::getLogratioFromList(Ratio = colnames(feature.df[,-1]),raMatrix = dat,Class = "test")
feature.df = data.frame(Status = dat[,1],xx)



# Permutation Result ------------------------------------------------------

### read SelEnergyPErm Resutls
fnames = dir("Results/Case_Study/")
bool = str_detect(fnames,pattern = "ALL")
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

pdf(file = "Figures/caseStudy_fecalcalprotecin_sepResults.pdf",width = 3 ,height = 3)
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







# AUC Comparison ---------------------------------------------------------

## Train Control
tc1  <- trainControl(method="repeatedcv",
                     repeats = 20,
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

## All Ratios
classes = unique(feature.df$Status)
glm.mdl1=train(x = lrs[,-1],
               y = lrs[,1],
               probMethod = "softmax",
               method = "pls",
               tuneGrid = expand.grid(ncomp = 2),
               metric = "AUC",
               trControl = tc1
)
preds = glm.mdl1$pred
mroc = pROC::multiclass.roc(as.character(preds$obs),preds[,as.character(classes)])
all.roc = DiCoVarML::rocPlot(mroc)
all.roc$Comp = str_replace_all(string = all.roc$Comp,pattern = "abnormal/normal",replacement = "All Logratios")
all = glm.mdl1$results
all$type = "all logratios"

## Signature
glm.mdl.signature=train(x = feature.df[,-1],
                        y = feature.df[,1],
                        probMethod = "softmax",
                        method = "pls",
                        tuneGrid = expand.grid(ncomp = 2),
                        metric = "AUC",
                        trControl = tc1
)
preds = glm.mdl.signature$pred
mroc = pROC::multiclass.roc(as.character(preds$obs),preds[,as.character(classes)])
signature.roc = DiCoVarML::rocPlot(mroc)
signature.roc$Comp = str_replace_all(string = signature.roc$Comp,pattern = "abnormal/normal",replacement = "Signature")
signature.roc = glm.mdl.signature$results
signature.roc$type = "signature"


cc = rbind(all,signature.roc)
cc$ymin = cc$AUC - (1.96*cc$AUCSD)/sqrt(20)
cc$ymax = cc$AUC + (1.96*cc$AUCSD)/sqrt(20)

## Read Results
pdf(file = "Figures/caseStudy_fecal_rocComp.pdf",width = 2.5 ,height = 1.5)
ggplot(cc,aes(type,AUC,col = type))+
  geom_errorbar(aes(ymin = ymin,ymax = ymax),position = "dodge",width = .5)+
  geom_point(position = position_dodge2(width = .9))+
  theme_bw()+
  ylab("AUC")+
  scale_y_continuous(limits = c(.75,.9))+
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








# PLS Model ---------------------------------------------------------------
vi = varImp(glm.mdl.signature,useModel = T,scale = F)
vi  =data.frame(Ratio = rownames(vi$importance), Importance = vi$importance[,1])
vi = vi %>%
  arrange(desc(-Importance))
vi$Ratio = str_replace(vi$Ratio,pattern = "___"," / ")
vi$Ratio = factor(vi$Ratio,levels = vi$Ratio)
ggplot(vi,aes(Ratio,Importance))+
  geom_col(width = .65,col = "black")+
  theme_bw()+
  coord_flip()+
  ggtitle("PLS-DA Variable Importance")+
  ylab("")+
  scale_fill_distiller(palette = "Purples",direction = 1)+
  theme(legend.position = c(.2,.9),plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(face = "bold",size = 7),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size  = 7,hjust = 1),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 7),
        #panel.grid = element_blank(),
        legend.key.size = unit(.075,units = "in"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        #legend.background = element_rect(colour = "black")
  )








# Log ratio means ---------------------------------------------------------

library(rstatix)
rp = data.frame(id = 1:nrow(feature.df),Status = feature.df[,1],dx = md$diagnosis,feature.df[,-1])
rp = rp %>%
  gather("taxa","log_ratio",4:ncol(rp))
rp$taxa = str_replace(rp$taxa,pattern = "___"," / ")
dd = rp %>%
  group_by(taxa,dx) %>%
  wilcox_test(data =., log_ratio ~ Status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  dplyr::arrange(dplyr::desc(-p))

## specific signtature
dd1 = dd %>%
  mutate(signf = if_else(p.adj.signif=="ns",0,1)) %>%
  group_by(taxa) %>%
  summarise(n = sum(signf)) %>%
  arrange(desc(n))
dd$taxa = factor(dd$taxa,levels = dd1$taxa)

## Compute Means
rp = data.frame(id = 1:nrow(feature.df),Status = feature.df[,1],dx = md$diagnosis,feature.df[,-1])
rp = rp %>%
  gather("taxa","log_ratio",4:ncol(rp)) %>%
  group_by(Status,taxa,dx) %>%
  summarise(n = n(),
            lb =  mean(log_ratio) - ((1.96*mean(log_ratio))/sqrt(n)),
            ub =  mean(log_ratio) + ((1.96*mean(log_ratio))/sqrt(n)),
            log_ratio = mean(log_ratio),
  )
rp$taxa = str_replace(rp$taxa,pattern = "___"," / ")

rp = left_join(rp,dd)
rp$taxa = factor(rp$taxa,levels = rev(dd1$taxa))
rp$col = if_else(rp$p.adj.signif=="ns","black","red")
rp$Status = factor(rp$Status,levels = c("normal","abnormal"))



pdf(file = "Figures/caseStudy_fecal_lrMeans.pdf",width = 8 ,height = 5)
ggplot(rp,aes(taxa,log_ratio,fill  = Status,label = p.adj.signif))+
  geom_col(position = position_identity(),alpha = .5,col = "black",width = .5)+
  geom_errorbar(aes(ymin = lb,ymax = ub),width = .15,size = .5)+
  geom_point(aes(fill = Status),pch = 21,size = 1)+
  theme_bw()+
  coord_flip()+
  facet_grid(.~dx)+
  ggtitle("Log Ratio Means")+
  geom_text(nudge_y = 1,size = 2.5,fontface = "bold",y = -8,col = rp$col)+
  scale_y_continuous(limits = c(-9,10))+
  ylab("Denominator enriched <-- Log Ratio --> Numerator Enriched")+
  ggsci::scale_fill_lancet()+
  ggsci::scale_color_lancet()+
  theme(legend.position = "right",
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        strip.text = element_text(size = 7,hjust = 0.5,face = "bold"),strip.background = element_blank(),
        axis.title = element_text(size = 7,face = "bold"),legend.justification = "center",
        axis.title.y = element_blank(),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-5,-5,-5),
        axis.text = element_text(size = 7),legend.box.just = "top",
        legend.key.size = unit(.15,units = "in"),legend.spacing = unit(1,units = "in"),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
  )
dev.off()




# PLS Projection ----------------------------------------------------------
sc= glm.mdl.signature$finalModel$scores[,]
score.df = data.frame(md,sc)
score.df$level = factor(score.df$level,levels = c("normal","abnormal"))
pcv = glm.mdl.signature$finalModel$Xvar / glm.mdl.signature$finalModel$Xtotvar

pdf(file = "Figures/caseStudy_fecal_pls.pdf",width = 2.5 ,height = 1.7)
ggplot(score.df,aes(Comp.1,Comp.2,col = level))+
  geom_point(size = 1.5,alpha = .75)+
  #stat_ellipse()+
  theme_bw()+
  xlab(paste0("Comp.1 (",round(pcv[1],4)*100,"%)"))+
  ylab(paste0("Comp.2 (",round(pcv[2],4)*100,"%)"))+
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





# LR Network --------------------------------------------------------------
vi = varImp(glm.mdl.signature,useModel = T,scale = F)
vi  =data.frame(Ratio = rownames(vi$importance), Importance = vi$importance[,1])
vi = vi %>%
  arrange(desc(-Importance))
vi$Ratio = str_replace(vi$Ratio,pattern = "___"," / ")
imp.df = data.frame(Ratio = vi$Ratio,Imp = vi$Importance)
keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = " / ",remove = F)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g)$weight = imp.df$Imp
#g = minimum.spanning.tree(graph = g,weights = -E(g)$weight)
w_degree  = strength(g,mode = "total")
w_degree = data.frame(Taxa = names(w_degree),Str = as.numeric(w_degree))
w_degree = w_degree %>%
  arrange(desc(Str)) %>%
  top_n(n = 10,wt = Str)
plot(g,layout = layout_with_fr)
toGephi(g,"Output/fecalcal")






# Traditional Methods -----------------------------------------------------

blocks = factor(md$perm_blocks)


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

