rm(list = ls())
gc()

## Required Packages
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(foreach)




# Exp Data Results ------------------------------------------------------------


## Not Run unless needed

# # ## Extract Scenarios
# scenarios = dir("Results/powerExpData/")
#
# ## Read data
# df = data.frame()
# for(j in 1:length(scenarios)){
#
#   path_ = paste0("Results/powerExpData/",scenarios[j])
#
#   fnames = dir(path_)
#   s1 = data.frame()
#   l = length(fnames)
#   s1 = foreach(i = 1:l,.combine = rbind)%do%{
#     message(i)
#     readr::read_csv(file = paste0(path_,"/",fnames[i]))
#   }
#   df = rbind(df,s1)
# }
#
# ## Store Results
# write_csv(df,"Results/Merged Results/expData_results.csv")




## Load Files Merged Results
df = read_csv("Results/Merged Results/expData_results.csv")

df = df %>%
  mutate(methStat = paste0(test,"_",Statistic))

## Define Statistics to keep for anlaysis
keepStats = c("ANOSIM_R",
              "Energy_Aitch_ALL_estat" ,
              "PERMANOVA_permanovaF",
              "PermDisp2_permDispF",
              "selProPerm_testStat"
              )
alpha = 0.05
df = df %>%
  mutate(reject = if_else(p<alpha,1,0)) %>%
  filter(methStat %in% keepStats)
nameMap = data.frame(methStat = keepStats,Method = c("ANOSIM",
                                                     "Energy",
                                                     "PERMANOVA",
                                                     "PERMDISP2",
                                                     "selEnergyPerm"
))

df = left_join(df,nameMap) %>%
  mutate(methStat = Method) %>%
  dplyr::select(-Method)



## Check for missing data
str(df)
status = df %>%
  group_by(methStat,scenario,sampleDesign,shiftParm,Seed) %>%
  summarise(n = n())
status = spread(status,"Seed","n")
bool = is.na(rowSums(status[,5:ncol(status)]))
status = status[bool,]
status = data.frame(status[,1:4],numMissing =  rowSums(is.na(status[,5:ncol(status)])))


## Sim Scenarios
test_ = unique(df$test)
simScenarios = unique(df$scenario)
sampleDesign = unique(df$sampleDesign)
shiftParm = unique(df$shiftParm)
df$binaryVote = if_else(df$p<=alpha,"different","same")

auc.df = data.frame()
BS_mean = F
bootstrapReps = 1
innerLoop = 1
for(r in 1:bootstrapReps){
  set.seed(r)
  for(sd in sampleDesign){
    for(ss in simScenarios){
      for(sp  in c(1,2,3,4)){
        for(l in 1:innerLoop){
          ph = df %>%
            filter(scenario == ss) %>%
            filter(sampleDesign == sd) %>%
            filter(shiftParm == sp)
          ph$trueDistr = factor(ph$trueDistr,levels = c("same","different"))

          if(length(unique(ph$trueDistr))>1){
            rr = length(unique(ph$Seed))
            tt = ph %>%
              filter(methStat=="selProPerm: testStat")
            numTrue = sum(tt$trueDistr=="same")
            numDiff = sum(tt$trueDistr=="different")
            seeds = sample(c(unique(ph$Seed)),size = length(unique(ph$Seed)),replace = BS_mean)
            ph = ph %>%
              filter(Seed %in% seeds) %>%
              group_by(methStat) %>%
              summarise(n = n(), num_same = numTrue,num_diff = numDiff,auc =as.numeric( pROC::auc(trueDistr,p)),
                        th = pROC::coords(pROC::roc(trueDistr,p,levels = levels(ph$trueDistr )),x=0.05, transpose = FALSE,input = "threshold",ret = c("tp","tn","fp","fn","fpr",
                                                                                                                       "fnr","ppv","npv","fnr",
                                                                                                                       "specificity",
                                                                                                                       "sensitivity","youden"))[1,],
                        MCC = mltools::mcc( actuals = factor(trueDistr),preds = factor(binaryVote,levels = c("different","same"))),
                        Sen = caret::confusionMatrix(factor(trueDistr),factor(binaryVote,levels = c("different","same")))$byClass[1],
                        Spe = caret::confusionMatrix(factor(trueDistr),factor(binaryVote,levels = c("different","same")))$byClass[2],
                        PPV = caret::confusionMatrix(factor(trueDistr),factor(binaryVote,levels = c("different","same")))$byClass[3],
                        NPV = caret::confusionMatrix(factor(trueDistr),factor(binaryVote,levels = c("different","same")))$byClass[4],
                        F1 = caret::confusionMatrix(factor(trueDistr),factor(binaryVote,levels = c("different","same")))$byClass[7],
                        posClass = caret::confusionMatrix(factor(trueDistr),factor(binaryVote,levels = c("different","same")))$positive,

              )
            phh = data.frame(Deisgn = sd,Scenario = ss,shiftParm = sp,ph)
            phh = cbind(phh,ph$th) %>%
              dplyr::select(-th)
            phh = phh %>%
              group_by(posClass,Deisgn,Scenario,shiftParm,methStat) %>%
              summarise_all(.funs = mean)
            phh$BootstrapRep = r
            auc.df = rbind(auc.df,phh)
          }
        }
      }
    }
  }
}






# Process Bootstrap Data --------------------------------------------------

auc.df[is.na(auc.df)] = 0
auc.df$t1Error = auc.df$fp/auc.df$n
auc.df$Youden = auc.df$Sen + auc.df$Spe - 1





# Overall Metrics ---------------------------------------------------------
pp = auc.df %>%
  select(Deisgn,shiftParm,methStat,MCC,sensitivity,specificity,ppv,npv,Youden,fpr)
pp =gather(pp,key = "metric","value",6:ncol(pp))

pp = pp %>%
  filter(methStat!="PERMDISP2")
pp$metric = factor(pp$metric,levels = unique(pp$metric),labels = c("MCC","Sensitivity","Specificity","PPV","NPV","Youden","FPR"))


pdf(file = "Figures/Additional File X.pdf",width =5 ,height = 7)
ggplot(pp,aes(methStat,value,fill = methStat))+
  stat_summary(fun = mean, geom = "col", aes(group = methStat,fill = methStat),size=.5,col = "black",width = .75)+
  stat_summary(fun.data = mean_se, aes(fill = methStat), geom = "errorbar",size = .5,width = .2)+
  #stat_summary(fun = mean, geom = "point", aes(group = methStat,fill = methStat),size=1,pch = 21)+
  coord_flip()+
  ggsci::scale_fill_futurama()+
  facet_grid(metric~Deisgn)+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(face = "bold",size = 8),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        axis.title.x = element_blank(),
        legend.margin=margin(-10,-10,-10,-10),
        strip.text = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
        strip.background = element_blank(),
        axis.title.y = element_blank())

dev.off()




# MCC ---------------------------------------------------------------------

auc.df1 = auc.df %>%
  group_by(methStat,Deisgn,shiftParm,Scenario) %>%
  summarise(mcc = mean(MCC),n  = n(),
            lb = quantile(MCC,probs = .025),
            ub = quantile(MCC,probs = .975)
  ) %>%
  mutate(dimensions = choose(shiftParm,2))

ovr = auc.df1 %>%
  group_by(methStat) %>%
  summarise_all(.funs = mean) %>%
  arrange(desc(-mcc))
auc.df1$methStat = factor(auc.df1$methStat,levels = ovr$methStat)

keepStats1  = c("ANOSIM",
                "Energy",
                "PERMANOVA",
                "PERMDISP2",
                "selEnergyPerm"
)
methStatOrder  = c("PERMANOVA","Energy","ANOSIM","PERMDISP2","selEnergyPerm")
pp = auc.df1 %>%
  filter(methStat%in% keepStats1) %>%
  mutate(methStat = factor(methStat,levels = methStatOrder))
parm = 1-seq(.5,.95,length.out = 4)[shiftParm];
cc = seq(1,4,length.out = 4)
ppp  = pp %>%
  group_by(Scenario,shiftParm) %>%
  summarise(n = n())
ppp$TrueParm = rep(c(cc,parm),2)
pp = left_join(pp,ppp,by = c("shiftParm", "Scenario"))
pp$Scenario = factor(pp$Scenario,
                      labels =  c("16S: Increasing Covariance Diff.",
                                  "16S: Increasing Signal Density",
                                  "WGS: Increasing Covariance Diff.",
                                  "WGS: Increasing Signal Density" ))
pp$Deisgn = factor(pp$Deisgn,
                   labels =  c("Balanced (n1= n2)","Unbalanced (n1 < n2)"))


pdf(file = "Figures/fig4.pdf",width =7.5 ,height = 5)
ggplot(pp,aes(TrueParm,mcc))+
  geom_point(aes(shape=methStat,col = methStat,fill = methStat),size = 2,alpha = 1)+
  geom_line(aes(group = methStat,col = methStat,lty = methStat),size = .75,alpha = .75)+
  geom_linerange(aes(ymin = lb,ymax =ub,col = methStat ),size = .5)+
  scale_fill_manual(values = c(ggsci::pal_d3()(3),"darkred","blue"))+
  scale_color_manual(values = c(ggsci::pal_d3()(3),"darkred","blue"))+
  scale_shape_manual(values = c(24,23,21,22,25))+
  scale_linetype_manual(values = c("solid","solid","solid","dashed","solid"))+
  scale_y_continuous(limits = c(-.5,1))+
  facet_grid(Deisgn~Scenario,scales = "free")+
  theme_bw()+
  xlab("Effect Level")+
  ylab(expression(paste('Matthews Correlation Coefficient (', alpha," = 0.05)")))+
  theme(legend.position = "top")+
  theme(legend.position = "top",
        #panel.grid = element_blank(),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(face = "bold",size = 8),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-10,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        #legend.key.size = unit(.05,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.background = element_rect(colour = "black")
  )+
  theme(legend.position = "top",
        strip.text = element_text(face = "bold",size = 7),
        strip.background = element_blank(),
        legend.title = element_blank(),
        #axis.title.y = element_blank(),
        #axis.text.y = element_blank()
  )
dev.off()










