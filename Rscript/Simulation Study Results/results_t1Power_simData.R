rm(list = ls())
gc()

## Required Packages
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(foreach)




## Not Run unless needed ; initial merge of results from cluster

# # Null Results ------------------------------------------------------------
#
# ## Extract Scenarios
# scenarios = dir("Results/powerSimData_Null")
#
# ## Read data
# df = data.frame()
# for(j in 1:length(scenarios)){
#
#   path_ = paste0("Results/powerSimData_Null/",scenarios[j])
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
# ## PERMAOVA Correction
# pmv_ = read_csv(file = "Results/permanovaCorrection_null.csv")
#
# ## Remove old PERMANOVA data
# df = df %>%
#   filter(test!="PERMANOVA")
#
# ## APpend corrected data
# df = rbind(df,pmv_)
# df = df %>%
#   mutate(methStat = paste0(test,"_",Statistic))
# unique(df$methStat)
#
# ## List of Statistics to keep for the final analysis
# keepStats = c("ANOSIM_R",
#               "Energy_Aitch_ALL_estat" ,
#               "PERMANOVA_permanovaF",
#               "PermDisp2_permDispF",
#               "selProPerm_testStat"
# )
#
# ## Define Alpha Significance
# alpha = 0.05
# df = df %>%
#   mutate(reject = if_else(p<alpha,1,0)) %>%
#   filter(methStat %in% keepStats)
# nameMap = data.frame(methStat = keepStats,Method = c("ANOSIM",
#                                                      "Energy",
#                                                      "PERMANOVA",
#                                                      "PERMDISP2",
#                                                      "selEnergyPerm"
# ))
# df = left_join(df,nameMap) %>%
#   mutate(methStat = Method) %>%
#   dplyr::select(-Method)
# df$Distr = "same"
#
# ## Store Null Results (Results where there were no true differences as ground truth)
# null.df = df
#
#
#
# # True Difference Results -------------------------------------------------
#
# scenarios = dir("Results/powerSimData/")
#
# ## Read data
# df = data.frame()
# for(j in 1:length(scenarios)){
#
#   path_ = paste0("Results/powerSimData/",scenarios[j])
#   fnames = dir(path_)
#   if(length(fnames>0)){
#     s1 = data.frame()
#     l = length(fnames)
#     s1 = foreach(i = 1:l,.combine = rbind)%do%{
#       message(i)
#       readr::read_csv(file = paste0(path_,"/",fnames[i]))    }
#     df = rbind(df,s1)
#   }
#
# }
#
# ## PERMAOVA Correction
# pmv = read_csv(file = "Results/permanovaCorrection.csv")
#
# ## Remove old PERMANOVA data
# df = df %>%
#   filter(test!="PERMANOVA")
#
# ## Append corrected data
# df = rbind(df,pmv)
# df = df %>%
#   mutate(methStat = paste0(test,"_",Statistic))
# alpha = 0.05
# df = df %>%
#   mutate(reject = if_else(p<alpha,1,0)) %>%
#   filter(methStat %in% keepStats)
# df = left_join(df,nameMap) %>%
#   mutate(methStat = Method) %>%
#   dplyr::select(-Method)
# df$Distr = "different"
#
# ## merge true and null data
# df = rbind(df,null.df)
#
# ## Write Merged Results to single File for faster loading
# write_csv(df,"Results/Merged Results/simData_results.csv")



df = read_csv("Results/Merged Results/simData_results.csv")

## Sim Scenarios
simScenarios = unique(df$scenario)
sampleDesign = unique(df$sampleDesign)
dims = unique(df$dimensions)
alpha = 0.05
df$binaryVote = if_else(df$p<=alpha,"different","same")


## Compute Mean with bootstrap confidence interval
auc.df = data.frame()
BS_mean = F
bootstrapReps = 1
innerLoop = 1
for(r in bootstrapReps){
  set.seed(r)
  for(sd in sampleDesign){
    for(ss in simScenarios){
      for(sp  in dims){
        for(l in 1:innerLoop){
          ph = df %>%
            filter(scenario == ss) %>%
            filter(sampleDesign == sd) %>%
            filter(dimensions == sp)
          ph$Distr = factor(ph$Distr,levels = c("same","different"))
          if(length(unique(ph$Distr))>1){
            seeds = sample(c(unique(ph$Seed)),size = length(unique(ph$Seed)),replace = BS_mean)
            ph = ph %>%
              filter(Seed %in% seeds) %>%
              group_by(methStat) %>%
              summarise(n = n(),auc =as.numeric( pROC::auc(Distr,p)),
                        th = pROC::coords(pROC::roc(Distr,p,levels = levels(ph$Distr)),x=0.05, transpose = FALSE,input = "threshold",ret = c("tp","tn","fp","fn","fpr",
                                                                                                                       "fnr","ppv","npv","fnr",
                                                                                                                       "specificity",
                                                                                                                       "sensitivity","youden"))[1,],
                        MCC = mltools::mcc( actuals = factor(Distr),preds = factor(binaryVote,levels = c("different","same"))),
                        Sen = caret::confusionMatrix(factor(Distr),factor(binaryVote,levels = c("different","same")))$byClass[1],
                        Spe = caret::confusionMatrix(factor(Distr),factor(binaryVote,levels = c("different","same")))$byClass[2],
                        PPV = caret::confusionMatrix(factor(Distr),factor(binaryVote,levels = c("different","same")))$byClass[3],
                        NPV = caret::confusionMatrix(factor(Distr),factor(binaryVote,levels = c("different","same")))$byClass[4],
                        F1 = caret::confusionMatrix(factor(Distr),factor(binaryVote,levels = c("different","same")))$byClass[7],
                        posClass = caret::confusionMatrix(factor(Distr),factor(binaryVote,levels = c("different","same")))$positive,

              )
            phh = data.frame(Deisgn = sd,Scenario = ss,dimensions = sp,ph)
            phh = cbind(phh,ph$th) %>%
              dplyr::select(-th)
          }
        }
        phh = phh %>%
          group_by(Deisgn,Scenario,dimensions,methStat) %>%
          summarise_all(.funs = mean)
        phh$BootstrapRep = r
        auc.df = rbind(auc.df,phh)
        message(r)
      }
    }
  }
}





## Process Boostrap Data
auc.df[is.na(auc.df)] = 0
auc.df$t1Error = auc.df$fp/auc.df$n
auc.df$Youden = auc.df$Sen + auc.df$Spe - 1





# Overall Metrics ---------------------------------------------------------
pp = auc.df %>%
  select(Deisgn,dimensions,methStat,MCC,sensitivity,specificity,ppv,npv,Youden,fpr)
pp =gather(pp,key = "metric","value",5:ncol(pp))

pp = pp %>%
  filter(methStat!="PERMDISP2")
pp$metric = factor(pp$metric,levels = unique(pp$metric),labels = c("MCC","Sensitivity","Specificity","PPV","NPV","Youden","FPR"))


pdf(file = "Figures/Additional File X1.pdf",width =5 ,height = 7)
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
  group_by(methStat,Deisgn,dimensions,Scenario) %>%
  summarise(mcc = mean(MCC),n  = n(),
            lb = quantile(MCC,probs = .025),
            ub = quantile(MCC,probs = .975)
  ) %>%
  mutate(dimensions = choose(dimensions,2))

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

pp = auc.df1 %>%
  filter(methStat%in% keepStats1)
pp$dimensions = factor(pp$dimensions,
                       labels =  c("lr = 1,225","lr = 11,175","lr = 31,125"))


pp = auc.df1 %>%
  filter(methStat%in% keepStats1) %>%
  mutate(methStat = factor(methStat,levels = ovr$methStat)) %>%
  mutate(label = if_else(methStat=="selEnergyPerm",3,2))
pp$Scenario = factor(pp$Scenario,
                       labels =  c("Scenario 1","Scenario 2","Scenario 3","Scenario 4"))
pp$Deisgn = factor(pp$Deisgn,
                     labels =  c("Balanced (n1= n2)","Unbalanced (n1 < n2)"))
pp$methStat = factor(pp$methStat,levels = c("PERMANOVA","Energy","ANOSIM","PERMDISP2","selEnergyPerm"))

pdf(file = "Figures/fig3.pdf",width =7.5 ,height = 5)
ggplot(pp,aes(dimensions,mcc,col = methStat,group = methStat))+
  geom_point(aes(shape=methStat,col = methStat,fill = methStat),size = 2,alpha = 1)+
  geom_line(aes(group = methStat,col = methStat,lty = methStat),size = .75,alpha = .75)+
  geom_errorbar(aes(ymin = lb,ymax =ub,col = methStat ),width = .07,size = .5)+
  scale_fill_manual(values = c(ggsci::pal_d3()(3),"darkred","blue"))+
  scale_color_manual(values = c(ggsci::pal_d3()(3),"darkred","blue"))+
  scale_shape_manual(values = c(24,25,21,22,23))+
  scale_linetype_manual(values = c("solid","solid","solid","dashed","solid"))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(limits = c(-.5,1))+
  facet_grid(Deisgn~Scenario)+
  theme_bw()+
  xlab("Number of Log Ratios")+
  ylab(expression(paste('Matthews Correlation Coefficient ( ' , alpha ," = 0.05 )")))+
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
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.title = element_blank()
        )
dev.off()


