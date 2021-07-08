rm(list = ls())
gc()

## Required Packages
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(foreach)



## Load Performance Data
fnames = dir("Results/featureSelection/")
fnames = paste0("Results/featureSelection/",fnames)
t1.df = data.frame()
l = length(fnames)
t1.df = foreach(i = 1:l,.combine = rbind)%do%{
  message(i)
  des= str_split(fnames[i],"-",simplify = T)[,2]
  ph = readr::read_csv(file = fnames[i])
  ph$design = des
  ph
}
dims = unique(t1.df$Dims)
sc = unique(t1.df$Scenario)
t1.df$CluseringCoefficient[is.na(t1.df$CluseringCoefficient)]=0

## Fix error in code from cluster runs where label was flipped
t1.df$design = if_else(t1.df$design=="Balanced","Unbalanced","Balanced")


# Compare Number of log rtios ----------------------------------------------------

ph.anova = t1.df %>%
  filter(design=="Balanced") %>%
  dplyr::select(design,c(colnames(t1.df)[1:7],"CluseringCoefficient","combinedF")) %>%
  dplyr::select(-totalRatios) %>%
  gather(key = "Metric","Value",6:9) %>%
  filter(Metric=="numRatios")
ph.anova$method = factor(ph.anova$method)
ph.anova$Scenario = factor(ph.anova$Scenario)
ph.anova$Scenario = factor(ph.anova$Scenario,label = paste("Scenario",1:5))
ph.anova$Dims = choose(ph.anova$Dims,2)


pdf(file = "Figures/featSelection_numRatios_balanced.pdf",width =7.5 ,height = 2.5)
ggplot(ph.anova,aes((Dims),(Value),fill = method,color = method))+
  stat_summary(fun = mean, geom = "line", aes(group = method,color = method),size=1)+
  stat_summary(fun.data = mean_cl_normal, aes(fill = method,color = method), geom = "pointrange",size = .5,pch = 23,)+
  scale_shape_manual(values = c(21,22,23,24,25))+
  facet_wrap(.~Scenario,nrow = 1,scales = "free")+
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(trans = "log10")+
  xlab("Log-ratio Dimensions")+
  ylab("# Log-Ratio Selected")+
  scale_color_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  scale_fill_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  theme(legend.position = "top",legend.margin = margin(0,0,0,0),legend.title = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 10),
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.box.margin = margin(t=0, r=0, b=-0.5, l=0, unit="cm"),
        #plot.margin = unit(x = c(-.2, .2, .2, .2), units = "cm"),
        axis.title.y =element_text(face = "bold",size = 10),
        axis.title.x =element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
  )

dev.off()


# Compare Association Strength ----------------------------------------------------

ph.anova = t1.df %>%
  filter(design=="Balanced") %>%
  dplyr::select(design,c(colnames(t1.df)[1:7],"CluseringCoefficient","combinedF")) %>%
  dplyr::select(-totalRatios) %>%
  gather(key = "Metric","Value",6:9) %>%
  filter(Metric=="combinedF")
ph.anova$method = factor(ph.anova$method)
ph.anova$Scenario = factor(ph.anova$Scenario)
ph.anova$Scenario = factor(ph.anova$Scenario,label = paste("Scenario",1:5))
ph.anova$Dims = choose(ph.anova$Dims,2)


pdf(file = "Figures/featSelection_AssocStr_balanced.pdf",width =7.5 ,height = 2.5)
# tiff(filename = "Figures/featSel_ovr_balanced_numLogRatios.tiff",width = 9,height = 2.5,units = "in",res = 300)
ggplot(ph.anova,aes((Dims),(Value),fill = method,color = method))+
  stat_summary(fun = mean, geom = "line", aes(group = method,color = method),size=1)+
  stat_summary(fun.data = mean_cl_normal, aes(fill = method,color = method), geom = "pointrange",size = .5,pch = 23)+
  scale_shape_manual(values = c(21,22,23,24,25))+
  facet_wrap(.~Scenario,nrow = 1,scales = "free")+
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab("Log-ratio Dimensions")+
  ylab("cF Statistic")+
  #annotation_logticks(sides = "l")+
  scale_color_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  scale_fill_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  theme(legend.position = "top",legend.margin = margin(0,0,0,0),legend.title = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 10),
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.box.margin = margin(t=0, r=0, b=-0.5, l=0, unit="cm"),
        #plot.margin = unit(x = c(-.2, .2, .2, .2), units = "cm"),
        axis.title.y =element_text(face = "bold",size = 10),
        axis.title.x =element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
  )

dev.off()




# Clustering Coefificents -------------------------------------------------

ph.anova = t1.df %>%
  filter(design=="Balanced") %>%
  dplyr::select(design,c(colnames(t1.df)[1:7],"CluseringCoefficient","combinedF")) %>%
  dplyr::select(-totalRatios) %>%
  gather(key = "Metric","Value",6:9) %>%
  filter(Metric=="CluseringCoefficient")
ph.anova$method = factor(ph.anova$method)
ph.anova$Scenario = factor(ph.anova$Scenario)
ph.anova$Scenario = factor(ph.anova$Scenario,label = paste("Scenario",1:5))
ph.anova$Dims = choose(ph.anova$Dims,2)

pdf(file = "Figures/featSelection_ClusCoef_balanced.pdf",width =7.5 ,height = 2.5)
# tiff(filename = "Figures/featSel_ovr_balanced_numLogRatios.tiff",width = 9,height = 2.5,units = "in",res = 300)
ggplot(ph.anova,aes((Dims),(Value),fill = method,color = method))+
  stat_summary(fun = mean, geom = "line", aes(group = method,color = method),size=1)+
  stat_summary(fun.data = mean_cl_normal, aes(fill = method,color = method), geom = "pointrange",size = .5,pch = 23)+
  scale_shape_manual(values = c(21,22,23,24,25))+
  facet_wrap(.~Scenario,nrow = 1,scales = "free")+
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #               labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  #scale_y_continuous(trans = "log10")+
  xlab("Log-ratio Dimensions")+
  ylab("Clustering Coeff.")+
  #annotation_logticks(sides = "l")+
  scale_color_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  scale_fill_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  theme(legend.position = "top",legend.margin = margin(0,0,0,0),legend.title = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 10),
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.box.margin = margin(t=0, r=0, b=-0.5, l=0, unit="cm"),
        #plot.margin = unit(x = c(-.2, .2, .2, .2), units = "cm"),
        axis.title.y =element_text(face = "bold",size = 10),
        axis.title.x =element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
  )

dev.off()



# Computational Time -------------------------------------------------

ph.anova = t1.df %>%
  filter(design=="Balanced") %>%
  dplyr::select(design,c(colnames(t1.df)[1:7],"CluseringCoefficient","combinedF")) %>%
  dplyr::select(-totalRatios) %>%
  gather(key = "Metric","Value",6:9) %>%
  filter(Metric=="compTime_sec")
ph.anova$method = factor(ph.anova$method)
ph.anova$Scenario = factor(ph.anova$Scenario)
ph.anova$Scenario = factor(ph.anova$Scenario,label = paste("Scenario",1:5))
ph.anova$Dims = choose(ph.anova$Dims,2)

pdf(file = "Figures/featSel_comptimeBalanaced.pdf",width =7.5 ,height =2.5 )
#tiff(filename = "Figures/.5featSel_ovr_balanced_numLogRatios.tiff",width = 9,height = 2.5,units = "in",res = 300)
ggplot(ph.anova,aes((Dims),(Value),fill = method,color = method))+
  stat_summary(fun = mean, geom = "line", aes(group = method,color = method),size=1)+
  stat_summary(fun.data = mean_cl_normal, aes(fill = method,color = method), geom = "pointrange",size = .5,pch = 23)+
  scale_shape_manual(values = c(21,22,23,24,25))+
  facet_wrap(.~Scenario,nrow = 1,scales = "free")+
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #               labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(trans = "log2")+
  xlab("Log-ratio Dimensions")+
  ylab("Computational Time (sec.)")+
  #annotation_logticks(sides = "l")+
  scale_color_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  scale_fill_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  theme(legend.position = "top",legend.margin = margin(0,0,0,0),legend.title = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 10),
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.box.margin = margin(t=0, r=0, b=-0.5, l=0, unit="cm"),
        #plot.margin = unit(x = c(-.2, .2, .2, .2), units = "cm"),
        axis.title.y =element_text(face = "bold",size = 10),
        axis.title.x =element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
  )


dev.off()







# Compare Number of log rtios ----------------------------------------------------

ph.anova = t1.df %>%
  filter(design=="Unbalanced") %>%
  dplyr::select(design,c(colnames(t1.df)[1:7],"CluseringCoefficient","combinedF")) %>%
  dplyr::select(-totalRatios) %>%
  gather(key = "Metric","Value",6:9) %>%
  filter(Metric=="numRatios")
ph.anova$method = factor(ph.anova$method)
ph.anova$Scenario = factor(ph.anova$Scenario)
ph.anova$Scenario = factor(ph.anova$Scenario,label = paste("Scenario",1:5))
ph.anova$Dims = choose(ph.anova$Dims,2)


pdf(file = "Figures/featSelection_numRatios_Unbalanced.pdf",width =7.5 ,height = 2.5)
ggplot(ph.anova,aes((Dims),(Value),fill = method,color = method))+
  stat_summary(fun = mean, geom = "line", aes(group = method,color = method),size=1)+
  stat_summary(fun.data = mean_cl_normal, aes(fill = method,color = method), geom = "pointrange",size = .5,pch = 23,)+
  scale_shape_manual(values = c(21,22,23,24,25))+
  facet_wrap(.~Scenario,nrow = 1,scales = "free")+
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(trans = "log10")+
  xlab("Log-ratio Dimensions")+
  ylab("# Log-Ratio Selected")+
  scale_color_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  scale_fill_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  theme(legend.position = "top",legend.margin = margin(0,0,0,0),legend.title = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 10),
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.box.margin = margin(t=0, r=0, b=-0.5, l=0, unit="cm"),
        #plot.margin = unit(x = c(-.2, .2, .2, .2), units = "cm"),
        axis.title.y =element_text(face = "bold",size = 10),
        axis.title.x =element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
  )

dev.off()


# Compare Association Strength ----------------------------------------------------

ph.anova = t1.df %>%
  filter(design=="Unbalanced") %>%
  dplyr::select(design,c(colnames(t1.df)[1:7],"CluseringCoefficient","combinedF")) %>%
  dplyr::select(-totalRatios) %>%
  gather(key = "Metric","Value",6:9) %>%
  filter(Metric=="combinedF")
ph.anova$method = factor(ph.anova$method)
ph.anova$Scenario = factor(ph.anova$Scenario)
ph.anova$Scenario = factor(ph.anova$Scenario,label = paste("Scenario",1:5))
ph.anova$Dims = choose(ph.anova$Dims,2)


pdf(file = "Figures/featSelection_AssocStr_Unbalanced.pdf",width =7.5 ,height = 2.5)
# tiff(filename = "Figures/featSel_ovr_balanced_numLogRatios.tiff",width = 9,height = 2.5,units = "in",res = 300)
ggplot(ph.anova,aes((Dims),(Value),fill = method,color = method))+
  stat_summary(fun = mean, geom = "line", aes(group = method,color = method),size=1)+
  stat_summary(fun.data = mean_cl_normal, aes(fill = method,color = method), geom = "pointrange",size = .5,pch = 23)+
  scale_shape_manual(values = c(21,22,23,24,25))+
  facet_wrap(.~Scenario,nrow = 1,scales = "free")+
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab("Log-ratio Dimensions")+
  ylab("cF Statistic")+
  #annotation_logticks(sides = "l")+
  scale_color_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  scale_fill_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  theme(legend.position = "top",legend.margin = margin(0,0,0,0),legend.title = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 10),
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.box.margin = margin(t=0, r=0, b=-0.5, l=0, unit="cm"),
        #plot.margin = unit(x = c(-.2, .2, .2, .2), units = "cm"),
        axis.title.y =element_text(face = "bold",size = 10),
        axis.title.x =element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
  )

dev.off()




# Clustering Coefificents -------------------------------------------------

ph.anova = t1.df %>%
  filter(design=="Unbalanced") %>%
  dplyr::select(design,c(colnames(t1.df)[1:7],"CluseringCoefficient","combinedF")) %>%
  dplyr::select(-totalRatios) %>%
  gather(key = "Metric","Value",6:9) %>%
  filter(Metric=="CluseringCoefficient")
ph.anova$method = factor(ph.anova$method)
ph.anova$Scenario = factor(ph.anova$Scenario)
ph.anova$Scenario = factor(ph.anova$Scenario,label = paste("Scenario",1:5))
ph.anova$Dims = choose(ph.anova$Dims,2)

pdf(file = "Figures/featSelection_ClusCoef_Unbalanced.pdf",width =7.5 ,height = 2.5)
# tiff(filename = "Figures/featSel_ovr_balanced_numLogRatios.tiff",width = 9,height = 2.5,units = "in",res = 300)
ggplot(ph.anova,aes((Dims),(Value),fill = method,color = method))+
  stat_summary(fun = mean, geom = "line", aes(group = method,color = method),size=1)+
  stat_summary(fun.data = mean_cl_normal, aes(fill = method,color = method), geom = "pointrange",size = .5,pch = 23)+
  scale_shape_manual(values = c(21,22,23,24,25))+
  facet_wrap(.~Scenario,nrow = 1,scales = "free")+
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #               labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  #scale_y_continuous(trans = "log10")+
  xlab("Log-ratio Dimensions")+
  ylab("Clustering Coeff.")+
  #annotation_logticks(sides = "l")+
  scale_color_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  scale_fill_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  theme(legend.position = "top",legend.margin = margin(0,0,0,0),legend.title = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 10),
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.box.margin = margin(t=0, r=0, b=-0.5, l=0, unit="cm"),
        #plot.margin = unit(x = c(-.2, .2, .2, .2), units = "cm"),
        axis.title.y =element_text(face = "bold",size = 10),
        axis.title.x =element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
  )

dev.off()



# Computational Time -------------------------------------------------

ph.anova = t1.df %>%
  filter(design=="Unbalanced") %>%
  dplyr::select(design,c(colnames(t1.df)[1:7],"CluseringCoefficient","combinedF")) %>%
  dplyr::select(-totalRatios) %>%
  gather(key = "Metric","Value",6:9) %>%
  filter(Metric=="compTime_sec")
ph.anova$method = factor(ph.anova$method)
ph.anova$Scenario = factor(ph.anova$Scenario)
ph.anova$Scenario = factor(ph.anova$Scenario,label = paste("Scenario",1:5))
ph.anova$Dims = choose(ph.anova$Dims,2)

pdf(file = "Figures/featSel_Unbalanced_comptime.pdf",width =7.5 ,height =2.5 )
#tiff(filename = "Figures/.5featSel_ovr_balanced_numLogRatios.tiff",width = 9,height = 2.5,units = "in",res = 300)
ggplot(ph.anova,aes((Dims),(Value),fill = method,color = method))+
  stat_summary(fun = mean, geom = "line", aes(group = method,color = method),size=1)+
  stat_summary(fun.data = mean_cl_normal, aes(fill = method,color = method), geom = "pointrange",size = .5,pch = 23)+
  scale_shape_manual(values = c(21,22,23,24,25))+
  facet_wrap(.~Scenario,nrow = 1,scales = "free")+
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #               labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(trans = "log2")+
  xlab("Log-ratio Dimensions")+
  ylab("Computational Time (sec.)")+
  #annotation_logticks(sides = "l")+
  scale_color_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  scale_fill_manual(values =  c(ggsci::pal_jco()(5)[2:5],"blue"))+
  theme(legend.position = "top",legend.margin = margin(0,0,0,0),legend.title = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 10),
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.box.margin = margin(t=0, r=0, b=-0.5, l=0, unit="cm"),
        #plot.margin = unit(x = c(-.2, .2, .2, .2), units = "cm"),
        axis.title.y =element_text(face = "bold",size = 10),
        axis.title.x =element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
  )


dev.off()







