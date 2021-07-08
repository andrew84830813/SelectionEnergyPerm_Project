

sparsePercent = 0.85



allPC = data.frame()

for(p in c(2,3,4)){
  for(s in 1:4){
    
    seed_ = 1# random seed selection
    sampleDesign = 1; design = if_else(sampleDesign==1,"Balanced","Unbalanced")
    shiftParm = p
    scenario = s
    label_parm = 2
    
    ###----------------------------------------*
    ## Power  Analysis Simulation Scenarios
    ###----------------------------------------*
    if(sampleDesign==1){
      g1 = 40
      g2 = 40
    }else{
      g1 = 20
      g2 = 60
    }
    
    
    switch(scenario,
           {
             ## Scenario-1 = Sim from 16S Data with mean shift
             mdwgs = readRDS("Output/16sModel.RDS");
             set.seed(seed_);
             dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
             dat = sample_n(dat,size = g1+g2,replace = F);
             labels = sample(c(rep("S1",g1),rep("S2",g2)));
             dat = data.frame(Status = labels,dat);
             fname = "16S_meanShift";
             path_ = "/nas/longleaf/home/andrew84/selPermEnergy/results2/powerExpData/s1/";
             
             if(label_parm==2){
               #process shift
               procData = processCompData(dat,minPrevalence = sparsePercent);
               dat = procData$processedData;
               impFact = procData$impFactor;
               y = dat[,-1];
               bool = colSums(y)==0;
               y = y[,!bool];
               dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));

               parm = seq(1,1.75,length.out = 4)[shiftParm];
               dat = simFromExpData.largeMeanShft(raMatrix = dat[,-1],n1 = g1,n2 = g2,featureShiftPercent =  parm, perFixedFeatures = .95)
             }
             
           },
           
           {
             ## Scenario-2 = Sim from 16S Data with cov shift
             mdwgs = readRDS("Output/16sModel.RDS");
             set.seed(seed_);
             dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
             dat = sample_n(dat,size = g1+g2,replace = F);
             labels = sample(c(rep("S1",g1),rep("S2",g2)));
             dat = data.frame(Status = labels,dat);
             fname = "16S_covShift";
             path_ = "/nas/longleaf/home/andrew84/selPermEnergy/results2/powerExpData/s2/";
             
             if(label_parm==2){
               procData = processCompData(dat,minPrevalence = sparsePercent);
               dat = procData$processedData;
               impFact = procData$impFactor;
               y = dat[,-1];
               bool = colSums(y)==0;
               y = y[,!bool];
               dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));

               parm = seq(1,8,length.out = 4)[shiftParm];
               dat = simFromExpData.covarianceShift(raMatrix = dat[,-1],n1 = g1,n2 = g2,maxCov =parm) 
             }
             
           },
           
           {
             ## Scenario-3 = Sim from WGS Data with mean shift
             mdwgs = readRDS("Output/wgsModel.RDS");
             set.seed(seed_);
             dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
             dat = sample_n(dat,size = g1+g2);
             labels = sample(c(rep("S1",g1),rep("S2",g2)));
             dat = data.frame(Status = labels,dat);
             fname = "WGS_meanShift";
             path_ = "/nas/longleaf/home/andrew84/selPermEnergy/results2/powerExpData/s3/";
             
             if(label_parm==2){
               procData = processCompData(dat,minPrevalence = sparsePercent);
               dat = procData$processedData;
               impFact = procData$impFactor;
               y = dat[,-1];
               bool = colSums(y)==0;
               y = y[,!bool];
               dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));

               parm = seq(1,1.75,length.out = 4)[shiftParm];
               dat = simFromExpData.largeMeanShft(raMatrix = dat[,-1],n1 = g1,n2 = g2,featureShiftPercent =  parm, perFixedFeatures = .95)
             }
           },
           
           
           {
             ## Scenario-4 = Sim from WGS Data with cov shift
             mdwgs = readRDS("Output/wgsModel.RDS");
             set.seed(seed_);
             dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
             dat = sample_n(dat,size = g1+g2);
             labels = sample(c(rep("S1",g1),rep("S2",g2)));
             dat = data.frame(Status = labels,dat);
             fname = "WGS_covShift";
             path_ = "/nas/longleaf/home/andrew84/selPermEnergy/results2/powerExpData/s4/";
             
             if(label_parm==2){
               procData = processCompData(dat,minPrevalence = sparsePercent);
               dat = procData$processedData;
               impFact = procData$impFactor;
               y = dat[,-1];
               bool = colSums(y)==0;
               y = y[,!bool];
               dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));
               parm = seq(1,10,length.out = 4)[shiftParm];
               dat = simFromExpData.covarianceShift(raMatrix = dat[,-1],n1 = g1,n2 = g2,maxCov =parm) 
             }
             
           }
           
    )
    
    
    
    message("seed",seed_,"-",design,"-dims",dims,"-",fname)
    
    
    ## Sparisty process
    procData = processCompData(dat,minPrevalence = .9)
    dat = procData$processedData
    impFact = procData$impFactor
    minorityClass = procData$minClss
    majorityClass = procData$majClass
    impClass = minorityClass
    Scenario = fname
    message("sparisty process Complete")
    
    
    ## Pre Process -  Removbe col with all 0's
    y = dat[,-1]
    bool = colSums(y)==0
    y = y[,!bool]
    dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact))
    message("preProcess Complete")
    
    ## Compute Ratios
    dat.f = data.frame(Status = factor(dat[,1]), fastImputeZeroes( compositions::clo(dat[,-1]) , impFactor  = impFact ) )
    lrs = calcLogRatio(dat.f)
    message("logratios computed")
    pc = prcomp(lrs[,-1],center = T,scale. = T)#Rtsne::Rtsne( dat[,-1] ,perplexity = 10)
    pc.df = data.frame(Shift = p,Scenario = s,Class = dat[,1],pc$x)
    allPC = rbind(allPC,pc.df)
    
  }
}










allPC$Shift = factor(allPC$Shift,labels = c(lr = c("Shift = 2","Shift = 3","Shift = 4")))
allPC$Scenario = factor(allPC$Scenario,labels = c("16S: Covariance shfit",
                                                  "16S: Sparse Mean Shift",
                                                  "WGS: Sparse Mean Shift",
                                                  "WGS: Covariance shfit" ))

tiff(filename = "Figures/dataSimulation_scenarios_expData.tiff",width = 9,height = 6,units = "in",res = 300)
ggplot(allPC,aes(PC1,PC2))+
  geom_point(aes(col = Class,shape = Class),alpha = .5,size = 2)+
  theme_bw()+
  #scale_y_continuous( position="right") + 
  facet_grid(Scenario~Shift,)+
  theme(
    axis.line = element_line(),panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0,face = "bold",size = 8,hjust = 0),
    axis.ticks = element_line(colour = "black",size = 1),
    #axis.title.y = element_blank(),
    axis.title = element_text(size = 8,face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 8),strip.background = element_blank(),
    legend.margin = margin(0,0,0,0,unit = "cm"),
    legend.box.spacing = unit(0,units = "in"), 
    legend.box.margin = margin(0,0,0.1,0,unit = "cm"),
    legend.key.height = unit(.1,units = "in"),
    axis.text =  element_text(size = 8),
    #axis.text.y = element_blank(),
    legend.title =element_text(size = 8),
    plot.caption = element_text(size = 8,face = "italic"))
dev.off()
