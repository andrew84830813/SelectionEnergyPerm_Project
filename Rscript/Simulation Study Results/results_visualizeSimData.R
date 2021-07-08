



allPC = data.frame()

for(d in c(50,150,250)){
  for(s in 1:4){
  
  args = c(1,1,250,3)
  seed_ = 1# random seed selection
  sampleDesign = 1; design = if_else(sampleDesign==1,"Balanced","Unbalanced")
  dims = d
  scenario = s
  
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
           ## Scenario-1 = difference in varaince via dirichlet distrubition
           dat = scenario1(seed = seed_,n1 = g1,n2 = g2,dms_ = dims);
           fname = "scenario1";
           path_ = "/nas/longleaf/home/andrew84/selPermEnergy/results2/powerSimData/s1/"
         },
         
         {
           ## Scenario-2 = single feature difference among sparse noisy signal
           dat = scenario2(seed = seed_,n1 = g1,n2 = g2,dms_ = dims);
           fname = "scenario2";
           path_ = "/nas/longleaf/home/andrew84/selPermEnergy/results2/powerSimData/s2/"
         },
         
         {
           ## Scenario-3 = Small Signal and Difference in Sparsity
           dat = scenario3(seed = seed_,n1 = g1,n2 = g2,dms_ = dims);
           fname = "scenario3";
           path_ = "/nas/longleaf/home/andrew84/selPermEnergy/results2/powerSimData/s3/"
         },
         
         {
           ## Scenario-4 = Difference in Mean same covariance additive log normal
           dat = scenario4(seed = seed_,n1 = g1,n2 = g2,dms_ = dims);
           fname = "scenario4";
           path_ = "/nas/longleaf/home/andrew84/selPermEnergy/results2/powerSimData/s4/"
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
  pc.df = data.frame(Dims = d,Scenario = s,Class = dat[,1],pc$x)
  allPC = rbind(allPC,pc.df)

  }
}










allPC$Dims = factor(allPC$Dims,labels = c(lr = c("lr = 1,225","lr = 11,175","lr = 31,125")))
allPC$Scenario = factor(allPC$Scenario,labels = c(lr = c("Scenario-1","Scenario-2","Scenario-3","Scenario-4")))

tiff(filename = "Figures/dataSimulation_scenarios.tiff",width = 9,height = 6,units = "in",res = 300)
ggplot(allPC,aes(PC1,PC2))+
  geom_point(aes(col = Class,shape = Class),alpha = .5,size = 2)+
  theme_bw()+
  #scale_y_continuous( position="right") + 
  facet_grid(Scenario~Dims,)+
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
