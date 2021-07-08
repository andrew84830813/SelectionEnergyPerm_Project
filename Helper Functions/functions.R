
benchMarkPerformance_features = function(tbl,plot_=F, ensemble = c("rf","lda","svmRadial","gbm","pam","nnet"),

                                         cvFolds = 5,cvRepeats = 20){

  #scale data
  tbl[,-1] = scale(tbl[,-1])

  message("1. Perform PCA Reduction")
  pc = prcomp(tbl[,-1])$x


  message("2. Train Model(s)")
  mdls = trainML_Models(trainLRs =  data.frame(pc),
                        testLRs = data.frame(pc),
                        ytrain = tbl[,1], y_test = tbl[,1],
                        cvMethod = "repeatedcv",mtry_ = 1,num_comp = 2,
                        numFolds = cvFolds,
                        numRepeats = cvRepeats,
                        testIDs = NULL,
                        models = ensemble)

  predsAll = mdls$predictionMatrix
  mdlNames = unique(predsAll$model)
  minorityClass = unique(tbl[,1])[1]
  trainPerf = apply(mdls$performance[,1:2], 2, mean)
  trainPerf = data.frame(trainPerf)

  ###------------------------------------------------------------------------------------------------------------------------###
  ## Get Predictions Cross Validation ####
  ###------------------------------------------------------------------------------------------------------------------------###
  cvPreds = data.frame()
  for(jj in 1:length(mdlNames)){
    mdlName =mdlNames[jj]
    mdl = mdls$models[[jj]]
    bestTune = data.frame( mdl$bestTune)
    idConversion = data.frame(rowIndex = 1:nrow(mdl$trainingData),
                              Sample.ID = as.numeric(str_split(rownames(mdl$trainingData),"_",2,T)[,2]))
    probs = left_join(bestTune,mdl$pred)
    probs = left_join(probs,idConversion)
    probs = probs[!is.na(probs$pred),]
    ll = c(which(colnames(probs)==majorityClass),which(colnames(probs)=="pred"),which(colnames(probs)=="rowIndex"))
    probs1 = separate(probs[,-ll],col = 'Resample',into = c("fold","rep"))
    ll = which(colnames(probs1)=="fold")
    probs1 = spread(probs1[,-ll],"rep",minorityClass)
    ll = str_detect(colnames(probs1),"Rep")
    probs1$p = rowMeans(probs1[,ll],na.rm = T)
    folds = unique(probs$Resample)
    dd = data.frame(Model = mdlName,Preds = probs1$p,ID = probs1$Sample.ID )
    cvPreds = rbind(cvPreds,dd)
  }

  cvPreds = cvPreds %>%
    group_by(ID,Model) %>%
    summarise_all(.funs = mean)
  cvPreds1 = spread(cvPreds,key = "Model",value = "Preds")
  preds = rowMeans(cvPreds1[,-1])
  LL = MLmetrics::LogLoss(preds,as.numeric(tbl[,1])-1)
  cvAUC = pROC::auc(tbl[,1],preds)



  message("4. Compute Association Measures")
  ## Energy disco
  d1 = parallelDist::parallelDist(as.matrix(pc))
  enf = energy::disco(x = d1,factors = tbl[,1],distance = T,R = 2)


  ## Combined F
  a.df = data.frame(Type = tbl[,1])
  pmv1 = vegan::adonis(d1~Type,data = a.df,permutations = 1)
  f1 = pmv1$aov.tab$F.Model[1]
  mod = vegan::betadisper((d1),group = tbl[,1])
  bd1 = permutest(mod,permutations = 1)
  f2 = as.numeric(bd1$statistic)


  ns = NA
  g = NA
  netStr = NA

  ##-----------------------------------------*
  ## LR Network ####
  ##-----------------------------------------*
  # feature.df = tbl
  # #### Kruskall Test
  # krus.test_df = data.frame()
  # # tbl =  mst.empirical$features[[1]]
  # lrs_ = data.frame(feature.df[,-1])
  # colnames(lrs_) = colnames(feature.df)[-1]
  # cnames = colnames(lrs_)
  # for(i in 1:ncol(lrs_)){
  #   ph = kruskal.test(x = lrs_[,i],g  = factor(feature.df[,1]))
  #   ph = data.frame(Ratio =cnames[i],pval = ph$p.value,Statistic = ph$statistic )
  #   krus.test_df = rbind(krus.test_df,ph)
  # }
  # #network construction
  # krus.test_df$p.adjust = p.adjust(krus.test_df$pval,method = "BH")
  # pval_level = 0.05
  # fdrLevel =  max(krus.test_df$p.adjust[krus.test_df$pval <= pval_level])
  # krus.test_df$signf = if_else(krus.test_df$p.adjust>0.05,F,T)
  #Importance df
  message("5. Construct log-ratio network")

  krus.test_df = data.frame(Ratio = colnames(tbl[,-1]))
  imp.df = data.frame(krus.test_df,Imp = 1)
  keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
  el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
  g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
  g = igraph::simplify(g, remove.loops = TRUE,
                       edge.attr.comb = igraph_opt("edge.attr.comb"))
  E(g)$weight = if_else((imp.df$Imp)<0,0,(imp.df$Imp))
  cc  = transitivity(g,type = "global")


  if(plot_){

    plot(g,layout = layout_with_fr,vertex.size = log(1/clo(strength(g)))+1,
         vertex.label.cex = .75,edge.curved = .2,edge.width = E(g)$weight*.15,
         edge.arrow.size = .25,edge.arrow.width = 1)
    plot(mod,label = F)
  }

  cn = str_split(colnames(tbl[,-1]),pattern = "___",n = 2,simplify = T)
  cn = n_distinct(c(cn[,1],cn[,2]))


  return(list(
    performance =  data.frame(NumParts = cn,CluseringCoefficient = cc,LogLoss = LL,t(trainPerf),
                              crossValAUC = cvAUC,energyF = enf$statistic,meansF = f1,dispF = f2,combinedF =(f1+f2)),
    graph = g
  )
  )

}




dcvStrength = function (dcv_mat)
{
  Str = NULL
  dcvScores = dcv_mat$dcv
  dcvScores$rowmean[dcvScores$rowmean < 0] = 0
  trainData.md = caret::preProcess(data.frame(dcvScores$rowmean),
                                   method = "range", rangeBounds = c(0, 1))
  scaledScore = stats::predict(trainData.md, data.frame(dcvScores$rowmean))
  el = data.frame(Ratio = dcvScores$Ratio, Score = scaledScore[,
                                                               1])
  el = tidyr::separate(data = el, col = 1, into = c("num",
                                                    "denom"), sep = "___", remove = F)
  g = igraph::graph_from_edgelist(as.matrix(el[, 2:3]))
  igraph::E(g)$weight = el$Score
  nodeStrength = data.frame(Node = names(igraph::strength(g)),
                            Str = igraph::strength(g)) %>% dplyr::arrange(dplyr::desc(Str))
  nodeStrength
}



toGephi = function(Graph,Name){
  el = igraph::as_edgelist(Graph)
  colnames(el) = c("Source","Target")
  el = as.data.frame(el)

  edgeWeight = as.matrix(E(Graph)$weight)
  edgeWeight = as.data.frame(edgeWeight)
  edgeWeight = cbind(el,edgeWeight)
  colnames(edgeWeight)[3] = "Weight"
  write_csv(edgeWeight,paste(Name,"_edgeWeights_",".csv",sep= ""))

  ddf = (Graph)
  ddf = V(ddf)$name
  ddf = as.data.frame(ddf)
  colnames(ddf) = "KO"
  # ex = merge(ddf, KO_Ref_Hierarchy,by = c("KO"))
  ex = cbind.data.frame(ID = ddf[,1],Label = ddf[,1])
  write_csv(ex,paste(Name,"_attributes_",".csv",sep=""))


}


getLogratioFromList = function (Ratio, raMatrix, Class)
{
  Ratio = data.frame(Ratio)
  keyRats = tidyr::separate(Ratio, 1, into = c("Num", "Denom"),
                            sep = "___", remove = F)
  el_ = data.frame(keyRats$Num, keyRats$Denom, keyRats$Ratio)
  ad = selEnergyPermR::stageData(raMatrix, labels = Class)
  selEnergyPermR::getLogRatios(ad$allData, el_)
}
