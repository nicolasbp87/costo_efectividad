# Functions to conduct sensitivity analysis using linear regression metamodeling
##################################################################################
# Code to convert rmd to R scripts
#' ## script for converting .Rmd files to .R scripts.

#' #### Kevin Keenan 2014
#' http://rstudio-pubs-static.s3.amazonaws.com/12734_0a38887f19a34d92b7311a2c9cb15022.html
#' This function will read a standard R markdown source file and convert it to 
#' an R script to allow the code to be run using the "source" function.
#' 
#' The function is quite simplisting in that it reads a .Rmd file and adds 
#' comments to non-r code sections, while leaving R code without comments
#' so that the interpreter can run the commands.
#' 
#' 
rmd2rscript <- function(infile){
  # read the file
  flIn <- readLines(infile)
  # identify the start of code blocks
  cdStrt <- which(grepl(flIn, pattern = "```{r*", perl = TRUE))
  # identify the end of code blocks
  cdEnd <- sapply(cdStrt, function(x){
    preidx <- which(grepl(flIn[-(1:x)], pattern = "```", perl = TRUE))[1]
    return(preidx + x)
  })
  # define an expansion function
  # strip code block indacators
  flIn[c(cdStrt, cdEnd)] <- ""
  expFun <- function(strt, End){
    strt <- strt+1
    End <- End-1
    return(strt:End)
  }
  idx <- unlist(mapply(FUN = expFun, strt = cdStrt, End = cdEnd, 
                       SIMPLIFY = FALSE))
  # add comments to all lines except code blocks
  comIdx <- 1:length(flIn)
  comIdx <- comIdx[-idx]
  for(i in comIdx){
    flIn[i] <- paste("#' ", flIn[i], sep = "")
  }
  # create an output file
  nm <- strsplit(infile, split = "\\.")[[1]][1]
  flOut <- file(paste(nm, "[rmd2r].R", sep = ""), "w")
  for(i in 1:length(flIn)){
    cat(flIn[i], "\n", file = flOut, sep = "\t")
  }
  close(flOut)
}

#Code to determine whether a particular package is installed and load it, or download it and load it.
#Source: http://stackoverflow.com/questions/5595512/what-is-the-difference-between-require-and-library
if(require("lme4")){
  print("lme4 is loaded correctly")
} else {
  print("trying to install lme4")
  install.packages("lme4")
  if(require(lme4)){
    print("lme4 installed and loaded")
  } else {
    stop("could not install lme4")
  }
}

#Function to determine the number of ticks in ggplot2 plots
number_ticks <- function(n) {function(limits) pretty(limits, n)} #function to determine limits of axis' ticks Source:http://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks-in-ggplot2
number_ticks_x0 <- function(n,change) {function(limits) c(pretty(limits, n),change)} #function to determine limits of axis' ticks Source:http://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks-in-ggplot2
fmt <- function(){function(x) format(x,nsmall = 2,scientific = FALSE)}

#Define offset as a new axis transformation. Source: http://blog.ggplot2.org/post/25938265813/defining-a-new-transformation-for-ggplot2-scales  
offset_trans <- function(offset=0) {
  trans_new(paste0("offset-", format(offset)), function(x) x-offset, function(x) x+offset)
}

OneWaySA<-function(Strategies,Parms,Outcomes,parm,range){
  #Extract parameter column number in Parms matrix
  x<-which(colnames(Parms)==parm)
  dep<-length(Strategies) #Number of dependent variables, i.e., strategies outcomes
  indep<-ncol(Parms) #Number of independent variables, i.e., parameters
  Sim <- data.frame(Outcomes,Parms)
  #Determine range of of the parameer to be plotted
  if (!missing("range")){ #If user defines a range
    vector<-seq(range[1],range[2],length.out=400)
  }
  else{ #Default range given by the domanin of the parameter's sample
    #vector to define 400 samples between the 2.5th and 97.5th percentiles
    y = seq(2.5,97.5,length=400) 
    j = round(y*(length(Parms[,x])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
    vector<-sort(Parms[j,x])    
  }
  
  #Generate a formula by pasting column names for both dependent and independent variables. Imposes a 1 level interaction
  f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep],collapse=','), ') ~ (','poly(',parm,',2,raw=TRUE)+' ,paste(colnames(Parms)[-x], collapse='+'),')'))
  #f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep],collapse=','), ') ~ (', paste(colnames(Sim)[(dep+1):(dep+indep)], collapse='+'),')^2'))
  #f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep],collapse=','), ') ~', paste(colnames(Sim)[(dep+1):(dep+indep)], collapse='+')))
  #Run Multiple Multivariate Regression (MMR) Metamodel
  Oway.mlm = lm(f,data=Sim)

  #Generate matrix to use for prediction 
  Sim.fit<-matrix(rep(colMeans(Parms)),nrow=length(vector),ncol=ncol(Parms), byrow=T)
  Sim.fit[,x]<-vector
  Sim.fit<-data.frame(Sim.fit) #Transform to data frame, the format required for predict
  colnames(Sim.fit)<-colnames(Parms) #Name data frame's columns with parameters' names
  
  #Predict Outcomes using MMMR Metamodel fit
  plotdata = data.frame(predict(Oway.mlm, newdata = Sim.fit))
  colnames(plotdata) <- Strategies #Name the predicted outcomes columns with strategies names
  
#   #Determine intersection points  
#   Optimal <- max.col(plotdata)
#   Opt.factors <- unique(factor(Optimal))
#   Opt.lev<-rev(destring(levels(Opt.factors)))
#   sections<-nlevels(OptStr)-1
#   change<-rep(0,sections)
#   for (i in 1:sections){
#     ind<-which(plotdata[,Opt.lev[i+1]]<plotdata[,Opt.lev[i]])
#     change[i] <- round(vector[ind[i]-1],2)
#   }
#   out<-plotdata[change,]
  
  #Reshape dataframe for ggplot
  plotdata = stack(plotdata, select=Strategies) #
  plotdata = cbind(Sim.fit, plotdata) #Append parameter's dataframe to predicted outcomes dataframe
  
  #A simple trick I use to define my variables in my functions environment
  #Borrowed from http://stackoverflow.com/questions/5106782/use-of-ggplot-within-another-function-in-r
  plotdata$parm<-plotdata[,parm];
  
  txtsize<-12 #Text size for the graphs
  ggplot(data = plotdata, aes(x = parm, y = values, color = ind)) +
    geom_line() +
    ggtitle("One-way sensitivity analysis \n Net Health Benefit") + 
    xlab(parm) +
    ylab("E[NHB]") +
    scale_colour_hue("Strategy", l=50) +
    scale_x_continuous(breaks=number_ticks(6)) + #Adjust for number of ticks in x axis
    scale_y_continuous(breaks=number_ticks(6)) +
#     geom_segment(aes(x = change, y = 0, xend = change, yend = healthy[change]), 
#                  arrow = arrow(length = unit(0.2, "cm"),type = "closed"), 
#                  colour="gray54",
#                  linetype="dotted") +
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

TwoWaySA<-function(Strategies,Parms,Outcomes,parm1,parm2,range1,range2){
  #Extract parameter column number in Parms matrix
  x1<-which(colnames(Parms)==parm1)
  x2<-which(colnames(Parms)==parm2)
  dep<-length(Strategies) #Number of dependent variables, i.e., strategies
  indep<-ncol(Parms) #Number of independent variables, i.e., parameters
  
  Sim <- data.frame(Outcomes,Parms)
  
  #Determine range of of the parameer to be plotted
  if (!missing("range1")&!missing("range2")){ #If user defines a range
    vector1<-seq(from=range1[1],to=range1[2],length.out=301)
    vector2<-seq(from=range2[1],to=range2[2],length.out=301)
  }
  else if (!missing("range1")&missing("range2")){ #Default range given by the domanin of the parameter's sample
    #vector to define 400 samples between the 2.5th and 97.5th percentiles
    vector1<-seq(from=range1[1],to=range1[2],length.out=301)
    y2 = seq(2.5,97.5,length.out=301)
    j2 = round(y2*(length(Parms[,x2])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
    vector2<-sort(Parms[j2,x2])    
  }
  else if (missing("range1")&!missing("range2")){ #Default range given by the domanin of the parameter's sample
    #vector to define 400 samples between the 2.5th and 97.5th percentiles
    vector2<-seq(from=range2[1],to=range2[2],length.out=301)
    y1 = seq(2.5,97.5,length.out=301)
    j1 = round(y1*(length(Parms[,x1])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
    vector1<-sort(Parms[j1,x1])    
  }
  else{
    y1 = seq(2.5,97.5,length.out=301) 
    y2 = seq(2.5,97.5,length.out=301) 
    j1 = round(y1*(length(Parms[,x1])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
    j2 = round(y2*(length(Parms[,x2])/100))
    vector1<-sort(Parms[j1,x1])
    vector2<-sort(Parms[j2,x2])
  }
  #Generate a formula by pasting column names for both dependent and independent variables
  f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep],collapse=','), ') ~ (','poly(',parm1,',8)+','poly(',parm2,',8)+' ,paste(colnames(Parms)[c(-x1,-x2)], collapse='+'),')'))
  #f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep],collapse=','), ') ~ (', paste(colnames(Sim)[(dep+1):(dep+indep)], collapse='+'),')^2'))
  #f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep],collapse=','), ') ~', paste(colnames(Sim)[(dep+1):(dep+indep)], collapse='+')))
  #Run Multiple Multivariate Regression (MMR) Metamodel
  Tway.mlm = lm(f,data=Sim)
  
  #vector to define 400 samples between the 2.5th and 97.5th percentiles !!!HAS TO BE EQUALLY SPACED!!!
#   y1 = seq(2.5,97.5,length=300) 
#   y2 = seq(2.5,97.5,length=300) 
#   j1 = round(y1*(length(Parms[,x1])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
#   j2 = round(y2*(length(Parms[,x2])/100))
#   vector1<-sort(Parms[j1,x1])
#   vector2<-sort(Parms[j2,x2])
#   vector1<-seq(from=range1[1],to=range1[2],length.out=301)
#   vector2<-seq(from=range2[1],to=range2[2],length.out=301)
  TWSA <- expand.grid(parm1=vector1,parm2=vector2)
  
  #Generate matrix to use for prediction 
  Sim.fit<-matrix(rep(colMeans(Parms)),nrow=nrow(TWSA),ncol=ncol(Parms), byrow=T)
  Sim.fit[,x1]<-TWSA[,1]
  Sim.fit[,x2]<-TWSA[,2]
  Sim.fit<-data.frame(Sim.fit) #Transform to data frame, the format required for predict
  colnames(Sim.fit)<-colnames(Parms) #Name data frame's columns with parameters' names
  
  #Predict Outcomes using MMMR Metamodel fit
  Sim.TW = data.frame(predict(Tway.mlm, newdata = Sim.fit))
  #Find optimal strategy in terms of maximum Outcome
  Optimal <- max.col(Sim.TW)
  #Get Outcome of Optimal strategy
  OptimalOut<-apply(Sim.TW,1,max)
  
  #colnames(plotdata) <- Strategies #Name the predicted outcomes columns with strategies names
  #plotdata = stack(plotdata, select=Strategies) #
  plotdata = Sim.fit #Append parameter's dataframe to predicted outcomes dataframe
  
  #A simple trick I use to define my variables in my functions environment
  #Borrowed from http://stackoverflow.com/questions/5106782/use-of-ggplot-within-another-function-in-r
  plotdata$parm1<-plotdata[,parm1];
  plotdata$parm2<-plotdata[,parm2];
  
  plotdata$Strategy<-factor(Optimal,labels=Strategies)
  plotdata$value<-OptimalOut
  
  txtsize<-12
  ggplot(plotdata, aes(x=parm1,y=parm2))+ 
    geom_tile(aes(fill=Strategy)) +
    theme_bw() +
    ggtitle(expression(atop("Two-way sensitivity analysis", 
                            atop("Net Health Benefit")))) + 
    scale_fill_discrete("Strategy: ", l=50)+
    xlab(parm1)+
    ylab(parm2)+
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

TornadoOpt <-function(Parms,Outcomes){
  # Grouped Bar Plot
  # Determine the overall optimal strategy
  opt<-which.max(colMeans(Outcomes))
  # calculate min and max vectors of the parameters (e.g., lower 2.5% and 97.5%)
  X <- as.matrix(Parms)
  y <- as.matrix(Outcomes[,opt])
  ymean <- mean(y)
  n <- nrow(Parms)
  nParams <- ncol(Parms)
  #paramNames <- Names[seq(8,7+nParams)]
  paramNames <- colnames(Parms)
  Parms.sorted <- apply(Parms,2,sort,decreasing=F)#Sort in increasing order each column of Parms
  lb <- 2.5
  ub <- 97.5 
  Xmean <- rep(1,nParams) %*% t(colMeans(X))
  XMin <- Xmean
  XMax <- Xmean
  paramMin <- as.vector(Parms.sorted[round(lb*n/100),])
  paramMax <- as.vector(Parms.sorted[round(ub*n/100),])
  paramNames2 <- paste(paramNames, "[", round(paramMin,2), ",", round(paramMax,2), "]")
  
  diag(XMin) <- paramMin
  diag(XMax) <- paramMax
  
  XMin <- cbind(1, XMin)
  XMax <- cbind(1, XMax)
  
  X <- cbind(1,X)
  B <- solve(t(X) %*% X) %*% t(X) %*% y
  yMin <- XMin %*% B - ymean
  yMax <- XMax %*% B - ymean
  ySize <- abs(yMax - yMin) 
  
  rankY<- order(ySize)
  xmin <- min(c(yMin, yMax)) + ymean
  xmax <- max(c(yMin, yMax)) + ymean
  
  Tor <- data.frame(
    Parameter=c(paramNames2[rankY],paramNames2[rankY]),  
    Level=c(rep("Low",nParams),rep("High",nParams)),
    value=ymean+c(yMin[rankY],yMax[rankY]),
    sort=seq(1,nParams)
  )
  #re-order the levels in the order of appearance in the data.frame
  Tor$Parameter2 <- factor(Tor$Parameter, as.character(Tor$Parameter))
  #Define offset as a new axis transformation. Source: http://blog.ggplot2.org/post/25938265813/defining-a-new-transformation-for-ggplot2-scales  
  offset_trans <- function(offset=0) {
    trans_new(paste0("offset-", format(offset)), function(x) x-offset, function(x) x+offset)
  }
  #Plot the Tornado diagram. Source: http://comments.gmane.org/gmane.comp.lang.r.ggplot2/8531
  txtsize<-12
  ggplot(Tor[Tor$Level=="Low",], aes(x=Parameter2,y=value, fill=Level)) +
    geom_bar(stat="identity") +
    ggtitle("Tornado Diagram")+
    scale_fill_discrete("Parameter Level: ", l=50)+
    scale_y_continuous(name="Net Benefit",trans=offset_trans(offset=ymean)) +
    scale_x_discrete(name="Parameter") +
    geom_bar(data=Tor[Tor$Level=="High",], aes(x=Parameter2,y=value, fill=Level), stat="identity") +
    geom_hline(yintercept = ymean, linetype = "dotted", size=0.5) +
    #annotate("text", x = 0, y = ymean, label = "Mean", size=5) +
    theme_bw()+
#     theme(axis.text.x = element_text(angle = 0, hjust = 1),
#           axis.ticks.y = element_blank()) +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize,angle = 0, hjust = 1),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize),
          axis.ticks.y = element_blank())+
    coord_flip()  
}

TornadoAll <-function(Strategies,Parms,Outcomes){
  #Outcomes=NHB
  opt<-which.max(colMeans(Outcomes))
  # calculate min and max vectors of the parameters (e.g., lower 2.5% and 97.5%)
  X <- as.matrix(Parms)
  y <- as.matrix(Outcomes[,opt])
  Y <- as.matrix(Outcomes)
  ymean <- mean(y)
  n <- nrow(Parms)
  nParams <- ncol(Parms)
  #paramNames <- Names[seq(8,7+nParams)]
  paramNames <- colnames(Parms)
  Parms.sorted <- apply(Parms,2,sort,decreasing=F)#Sort in increasing order each column of Parms
  lb <- 2.5
  ub <- 97.5 
  Xmean <- rep(1,nParams) %*% t(colMeans(X))
  XMin <- Xmean
  XMax <- Xmean
  paramMin <- as.vector(Parms.sorted[round(lb*n/100),])
  paramMax <- as.vector(Parms.sorted[round(ub*n/100),])
  
  diag(XMin) <- paramMin
  diag(XMax) <- paramMax
  
  XMin <- cbind(1, XMin)
  XMax <- cbind(1, XMax)
  
  X <- cbind(1,X)
  B <- solve(t(X) %*% X) %*% t(X) %*% y
  #install.packages("matrixStats")
  library(matrixStats)
  bigBeta <- solve(t(X) %*% X) %*% t(X) %*% Y
  yMin <- rowMaxs(XMin %*% bigBeta - ymean)
  yMax <- rowMaxs(XMax %*% bigBeta - ymean)
  ySize <- abs(yMax - yMin) 
  
  rankY<- order(ySize)
  xmin <- min(c(yMin, yMax)) + ymean
  xmax <- max(c(yMin, yMax)) + ymean
  
  paramNames2 <- paste(paramNames, "[", round(paramMin,2), ",", round(paramMax,2), "]")
  
  strategyNames<-Strategies
  strategyColors <- c("red","darkgreen","blue")
  
  ## Polygon graphs:
  nRect <- 0
  x1Rect <- NULL
  x2Rect <- NULL
  ylevel <- NULL
  colRect <- NULL
  
  for (p in 1:nParams){
    xMean <- colMeans(X)
    xStart = paramMin[rankY[p]]
    xEnd = paramMax[rankY[p]]
    xStep = (xEnd-xStart)/1000
    for (x in seq(xStart,xEnd, by = xStep)){
      #for each point determine which one is the optimal strategy
      xMean[rankY[p] + 1] <- x 
      yOutcomes <- xMean %*% bigBeta
      yOptOutcomes <- max(yOutcomes)
      yOpt <- which.max(yOutcomes)
      if (x == xStart){
        yOptOld <- yOpt
        y1 <- yOptOutcomes
      }
      #if yOpt changes, then plot a rectangle for that region
      if (yOpt != yOptOld | x == xEnd){
        nRect <- nRect + 1
        x1Rect[nRect] <- y1
        x2Rect[nRect] <- yOptOutcomes
        ylevel[nRect] <- p
        colRect[nRect] <- strategyColors[yOptOld]
        yOptOld <- yOpt
        y1 <- yOptOutcomes
      }
    }
  }
  
  txtsize <- 12
  d=data.frame(x1=x2Rect, x2=x1Rect, y1=ylevel-0.4, y2=ylevel+0.4, t=colRect, r = ylevel)
  ggplot(d, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = t)) +
    ggtitle("Torando Diagram") + 
    xlab("Expected NHB") +
    ylab("Parameters") + 
    geom_rect()+
    theme_bw() + 
    scale_y_continuous(limits = c(0.5, nParams + 0.5),breaks=seq(1:8), labels=paramNames2[rankY]) +
    #scale_y_discrete(breaks=seq(1:8), labels=paramNames2[rankY]) + 
    scale_fill_discrete(name="Optimal\nStrategy",
                        #breaks=c("ctrl", "trt1", "trt2"),
                        labels=strategyNames,
                        l=50) + 
    geom_vline(xintercept=ymean, linetype="dotted") + 
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

PlaneCE<-function(Strategies,Outcomes){
  ndep<-length(Strategies)*2 #Determine number of outcomes for all starteges, i.e., cost and effectiveness
  ind_c<-seq(1,(ndep-1),by=2) #Index to extract the costs from matrix Outcomes
  ind_e<-seq(2,ndep,by=2) #Index to extract the effectiveness from matrix Outcomes
  Cost<-melt(Outcomes[,ind_c],variable_name="Strategy")
  levels(Cost$Strategy)<-Strategies
  Eff<-melt(Outcomes[,ind_e],variable_name="Strategy")
  levels(Eff$Strategy)<-Strategies
  CE<-cbind(Cost,Eff[,2])
  colnames(CE)<-c("Strategy","Cost","Effectiveness")
  
  #Dataframe with means of strategies. Source: http://stackoverflow.com/questions/18729724/ggplot2-scatter-plot-with-overlay-of-means-and-bidirectional-sd-bars
  Means <- ddply(CE,.(Strategy),summarise,
                 N = length(Cost),
                 Cost.mean = mean(Cost),
                 Eff.mean = mean(Effectiveness))
  
  #Define ggplot object
  txtsize<-12
  ggplot(Means, aes(x = Eff.mean, y = Cost.mean, color=Strategy)) + 
    geom_point(size=4, aes(shape=Strategy)) +
    ggtitle("Cost-Effectiveness Plane") +
    scale_colour_discrete(l=50) +  # Use a slightly darker palette than normal
    scale_y_continuous(labels = dollar)+
    scale_x_continuous(breaks=number_ticks(6), labels=comma)+
    xlab("Effectiveness")+
    ylab("Cost")+
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

ScatterCE<-function(Strategies,Outcomes){
  ndep<-length(Strategies)*2 #Determine number of outcomes for all starteges, i.e., cost and effectiveness
  ind_c<-seq(1,(ndep-1),by=2) #Index to extract the costs from matrix Outcomes
  ind_e<-seq(2,ndep,by=2) #Index to extract the effectiveness from matrix Outcomes
  Cost<-melt(Outcomes[,ind_c],variable_name="Strategy")
  levels(Cost$Strategy)<-Strategies
  Eff<-melt(Outcomes[,ind_e],variable_name="Strategy")
  levels(Eff$Strategy)<-Strategies
  CE<-cbind(Cost,Eff[,2])
    colnames(CE)<-c("Strategy","Cost","Effectiveness")

  # Ellipses code
  df_ell <- data.frame() #create an empty dataframe
  # for each level in df$groups 
  for(g in levels(CE$Strategy)){
    # create 100 points per variable around the mean of each group
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(CE[CE$Strategy==g,], 
                                                     ellipse(cor(Effectiveness, Cost), 
                                                     scale=c(sd(Effectiveness),sd(Cost)), 
                                                     centre=c(mean(Effectiveness),mean(Cost)))
    )),group=g))
  }

  #Dataframe with means of strategies. Source: http://stackoverflow.com/questions/18729724/ggplot2-scatter-plot-with-overlay-of-means-and-bidirectional-sd-bars
  Means <- ddply(CE,.(Strategy),summarise,
                 N = length(Cost),
                 Cost.mean = mean(Cost),
                 Eff.mean = mean(Effectiveness))
  
  #Define ggplot object
  txtsize<-12
  ggplot(CE, aes(x=Effectiveness, y=Cost, color=Strategy)) + 
    geom_point(size=0.7) +
    geom_point(data = Means,aes(x = Eff.mean, y = Cost.mean, shape=Strategy),size=8,fill="white") +
    #guides(shape=FALSE)+ #Source: http://stackoverflow.com/questions/14604435/turning-off-some-legends-in-a-ggplot
    geom_text(data = Means,aes(x = Eff.mean, y = Cost.mean, label=c(1,2,3)),size=5,colour="gray",alpha=1) +
    geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=1, linetype=2, alpha=1) + # draw ellipse lines
    #geom_abline(intercept=, slope = 10000) +
    #geom_line(data = Means,aes(x = Eff.mean, y = Cost.mean),colour="black")+
    ggtitle("Cost-Effectiveness Scatterplot") +
    scale_colour_discrete(l=50) +  # Use a slightly darker palette than normal
    scale_y_continuous(labels = dollar)+
    scale_x_continuous(breaks=number_ticks(6), labels=comma)+
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

ScatterICER<-function(Strategies,Outcomes){
  ndep<-length(Strategies)*2 #Determine number of outcomes for all starteges, i.e., cost and effectiveness
  ind_c<-seq(1,(ndep-1),by=2) #Index to extract the costs from matrix Outcomes
  ind_e<-seq(2,ndep,by=2) #Index to extract the effectiveness from matrix Outcomes
  Cost<-melt(Outcomes[,ind_c],variable_name="Strategy")
  levels(Cost$Strategy)<-Strategies
  Eff<-melt(Outcomes[,ind_e],variable_name="Strategy")
  levels(Eff$Strategy)<-Strategies
  CE<-cbind(Cost,Eff[,2])
  colnames(CE)<-c("Strategy","Cost","Effectiveness")
  
  # Ellipses code
  df_ell <- data.frame() #create an empty dataframe
  # for each level in df$groups 
  for(g in levels(CE$Strategy)){
    # create 100 points per variable around the mean of each group
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(CE[CE$Strategy==g,], 
                                                     ellipse(cor(Effectiveness, Cost), 
                                                             scale=c(sd(Effectiveness),sd(Cost)), 
                                                             centre=c(mean(Effectiveness),mean(Cost)))
    )),group=g))
  }
  
  #Dataframe with means of strategies. Source: http://stackoverflow.com/questions/18729724/ggplot2-scatter-plot-with-overlay-of-means-and-bidirectional-sd-bars
  Means <- ddply(CE,.(Strategy),summarise,
                 N = length(Cost),
                 Cost.mean = mean(Cost),
                 Eff.mean = mean(Effectiveness))
  
  #Define ggplot object
  txtsize<-12
  ggplot(CE, aes(x=Effectiveness, y=Cost, color=Strategy)) + 
    geom_point(size=0.7) +
    geom_point(data = Means,aes(x = Eff.mean, y = Cost.mean, shape=Strategy),size=8,fill="white") +
    #guides(shape=FALSE)+ #Source: http://stackoverflow.com/questions/14604435/turning-off-some-legends-in-a-ggplot
    geom_text(data = Means,aes(x = Eff.mean, y = Cost.mean, label=c(1,2)),size=5,colour="gray",alpha=1) +
    geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=1, linetype=2, alpha=1) + # draw ellipse lines
    #geom_abline(intercept=, slope = 10000) +
    #geom_line(data = Means,aes(x = Eff.mean, y = Cost.mean),colour="black")+
    ggtitle("Incremental Cost-Effectiveness Scatterplot") +
    xlab("Incremental Effectiveness")+
    ylab("Incremental Cost")+
    scale_colour_discrete(l=50) +  # Use a slightly darker palette than normal
    scale_y_continuous(labels = dollar)+
    scale_x_continuous(breaks=number_ticks(6), labels=comma)+
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

ScatterICER2<-function(Strategies,Outcomes){
  ind_c<-c(1,3) #Index to extract the costs from matrix Outcomes
  ind_e<-c(2,4) #Index to extract the effectiveness from matrix Outcomes
  Cost<-Outcomes[,3]-Outcomes[,1]
  Eff<-Outcomes[,4]-Outcomes[,2]
  CE<-data.frame(Strategy=paste(Strategies[2],"vs",Strategies[1],"(reference)"),Cost,Eff)
  colnames(CE)<-c("Strategy","Cost","Effectiveness")
  
  # Ellipses code
  df_ell <- data.frame() #create an empty dataframe
  # for each level in df$groups 
  for(g in levels(CE$Strategy)){
    # create 100 points per variable around the mean of each group
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(CE[CE$Strategy==g,], 
                                                     ellipse(cor(Effectiveness, Cost), 
                                                             scale=c(sd(Effectiveness),sd(Cost)), 
                                                             centre=c(mean(Effectiveness),mean(Cost)))
    )),group=g))
  }
  
  #Dataframe with means of strategies. Source: http://stackoverflow.com/questions/18729724/ggplot2-scatter-plot-with-overlay-of-means-and-bidirectional-sd-bars
  Means <- ddply(CE,.(Strategy),summarise,
                 N = length(Cost),
                 Cost.mean = mean(Cost),
                 Eff.mean = mean(Effectiveness))
  
  #Define ggplot object
  txtsize<-12
  ggplot(CE, aes(x=Effectiveness, y=Cost, color=Strategy)) + 
    geom_point(size=0.7) +
    geom_point(data = Means,aes(x = Eff.mean, y = Cost.mean, shape=Strategy),size=8,fill="white") +
    #geom_text(data = Means,aes(x = Eff.mean, y = Cost.mean, label=c(1)),size=5,colour="black",alpha=1) +
    geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=0.5, linetype=2, alpha=1) + # draw ellipse lines
    #geom_line(data = Means,aes(x = Eff.mean, y = Cost.mean),colour="black")+
    xlab("Incremental Effectiveness (# of preterm births averted)")+
    ylab("Incremental Cost (million $)")+
    ggtitle("Incremental Cost-Effectiveness Scatterplot") +
    scale_colour_discrete(l=50) +  # Use a slightly darker palette than normal
    scale_y_continuous(labels = dollar)+
    scale_x_continuous(breaks=number_ticks(6),labels = comma)+
    stat_abline(intercept=0, slope=0, linetype="dotted") + # add a reference line
    stat_vline(xintercept=0, slope=0, linetype="dotted") + # add a reference line
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

CEAC<-function(range,Strategies,Outcomes){
  # Outcomes must be ordered in a way that for each strategy the cost must appear first then the effectiveness
  lambda<-seq(range[1],range[2], length.out=range[3])
  NHB <- array(0, dim=c(dim(Outcomes)[1],length(Strategies))) # Matrix to store NHB for each strategy
    colnames(NHB)<-Strategies
  CEA<-array(0,dim=c(length(lambda),length(Strategies)))
  costInd<-seq(1,2*length(Strategies),by=2) #vector to index costs
  effInd<-seq(2,2*length(Strategies),by=2) #vector to index effectiveness
  
  for(l in 1:length(lambda)){
    NHB <-  Outcomes[,effInd]-Outcomes[,costInd]/lambda[l] # Effectiveness minus Costs, with vector indexing
    Max.NHB <- max.col(NHB)
    Optimal <- array(0, dim=c(dim(Outcomes)[1],length(Strategies))) # Matrix with dummy variables indicating if strategy is optimal
    for (j in 1:dim(Optimal)[1]){
      k<-Max.NHB[j]
      Optimal[j,k]<-1
    }
    CEA[l,]<-colMeans(Optimal)
  }
  CEA<-data.frame(cbind(lambda,CEA))
    colnames(CEA)<-c("Lambda", Strategies)
  #matplot(CEA[,1],CEA[,2:4])
  
  CEAC<-melt(CEA, id.vars = "Lambda") 
  
  txtsize<-12
  ggplot(data = CEAC, aes(x = Lambda, y = value, color = variable)) +
    geom_point() +
    geom_line() +
    ggtitle("Cost-Effectiveness Acceptability Curves") + 
    scale_colour_hue("Strategies: ",l=50) +
    scale_x_continuous(breaks=number_ticks(6))+
    xlab(expression("Willingness to Pay "*lambda*" ($/QALY)")) +
    ylab("Pr Cost-Effective") +
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
         legend.key = element_rect(colour = "black"),
         legend.text = element_text(size = txtsize),
         title = element_text(face="bold", size=15),
         axis.title.x = element_text(face="bold", size=txtsize),
         axis.title.y = element_text(face="bold", size=txtsize),
         axis.text.y = element_text(size=txtsize),
         axis.text.x = element_text(size=txtsize))
}

EVPI<-function(range,Strategies,Outcomes){
  # Outcomes must be ordered in a way that for each strategy the cost must appear first then the effectiveness
  lambda<-seq(range[1],range[2], length.out=range[3])
  NHB <- array(0, dim=c(dim(Outcomes)[1],length(Strategies))) # Matrix to store NHB for each strategy
  colnames(NHB)<-Strategies
  EVPI<-rep(0,length(lambda))
  costInd<-seq(1,2*length(Strategies),by=2) #vector to index costs
  effInd<-seq(2,2*length(Strategies),by=2) #vector to index effectiveness
  
  for(l in 1:length(lambda)){
    NHB <-  lambda[l]*Outcomes[,effInd]-Outcomes[,costInd] # Effectiveness minus Costs, with vector indexing
    maxEV<-which.max(colMeans(NHB)) #Compute expected NHB for all strategies and determine which one has the highest   
    VPI <- rep(0, nrow(Outcomes)) # Matrix with difference between optimal strategy and strategy with max NHB per simulation iteration
    maxNHB<-max.col(NHB)
    #VPI<-NHB[,max.col(NHB)]-NHB[,maxEV]
    #VPI<-apply(NHB,1,function(NHB,maxEV) NHB[,max.col(NHB)]-NHB[,maxEV])
    for (j in 1:length(VPI)){
      k<-which.max(NHB[j,])
      VPI[j]<-NHB[j,k]-NHB[j,maxEV]
    }
    EVPI[l]<-mean(VPI)
  }
  EVPI<-data.frame(cbind(lambda,EVPI))
  colnames(EVPI)<-c("Lambda", "EVPI")
  #matplot(CEA[,1],CEA[,2:4])
  
  txtsize<-12
  ggplot(data = EVPI, aes(x = Lambda, y = EVPI)) +
    geom_point() +
    geom_line() +  
    ggtitle("Expected Value of Perfect Information") + 
    scale_x_continuous(breaks=number_ticks(6))+
    scale_y_continuous(breaks=number_ticks(6))+
    xlab(expression("Willingness to Pay "*lambda*" ($/QALY)")) +
    ylab("EVPI ($)") +
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

MNL.SA <- function(parm,range,Strategies,Parms,Outcomes){
  #Extract parameter column number in Parms matrix
  x<-which(colnames(Parms)==parm)
  #Calculate the preferred strategy (i.e., optimal startegy) in terms of its cost-efectiveness 
  #for each of the simulations.
  Optimal <- data.frame(max.col(Outcomes)); names(Optimal)<-"Strategy"
  
  MNL <- data.frame(Optimal,Parms)
  MNL.vglm = vglm(Strategy ~ ., data = MNL, family=multinomial)
  
  y = seq(2.5,97.5,length=400) #vector to define 400 samples between the 2.5th and 97.5th percentiles
  j = round(y*(length(Parms[,x])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
  vector<-sort(Parms[j,x])
  #vector<-seq(from=range[1],to=range[2],length.out=201)
  
  #Generate matrix to use for prediction 
  Sim.fit<-matrix(rep(colMeans(Parms)),nrow=length(vector),ncol=ncol(Parms), byrow=T)
  Sim.fit[,x]<-vector
  MNL.fit<-data.frame(Sim.fit) #Transform to data frame, the format required for predict
  colnames(MNL.fit)<-colnames(Parms) #Name data frame's columns with parameters' names
  
  #Predict Outcomes using MMMR Metamodel fit
  plotdata = data.frame(predict(MNL.vglm, newdata = MNL.fit, type = "response"))
  
  colnames(plotdata) <- Strategies
  plotdata = stack(plotdata, select=Strategies)
  plotdata = cbind(MNL.fit, plotdata) 
  
  plotdata$parm<-plotdata[,parm];
  
  txtsize<-12 #Text size for the graphs
  ggplot(data = plotdata, aes(x = parm, y = values, color = ind)) +
    geom_point(size = 2) + #maybe shape=3;18;21;124; Other shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
    geom_line() +
    ggtitle("Multinomial sensitivity analysis") + 
    xlab(parm) +
    ylab("Probability of Strategy Being Optimal") +
    scale_colour_hue("Strategy", l=50) +
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
 
#   MNL.fit <- expand.grid(parm=vector)
#   Sim.MNL <- array(0,dim=c(dim(MNL.fit)[1],length(Parms)))
#   for (i in 1:(length(Parms))){
#     if (i==x){
#       Sim.MNL[,i] <-MNL.fit[,1]
#     } else {
#       Sim.MNL[,i] <-mean(Parms[,i])
#     }
#   }
#   
#   MNL.fit = data.frame(Sim.MNL)
#   colnames(MNL.fit)<-colnames(Parms)
#   plotdata = data.frame(predict(MNL.vglm, newdata = MNL.fit, type = "response"))
#   colnames(plotdata) <- Strategies
#   plotdata = stack(plotdata, select=Strategies)
#   plotdata = cbind(MNL.fit, plotdata) 
#   
#   txtsize<-12 #Text size for the graphs
#   ggplot(data = plotdata, aes(x = muDieCancer, y = values, color = ind)) +
#     geom_point(size = 2) + #maybe shape=3;18;21;124; Other shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
#     geom_line() +
#     ggtitle("Multinomial sensitivity analysis") + 
#     xlab("muDieCancer") +
#     ylab("Probability of Strategy Being Optimal") +
#     scale_colour_hue("Strategy", l=50) +
#     theme_bw() +
#     theme(legend.position="bottom",legend.title=element_text(size = txtsize),
#           legend.key = element_rect(colour = "black"),
#           legend.text = element_text(size = txtsize),
#           title = element_text(face="bold", size=15),
#           axis.title.x = element_text(face="bold", size=txtsize),
#           axis.title.y = element_text(face="bold", size=txtsize),
#           axis.text.y = element_text(size=txtsize),
#           axis.text.x = element_text(size=txtsize))
}