################################################################################
#
# Code to accompany: "The batched stepped wedge design: a design robust to 
#                     delays in cluster recruitment" 
#                     by J Kasza, R Bowden, R Hooper, AB Forbes.
#
# Questions or comments to J Kasza, jessica.kasza@monash.edu
#
#
# BatchedSW File 1: functions required to run the simulation study
#
# 2022-04-11
################################################################################

library("lme4")
library("geepack")
library("glmm")
library("swdpwr")


#NOTE that the Binary data generation requires downloading the binGEN function from
#https://github.com/lifanfrank/Li_SIM_SuppData_R_Code.git
#(see reference therein for accompanying paper by Fan Li et al.)

#SW schematic
SWdesmat <- function(Ts) {
  Xsw <- matrix(data=0, ncol = (Ts), nrow = (Ts-1))
  for(i in 1:(Ts-1)) {
    Xsw[i,(i+1):Ts] <- 1
  }
  return(Xsw)
}
#Batched SW schematic
batchSWscheme <- function(olap, S, Ts,  Kseq){
  #Creates a batched SW design schematic with
  #S batches of Ts-period stepped wedge designs, 
  #with overlap between batches of olap periods
  #and K clusters in each sequence
  
  g<- Ts- olap #gap between each batch and the next
  
  SWbatch <- matrix(data=NA, nrow=S*(Ts-1), ncol=Ts + (S-1)*(Ts-olap))  
  
  for(i in 1:S){
    SWbatch[((i-1)*(Ts-1) + 1):(i*(Ts-1)) , 
            ((i-1)*g +1):((i-1)*g + Ts)] <- SWdesmat(Ts)
    
  }
  
  return( SWbatch[sort(rep(1:nrow(SWbatch), Kseq)), ])
  
  
}

##################################################################################
#Functions for continuous outcomes
BatchedPowerSim_conts_wrap <- function(nrep, Ts, nbatch, olap, Kseq, m, TimeEffsInd, Rxeff, ICC, CAC){
  #Wrapper function forthe BatchedPowerSim function. 
  #nrep = number of replications in the simulation study
  #Ts = number of periods in each batch of the design
  #olap = number of periods of overlap between successive batches
  #Kseq = number of clusters assigned to each sequence
  #m = number of observations in each cluster in each period
  #TimeEffsInd = 0 if time effects are shared across batches; 
  #            = 1 if separate time effects for each batch
  #Rxeff = the treatment effect
  #ICC = intracluster correlation
  #CAC = cluster autocorrelation
  
  #Generate the design matrix
  DesMatrix <- batchSWscheme(olap, nbatch, Ts, 1)
  Batch <- rep(seq(1:nbatch), each = (Ts-1))
  
  
  if(TimeEffsInd == 0){
    #Non-null time effects, shared across batches
    TimeEffs <- 0.1*matrix(data=seq(1:ncol(DesMatrix)), nrow=nrow(DesMatrix), 
                           ncol=ncol(DesMatrix), byrow=TRUE)
    TimeEffs <- TimeEffs*(!is.na(DesMatrix))
  }
  else if(TimeEffsInd == 1){
    #Non-null time effects, NOT shared across batches
    TimeEffs <- matrix(data=rnorm(n=ncol(DesMatrix)*nbatch), nrow=nbatch, ncol=ncol(DesMatrix) )
    TimeEffs <- TimeEffs[sort(rep(1:nrow(TimeEffs), (Ts-1))), ]
    TimeEffs <- TimeEffs*(!is.na(DesMatrix))   
  }
  
  output <- replicate(nrep, BatchedPowerSim_conts(DesMatrix, Kseq, m, Batch, TimeEffs, Rxeff, ICC, CAC))
  
  #3 models have been fit to the simulated datasets:
  #1. Model with time by batch effects
  #2. Model with shared time effects
  #3. Model with time and batch effects
  #For each, want the following:
  #- bias
  #- empirical standard error
  #- MSE
  #- Average Model SE
  #- Rejection percentage

  myresults <- NULL
  #Bias
  myresults[1] <- mean(output[1,]-Rxeff)/nrep
  myresults[2] <- mean(output[3,]-Rxeff)/nrep
  myresults[3] <- mean(output[5,]-Rxeff)/nrep
  #Empirical SE
  myresults[4] <- sqrt(sum((output[1,]-mean(output[1,]))^2)/(nrep-1))
  myresults[5] <- sqrt(sum((output[3,]-mean(output[3,]))^2)/(nrep-1))
  myresults[6] <- sqrt(sum((output[5,]-mean(output[5,]))^2)/(nrep-1))
  #MSE
  myresults[7] <- (sum((output[1,]-Rxeff)^2)/(nrep))
  myresults[8] <- (sum((output[3,]-Rxeff)^2)/(nrep))
  myresults[9] <- (sum((output[5,]-Rxeff)^2)/(nrep))
  #Average model SE
  myresults[10] <- sqrt(sum((output[2,])^2)/(nrep))
  myresults[11] <- sqrt(sum((output[4,])^2)/(nrep))
  myresults[12] <- sqrt(sum((output[6,])^2)/(nrep))
  #Rejection percentage
  myresults[13] <-  sum(abs(output[1,])/output[2,]>1.96)
  myresults[14] <-  sum(abs(output[3,])/output[4,]>1.96)
  myresults[15] <-  sum(abs(output[5,])/output[6,]>1.96)
  
  return(myresults)
  
}

BatchedPowerSim_conts <- function(DesMatrix, Kseq, m, Batch, Teffs, Rxeff, ICC, CAC){
  #A function to generate and analyse a single dataset
  #DesMatrix = the design matrix
  #Kseq = number of clusters assigned to each sequence
  #m = number of observations in each cluster in each period
  #Teffs = thee matrix of time effects
  #Rxeff = the treatment effect
  #ICC = intracluster correlation
  #CAC = cluster autocorrelation
  
  #Simulate the data:
  temp <- LongitudinalCRT_onedataset_conts(DesMatrix, Kseq, m, Batch, Teffs, Rxeff, ICC, CAC)
  
  #Fit a model including time by batch effects:
  fit_BbyT <- lmer(YY1i~Xvec +timef*batchf + (1|clusterf) +(1|clusbytimef), temp)
  Rxteststat_BbyT <- abs(fixef(fit_BbyT)[2])/sqrt(vcov(fit_BbyT)[2,2])
  
  #Fit a model including shared time effects only
  fit_T <- lmer(YY1i~Xvec +timef + (1|clusterf) +(1|clusbytimef), temp)
  Rxteststat_T <- abs(fixef(fit_T)[2])/sqrt(vcov(fit_T)[2,2])
  
  #Fit a model including shared time effects and a batch effect
  fit_BT <- lmer(YY1i~Xvec +timef +batchf + (1|clusterf) +(1|clusbytimef), temp)
  Rxteststat_BT <- abs(fixef(fit_BT)[2])/sqrt(vcov(fit_BT)[2,2])
  
  return(c(fixef(fit_BbyT)[2], sqrt(vcov(fit_BbyT)[2,2]), 
           fixef(fit_T)[2], sqrt(vcov(fit_T)[2,2]), 
           fixef(fit_BT)[2],sqrt(vcov(fit_BT)[2,2])))
  
}

LongitudinalCRT_onedataset_conts <- function(DesMatrix, Kseq, m, Batch, Teffs, Rxeff, ICC, CAC){
  #A function to generate one dataset from a longitudinal CRT with design matrix
  #given by DesMatrix
  #Kseq: number of clusters assigned to each sequence of DesMatrix
  #m: number of observations in each cluster in each period
  #Batch: vector indicating which row of DesMatrix belongs to which batch
  #Teffs: the vector of time effects
  #Rxeff: the treatment effect
  #ICC: the intracluster correlation
  #CAC: the cluster autocorrelation. 
  
  #Assume total variance of 1
  sigma_eps2 = 1 -ICC # error variance
  sigmaA2 = CAC*ICC # variance of cluster random effect
  sigmaG2 = ICC*(1-CAC) # variance of cluster-period random effects
  
  #Teffs is a matrix of time effects, of dimension 
  # Ttot columns = total number of periods in the entire design
  # Nseq rows = total number of sequences in the design
  #Some of the elements of Teffs may be NA, when no clusters are observed in that
  #sequence in that period.
  Ttot = ncol(Teffs)
  Nseq = nrow(Teffs)
  #This should be of the same dimension as DesMatrix.
  #Generate a vector for the design and a vector for the time effects
  fulldesmat <- DesMatrix[sort(rep(1:nrow(DesMatrix), Kseq)), ]
  Xvec <- rep(as.vector(t(fulldesmat)), each = m)
  
  fulltimemat <- Teffs[sort(rep(1:nrow(Teffs), Kseq)), ]
  Tvec <- rep(as.vector(t(fulltimemat)), each = m)
  
  #Generate the data, first assuming a complete design
  #This will be reduced to an incomplete design later according to DesMatrix
  #Error terms:
  epsi = rnorm(Nseq*Kseq*Ttot*m,mean=0,sd=sqrt(sigma_eps2)) # epsilon (error term)
  #Cluster random effects:
  randomeffs = rnorm(Nseq*Kseq,mean=0, sd=sqrt(sigmaA2))
  randeffsi <- rep(randomeffs, each=Ttot*m) #one for each participant in each cluster
  #Cluster-period random effects
  CPrandomeffs = rnorm(Nseq*Kseq*Ttot,mean=0, sd=sqrt(sigmaG2))
  CPrandeffsi <- rep(CPrandomeffs, each=m)
  
  clusterVi <- rep(seq(1:(Nseq*Kseq)), each=Ttot*m)
  clusterf = factor(clusterVi)
  
  timeVi <- rep(seq(1:Ttot), each=m)
  timeVi <- rep(timeVi, times=(Kseq*Nseq))
  timef = factor(timeVi)
  
  clusbytime <- rep(1:(Ttot*Kseq*Nseq), each = m)
  clusbytimef = factor(clusbytime)
  
  #indicator for batch
  batch <- rep(Batch, each = m*Ttot*Kseq)
  batchf = factor(batch)
  
  
  #Put everything together to get the outcomes:
  YY1i = Rxeff*Xvec + Tvec + epsi+randeffsi
  Ydf1 = data.frame(YY1i, clusterf, timef, clusbytimef, Xvec, batchf)
  #  Ydf1 = data.frame(YY1i, clusterf, timef, clusbytimef,randeffsi, CPrandeffsi, Xvec, batchf)
  Ydf1 <- Ydf1[!is.na(Xvec),]
  
  return(Ydf1)  
  
  
}

##################################################################################
#Functions for binary outcomes


#A function that generates binary data, fits the models, and returns results
BatchedPowerSim_binary <- function(Ts, nbatch, olap, n, C, icc, cac, treateff, timeeffshared){
  
  #Simulate the data:
  temp <- LongitudinalCRT_onedataset_binaryALT(Ts, nbatch, olap, n, C, icc, cac, treateff, timeeffshared)
  
  #Include in a tryCatch block to allow the simulation to continue despite errors associated with
  #particular datasets:
  #Models fit via GEE
  #Fit a model including time by batch effects (i.e. no shared time effects):
  fit_gee_Tnoshare = NULL
  output1 = tryCatch({
    #NOTE: time in the GEE model must appear as a linear term to match the 
    #sample size calculation using swdpower
    fit_gee_Tnoshare <- geeglm(outcome~ as.factor(treateff) + as.factor(batch)*as.numeric(T_labels), id = C_labels, data=temp, family=binomial(link="logit"), corstr="exchangeable")
    list(fit_gee_Tnoshare, 0)
  },
  warning = function(war){
    print(paste("MY_WARNING1: ", war))
    return(list(fit_gee_Tnoshare,1))    
  },
  error = function(err){
    print(paste("MY_ERROR1: ", err))
    return(list(fit_gee_Tnoshare,2))    
  })
  fit_gee_Tnoshare = output1[[1]]
  status_fit_gee_Tnoshare = output1[[2]]
  
  #Fit a model including separate time effects (i.e. shared effects across batches)
  fit_gee_Tshare = NULL
  output2 = tryCatch({
    #NOTE: time in the GEE model must appear as a linear term to match the 
    #sample size calculation using swdpower
    fit_gee_Tshare <- geeglm(outcome~ as.factor(treateff) + as.numeric(T_labels), id = C_labels, data=temp, family=binomial(link="logit"), corstr="exchangeable")
    list(fit_gee_Tshare, 0)
  },
  warning = function(war){
    print(paste("MY_WARNING1: ", war))
    return(list(fit_gee_Tshare,1))    
  },
  error = function(err){
    print(paste("MY_ERROR1: ", err))
    return(list(fit_gee_Tshare,2))    
  })
  fit_gee_Tshare = output2[[1]]
  status_fit_gee_Tshare = output2[[2]]
  
  #GEE with identity link
  fit_geeid_Tnoshare <- NULL
  output3 = tryCatch({
    fit_geeid_Tnoshare <- geeglm(outcome~ as.factor(treateff) + as.factor(batchtime), id = C_labels, data=temp, family=binomial(link="identity"), corstr="exchangeable")
    list(fit_geeid_Tnoshare, 0)
  },
  warning = function(war){
    print(paste("MY_WARNING1: ", war))
    return(list(fit_geeid_Tnoshare,1))    
  },
  error = function(err){
    print(paste("MY_ERROR1: ", err))
    return(list(fit_geeid_Tnoshare,2))    
  })
  fit_geeid_Tnoshare = output3[[1]]
  status_fit_geeid_Tnoshare = output3[[2]]
  
  #Models fit via LMM (include the categorical effects for periods here)
  # Separate batch by period effects
  fit_lmm_Tnoshare <- NULL
  output4 = tryCatch({
    fit_lmm_Tnoshare <- lmer(outcome~ as.factor(treateff) + as.factor(batchtime) + (1|C_labels), temp)
    list(fit_lmm_Tnoshare, 0)
  },
  warning = function(war){
    print(paste("MY_WARNING1: ", war))
    return(list(fit_lmm_Tnoshare,1))    
  },
  error = function(err){
    print(paste("MY_ERROR1: ", err))
    return(list(fit_lmm_Tnoshare,2))    
  })
  fit_lmm_Tnoshare = output4[[1]]
  status_fit_lmm_Tnoshare = output4[[2]]
  
  #LMM with shared period effects across batches
  fit_lmm_Tshare<- NULL
  output5 = tryCatch({
    fit_lmm_Tshare <- lmer(outcome~ as.factor(treateff) + as.factor(T_labels) + (1|C_labels), temp)
    list(fit_lmm_Tshare, 0)
  },
  warning = function(war){
    print(paste("MY_WARNING1: ", war))
    return(list(fit_lmm_Tshare,1))    
  },
  error = function(err){
    print(paste("MY_ERROR1: ", err))
    return(list(fit_lmm_Tshare,2))    
  })
  fit_lmm_Tshare = output5[[1]]
  status_fit_lmm_Tshare = output5[[2]]
  
  rm(output1)
  rm(output2)
  rm(output3)
  rm(output4)
  rm(output5)
  
  if (is.null(fit_gee_Tnoshare)){
    geeTnosharepval <- NA
  }
  else {
    geeTnosharepval <- coef(summary(fit_gee_Tnoshare))[2,4] < 0.05
  }
  
  if (is.null(fit_gee_Tshare)){
    geeTsharepval <- NA
  }
  else {
    geeTsharepval <- coef(summary(fit_gee_Tshare))[2,4] < 0.05
  }
  
  if (is.null(fit_geeid_Tnoshare)){
    geeidTnosharepval <- NA
  }
  else {
    geeidTnosharepval <- coef(summary(fit_geeid_Tnoshare))[2,4] < 0.05  
  }
  
  if (is.null(fit_lmm_Tnoshare)){
    lmmTnosharepval <- NA
  }
  else {
    lmmTnosharepval <- abs(fixef(fit_lmm_Tnoshare)[2])/sqrt(vcov(fit_lmm_Tnoshare)[2,2])> 1.96  
  }
  
  if (is.null(fit_lmm_Tshare)){
    lmmTsharepval <- NA
  }
  else {
    lmmTsharepval <- abs(fixef(fit_lmm_Tshare)[2])/sqrt(vcov(fit_lmm_Tshare)[2,2])> 1.96 
  }
  
  rm(fit_gee_Tnoshare)
  rm(fit_gee_Tshare)
  rm(fit_geeid_Tnoshare)
  rm(fit_lmm_Tnoshare)
  rm(fit_lmm_Tshare)
  
  return(c(geeTsharepval, geeTnosharepval, lmmTsharepval, lmmTnosharepval, geeidTnosharepval))
  
}


BatchedPowerSim_bin_wrap <- function(nrep, Ts, nbatch, olap, n, C, icc, cac, treateff, timeeffshared){
  #A wrapper for the BatchedPowerSim_binary function
  
  output <- replicate(nrep, BatchedPowerSim_binary(Ts, nbatch, olap, n, C, icc, cac, treateff, timeeffshared))
  
  #4 models have been fit to the simulated datasets:
  #1. GEE with shared time effects
  #2. GEE with batch by time effects
  #3. LMM with shared time effects
  #4. LMM with batch by time effects
  #For each, want the following:
  #- Rejection percentage
  
  myresults <- NULL
  
  #Rejection percentage
  myresults[1] <-  sum(output[1,], na.rm=TRUE)
  myresults[2] <- sum(is.na(output[1,]))
  myresults[3] <-  sum(output[2,], na.rm=TRUE)
  myresults[4] <- sum(is.na(output[2,]))
  myresults[5] <-  sum(output[3,], na.rm=TRUE)
  myresults[6] <- sum(is.na(output[3,]))
  myresults[7] <-  sum(output[4,], na.rm=TRUE)
  myresults[8] <- sum(is.na(output[4,]))
  myresults[9] <-  sum(output[5,], na.rm=TRUE)
  myresults[10] <- sum(is.na(output[5,]))
  
  
  
  return(myresults)
  
}



#Wrapper function for Fan Li's binGEN function, to get the data
#in the correct format
binGEN_wrap <- function(n, m, t, delta, beta, alpha){
  y<-binGEN(n,m,t,delta,beta,alpha)
  y<-c(y)
  
  # marginal mean design matrix including period and treatment indicators
  X<-NULL
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  g<-n/(t-1) # number of clusters per step
  for(i in 1:(t-1)){
    for(j in 1:g){
      X<-rbind(X,kronecker(cbind(diag(t),trtSeq[i,]),rep(1,m)))}
  }
  cluster<-rep(1:n,each=t*m)         # create cluster id
  ind<-rep(rep(1:m,t),n)             # create individual id
  period<-rep(rep(1:t,each=m),n)     # create period label
  
  simdata_bin<-data.frame(cbind(y,ind,cluster,period,X))
  
  return(simdata_bin)
  
}

#A function that acts as a wrapper to binGEN_wrap to generate data from a batched SW
LongitudinalCRT_onedataset_binaryALT <- function(Ts, nbatch, olap, n, C, icc, cac, treateff, timeeffshared=0){
  #Ts = number of periods in the component blocks of the batched SW
  #nbatch = number of batches of SWs
  #olap = number of periods of overlap between consecutive stepped wedges
  #n = number of observations in each cluster in each  period
  #C = number of clusters assigned to each sequence
  #icc = within-period ICC
  #cac = between-period ICC
  #treateff = treatment effect
  #timeeffshared = 1 if time effects are shared across batches; 0 otherwise
  #ONLY ALLOW UNSHARED EFFECTS HERE
  
  #A base matrix for the time effects for a single batch. 
  timeeffmatBASE <- matrix(data=seq(0.4, 0.39+(Ts)*0.01, 0.01), nrow=(Ts-1), ncol=Ts, byrow =TRUE)
  deltaeff <- log((treateff + 0.39+(Ts)*0.01)/(1- treateff - 0.39-(Ts)*0.01))  -log((0.39+(Ts)*0.01)/(1-(0.39+(Ts)*0.01)))
  
  #Generate the design matrix for a single batch:
  desmat <- SWdesmat(Ts)
  
  #Generate data batch by batch
  mydata <- NULL
  for(b in 1:nbatch){
    dessmat <- SWdesmat(Ts)
    if(timeeffshared == 1){
      timeeffmat <- timeeffmatBASE + (b-1)*(Ts-olap)*0.01
    }
    else if(timeeffshared == 0){
      timeeffmat <-timeeffmatBASE
    }
    
    timeffmat_logit <- log(timeeffmatBASE[1,]/(1-timeeffmatBASE[1,]))
    
    nextbatch <- binGEN_wrap(n = C*(Ts-1), m=n, t=Ts, delta=deltaeff, beta=timeffmat_logit, alpha=c(icc, cac*icc, 0))
    
    #Need to make sure that the cluster labels and time labels are correct
    lob <- ((b-1)*(Ts-1)*C +1)
    hib <- (b*(Ts-1)*C) 
    nextbatch$C_labels <- rep(seq(from=lob, to=hib), each = Ts*n )
    nextbatch$T_labels <- nextbatch$period
    if(b>1){
      timetemp <- rep((b-1)*(Ts-olap), each=n)
      nextbatch$T_labels <-  nextbatch$T_labels + rep(timetemp, each = (Ts-1)*C)  
    }
    
    #Batch indicator:
    nextbatch$batch <- rep(b, times=nrow(nextbatch))
    mydata <- rbind(mydata, nextbatch)   
    
  }
  
  #A corrected cluster by time indicator 
  clusbytime <- rep(1:(Ts*(Ts-1)*C*nbatch), each = n)
  mydata$clusbytime <- clusbytime
  
  #include the treatment effect indicator 
  desmatexp <- desmat[rep(1:nrow(desmat), each = C), ]
  treateffind <- rep(as.vector(t(desmatexp)), each = n)
  mydata$treateff <- rep(treateffind, times = nbatch)
  
  #Batch by time factor
  mydata$batchtime <- with(mydata, interaction(as.factor(batch),  as.factor(T_labels)))
  mydata$batchtime <- factor(as.numeric( mydata$batchtime) - 1)
  
  mydata$outcome <- mydata$y
  mydata$treateff <- mydata$V11  
  return(mydata)
  
}
