################################################################################
#
# Code to accompany: "The batched stepped wedge design: a design robust to 
#                     delays in cluster recruitment" 
#                     by J Kasza, R Bowden, R Hooper, AB Forbes.
#
# Questions or comments to J Kasza, jessica.kasza@monash.edu
#
#
# BatchedSW File 2: code to run the simulation study 
#
# 2022-04-11
################################################################################

source("File01_BatchSW.R")


################################################################################
#Simulation of continuous outcomes
#Set-up for the continuous simulations
Tlist <- c(6)
nbatchlist <- c(2)
Klist <- c(4)
olaplist <- c(5,4,3,2,1,0)
mlist <- c(10)
icclist <- c(0.01, 0.05, 0.1)
caclist <- c(1, 0.95, 0.75)
sharedtimeeffs <- c(0, 1)
effsize <- c(0,0.15)

#Generate a matrix with all combos:
allcombos <- expand.grid(Tlist, nbatchlist, olaplist, Klist, mlist, sharedtimeeffs, effsize, icclist, caclist)
contsresults <- matrix(data= NA, nrow=nrow(allcombos), ncol = 15)

start_time <- Sys.time()
set.seed(3728)
for(i in 1:nrow(allcombos)){
  nrep<- 1000
  contsresults[i,] <- BatchedPowerSim_conts_wrap(nrep, allcombos[i,1], allcombos[i,2], allcombos[i,3],
                                                 allcombos[i,4], allcombos[i,5], allcombos[i,6], 
                                                 allcombos[i,7], allcombos[i,8],allcombos[i,9])
  
}
end_time <- Sys.time()
end_time - start_time

write.csv(contsresults, file="contsresults.csv")


################################################################################
#Simulation of binary outcomes
#Set-up for the binary simulations

Tlist <- c(6)
nbatchlist <- c(2)
Klist <- c(20) 
olaplist <- c(5,4,3,2,1,0)
mlist <- c(10)
icclist <- c(0.01, 0.05)
caclist <- c(1) 
sharedtimeeffs <- 0 
effsize <- c(0, 0.025)

allcombosBIN <- expand.grid(Tlist, nbatchlist, olaplist, Klist, mlist, sharedtimeeffs, effsize, icclist, caclist)
binaryresults_new <- matrix(data= NA, nrow=nrow(allcombosBIN), ncol = 10)

start_time <- Sys.time()
set.seed(7878)
for(i in 1:20){
  nrep<- 1000
  print(allcombosBIN[i,])
  binaryresults_new[i,] <- BatchedPowerSim_bin_wrap(nrep, allcombosBIN[i,1], allcombosBIN[i,2], allcombosBIN[i,3],
                                                    allcombosBIN[i,5], allcombosBIN[i,4], allcombosBIN[i,8], 
                                                    allcombosBIN[i,9], allcombosBIN[i,7] ,allcombosBIN[i,6])
  print(binaryresults_new[i,])
}

start_time <- Sys.time()
set.seed(3636)
for(i in 21:nrow(allcombosBIN)){
  nrep<- 1000
  print(allcombosBIN[i,])
  binaryresults_new[i,] <- BatchedPowerSim_bin_wrap(nrep, allcombosBIN[i,1], allcombosBIN[i,2], allcombosBIN[i,3],
                                                    allcombosBIN[i,5], allcombosBIN[i,4], allcombosBIN[i,8], 
                                                    allcombosBIN[i,9], allcombosBIN[i,7] ,allcombosBIN[i,6])
  print(binaryresults_new[i,])
}
end_time <- Sys.time()
end_time - start_time


#Theoretical power:
binaryresults_theor <- matrix(data= NA, nrow=nrow(allcombosBIN), ncol = 1)
testdesign <- batchSWscheme(6, 1, 6, 40)
for(i in 1:nrow(allcombosBIN)){
  binaryresults_theor[i,] <- swdpower(K=allcombosBIN[i,5], design = testdesign, family = "binomial", model="marginal",
                                      link="logit", type="cross-sectional",
                                      meanresponse_start = 0.4, meanresponse_end0 = 0.45,
                                      meanresponse_end1 = 0.45+allcombosBIN[i,7], typeIerror = 0.05,
                                      alpha0 = allcombosBIN[i,8], alpha1= allcombosBIN[i,8]*allcombosBIN[i,9])$Power
  
}

#only want the GEE results
binary_results <- cbind(allcombosBIN, binaryresults_new[,3], binaryresults_theor)
write.csv(binary_results, file="binaryresults_202204.csv")


