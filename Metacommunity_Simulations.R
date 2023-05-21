
#This script codes the metacommunity simulation and performs Bayesian analysis
#of the data resulting from the simulation.

################################################################################
#Part 1: Simulations
################################################################################

#Parameters
birth.rate.nat <- 0.5 #Birth rate of natives
birth.rate.non <- 0.5 #Birth rate of nonnatives
death.rate.non <- 0.1 #Death rate of nonnatives
death.rate.nat <- 0.1 #Death rate of natives
nat.wet.rate <- 1.6560/10 #Probability of movement of natives in wet environment
nat.dry.rate <-3.9031/10 #Probability of movement of natives in dry environment
non.wet.rate <- 6.3341/10 #Probability of movement of nonnatives in wet environment
non.dry.rate <- 8.2927/10 #Probability of movement of nonnatives in dry environment

#Starting values in starting patches
avg.natives <- 10 #Average initial value of natives in occupied patches
sd.natives <-2 #standard deviation of initial value of natives in occupied patches
avg.exotics <- 10 #Average initial value of nonnatives in occupied patches
sd.exotics <- 2 #standard deviation of initial value of nonnatives in occupied patches
num.patches.occupied <- 10 #How many of the 20 patches start out occupied

#Snail Simulation
#Keep track of all of the trials performed in wet conditions
All.wet <- data.frame(trial=numeric(0),patch=numeric(0),natives=numeric(0),exotics=numeric(0),total=numeric(0))
#Keep track of all of the trials performed in dry conditions
All.dry <- data.frame(trial=numeric(0),patch=numeric(0),natives=numeric(0),exotics=numeric(0),total=numeric(0))


for (replicate in 1:100){ #Perform 100 replicates of both dry and wet
  
  for (env in 1:2){#1 is wet environment, 2 is dry environment
    
    if (env==1){ #wet environment
      nat.rate <- nat.wet.rate
      non.rate <- non.wet.rate
    } else{ #dry environment 
      nat.rate <- nat.dry.rate
      non.rate <- non.dry.rate
    }
    
    #A data.frame to keep track of natives and exotics in each patch
    patches <- data.frame(patch=numeric(20),natives=numeric(20),exotics=numeric(20),total=numeric(20))
    patches$patch <- 1:20 # 20 patches)
    patches$natives <- 0
    patches$exotics <- 0
    
    #get patches occupied
    list.patches.i <- numeric(num.patches.occupied)
    for (i in 1:num.patches.occupied){
      proposed.patch.i <- sample(1:20,1)
      while (is.element(proposed.patch.i,list.patches.i)){
        proposed.patch.i <- sample(1:20,1)
      }
      list.patches.i[i] <- proposed.patch.i
    }
    
    #populate those random patches
    patches$natives[list.patches.i] <- round(rnorm(num.patches.occupied,avg.natives,sd.natives))
    patches$exotics[list.patches.i] <- round(rnorm(num.patches.occupied,avg.exotics,sd.exotics))
    patches$total <- patches$natives + patches$exotics
    
    while (sum(patches$total) < (100*20)){#Keep going until all patches at carrying capacity
      for (i in 1:20){
        #Native dynamics
        for (j in 1:patches$natives[i]){
          #Births
          if (patches$total[i] < 100 && runif(1)< birth.rate.nat){
            patches$natives[i] = patches$natives[i]+1
            patches$total[i] <- patches$natives[i] + patches$exotics[i]
          } #end Births loop
          #Deaths
          if (patches$natives[i]>0 && runif(1) < death.rate.nat){
            patches$natives[i] = patches$natives[i]-1
            patches$total[i] <- patches$natives[i] + patches$exotics[i]
          } #end Deaths loop
          #Movement
          k=1
          candidate.recips <- numeric(0)
          for (other.patch in 1:20){  #Get candidate recipient patches
            if (patches$total[i] > patches$total[other.patch]){
              candidate.recips[k] <- other.patch
              k = k+1
            }
          } #End getting candidate recipient patches
          if (runif(1) < nat.rate && length(candidate.recips)>0){ #Movement
            recipient.patch <- sample(candidate.recips,1)
            if (patches$total[recipient.patch] < 100 && patches$natives[i]>0){
              patches$natives[i] <- patches$natives[i]-1
              patches$natives[recipient.patch] <- patches$natives[recipient.patch]+1
              patches$total <- patches$natives + patches$exotics
            }
          } # End movement
        } #End native dynamics
        
        #Non-native Dynamics
        for (j in 1:patches$exotics[i]){
          #Births
          if (patches$total[i] < 100 && runif(1)< birth.rate.non){
            patches$exotics[i] = patches$exotics[i]+1
            patches$total[i] <- patches$natives[i] + patches$exotics[i]
          }
          #Deaths
          if (patches$exotics[i]>0 && runif(1) < death.rate.non){
            patches$exotics[i] = patches$exotics[i]-1
            patches$total[i] <- patches$natives[i] + patches$exotics[i]
          }
          #Movement
          k=1
          candidate.recips <- numeric(0)
          for (other.patch in 1:20){
            if (patches$total[i] > patches$total[other.patch]){
              candidate.recips[k] <- other.patch
              k = k+1
            }
          }
          if (runif(1) < non.rate && length(candidate.recips)>0){
            recipient.patch <- sample(candidate.recips,1)
            if (patches$total[recipient.patch] < 100 && patches$exotics[i]>0){
              patches$exotics[i] <- patches$exotics[i]-1
              patches$exotics[recipient.patch] <- patches$exotics[recipient.patch]+1
              patches$total <- patches$natives + patches$exotics
            }
          }
        }#End nonnative loop
        
        
        patches$total <- patches$natives + patches$exotics
      } #End going through all of the patches
    } #End of the while loop
    
    if (env==1){
      All.wet[(1+(replicate-1)*20):(20+(replicate-1)*20),1] <- replicate
      All.wet[(1+(replicate-1)*20):(20+(replicate-1)*20),2:5] <- patches
    } else {
      All.dry[(1+(replicate-1)*20):(20+(replicate-1)*20),1] <- replicate
      All.dry[(1+(replicate-1)*20):(20+(replicate-1)*20),2:5] <- patches
    }
  }#End different types of environments
  
}#End of the replicate

################################################################################
#Part 2: Analysis of Simulations
################################################################################

#This analysis requires Rtools

#For more information on installing rstan and rethinking, see
#https://github.com/rmcelreath/rethinking

library(rethinking)
library(rstan)
library(ggplot2)

par(mfrow=c(2,2), mai = c(1, 1, 0.1, 0.1)) 

hist(All.wet$natives, main="", xlab="", ylim=c(0,400), xlim=c(0,100), ylab="Frequency")
title("A) Natives, wet",line=-1, cex.main=1, adj=0.05 )
abline(v = mean(All.wet$natives), col = "darkgray", lwd = 2)
abline(v = 50, col = "darkgray", lwd = 2, lty="dotted")
box(which = "plot", lty = "solid")

hist(All.wet$exotics, main="", xlab="", ylim=c(0,400), xlim=c(0,100), ylab="")
title("B) Exotics, wet",line=-1, cex.main=1, adj=0.05 )
abline(v = mean(All.wet$exotics), col = "darkgray", lwd = 2)
abline(v = 50, col = "darkgray", lwd = 2, lty="dotted")
box(which = "plot", lty = "solid")

hist(All.dry$natives, main="", xlab="Number of natives in the patch", ylim=c(0,400), xlim=c(0,100), ylab="Frequency")
title("C) Natives, dry",line=-1, cex.main=1, adj=0.05 )
abline(v = mean(All.dry$natives), col = "darkgray", lwd = 2)
abline(v = 50, col = "darkgray", lwd = 2, lty="dotted")
box(which = "plot", lty = "solid")

hist(All.dry$exotics, main="", xlab="Number of exotics in the patch", ylim=c(0,400), xlim=c(0,100), ylab="")
title("D) Exotics, dry",line=-1, cex.main=1, adj=0.05 )
abline(v = mean(All.dry$exotics), col = "darkgray", lwd = 2)
abline(v = 50, col = "darkgray", lwd = 2, lty="dotted")
box(which = "plot", lty = "solid")

require(graphics)

ks.test(All.wet$natives, All.wet$exotics)



hist(All.dry$natives, main="Distribution of natives in patches", xlab="Number of natives in the patch")
hist(All.dry$exotics, main="Distribution of nonnatives in patches", xlab="Number of nonnatives in the patch")


require(graphics)
ks.test(All.dry$natives, All.dry$exotics)

############################################
#Combine all data
samples <- length(All.dry$natives) #to figure out how many samples there are
Data <- data.frame(species=numeric(4*samples),environment=numeric(4*samples),population=numeric(4*samples))

Data$species[1:samples] <- 0 #natives
Data$environment[1:(2*samples)] <- 0 #wet
Data$population[1:samples] <- All.wet$natives
Data$species[(samples+1):(2*samples)] <- 1 #nonnatives
Data$population[(samples+1):(2*samples)] <- All.wet$exotics
Data$environment[(2*samples + 1):(4*samples)] <- 1 #dry
Data$species[(2*samples+1):(3*samples)] <- 0 #natives
Data$population[(2*samples+1):(3*samples)] <- All.dry$natives
Data$species[(3*samples+1):(4*samples)] <- 1 #nonnatives
Data$population[(3*samples+1):(4*samples)] <- All.dry$exotics

###################################
#The effect of species, environment and their interaction
m1 <- ulam(
  alist(
    population ~ dnorm( mu , sigma ) ,
    mu <- a + be*environment + bs*species + bes*environment*species ,
    a ~ dnorm( 0 , 10 ) ,
    be ~ dnorm( 0 , 10 ) ,
    bs ~ dnorm( 0 , 10 ) ,
    bes ~ dnorm( 0 , 10 ) ,
    sigma ~ dcauchy( 0 , 1 )
  ) ,
  data=Data, log_lik=TRUE,
  start=list(a=0,be=0,bs=0, bes=0, sigma=1) ,
  warmup=1000 , iter=5000 )

###################################
#The effect of species and environment
m2 <- ulam(
  alist(
    population ~ dnorm( mu , sigma ) ,
    mu <- a + be*environment + bs*species,
    a ~ dnorm( 0 , 10 ) ,
    be ~ dnorm( 0 , 10 ) ,
    bs ~ dnorm( 0 , 10 ) ,
    sigma ~ dcauchy( 0 , 1 )
  ) ,
  data=Data, log_lik=TRUE,
  start=list(a=0,be=0,bs=0, sigma=1) ,
  warmup=1000 , iter=5000 )

#############################################
#Just the effect of environment
m3 <- ulam(
  alist(
    population ~ dnorm( mu , sigma ) ,
    mu <- a + be*environment,
    a ~ dnorm( 0 , 10 ) ,
    be ~ dnorm( 0 , 10 ) ,
    sigma ~ dcauchy( 0 , 1 )
  ) ,
  data=Data, log_lik=TRUE,
  start=list(a=0,be=0, sigma=1) ,
  warmup=1000 , iter=5000 )

###########################################
#Just the effect of species
m4 <- ulam(
  alist(
    population ~ dnorm( mu , sigma ) ,
    mu <- a + bs*species,
    a ~ dnorm( 0 , 10 ) ,
    bs ~ dnorm( 0 , 10 ) ,
    sigma ~ dcauchy( 0 , 1 )
  ) ,
  data=Data, log_lik=TRUE,
  start=list(a=0,bs=0, sigma=1) ,
  warmup=1000 , iter=5000 )
############################################
par(mfrow=c(1,1), mar = c(1, 1, 1, 1))

#Diagnose the MCMC
plot(m1) #interaction model
plot(m2) #env + sp
plot(m3) #env
plot(m4) #sp

#Parameter estimates
precis(m1)
precis(m2)
precis(m3)
precis(m4)

#Compare models
compare(m1,m2,m3,m4)

pairs(m1)

post <- extract.samples(m1)

##############################
m5 <- map(
  alist(
    population ~ dnorm( mu , sigma ) ,
    mu <- a + be*environment + bs*species + bes*environment*species ,
    a ~ dnorm( 40 , 10 ) ,
    be ~ dnorm( 5 , 10 ) ,
    bs ~ dnorm( 20 , 10 ) ,
    bes ~ dnorm( -10 , 10 ) ,
    sigma ~ dcauchy( 12 , 1 )
  ) ,
  data=Data, 
  start=list(a=40,be=-5,bs=20, bes=-10, sigma=12) )


############################
#Plot of model
par(mfrow=c(1,1))
# get minimum and maximum rugged values
q.species <- range(Data$species)
# compute lines and confidence intervals
mu.specieslo <- link( m1 , data=list(species=c(0,0),environment=(0:1) ))
mu.specieslo.mean <- apply( mu.specieslo , 2 , mean )
mu.specieslo.PI <- apply( mu.specieslo , 2 , PI )
mu.specieshi <- link( m1 , data=list(species=c(1,1),environment=0:1) )
mu.specieshi.mean <- apply( mu.specieshi , 2 , mean )
mu.specieshi.PI <- apply( mu.specieshi , 2 , PI )
# plot it all, splitting points at mean
med.r <- mean(Data$species)
ox <- ifelse( Data$species > med.r , 0.05 , -0.05 )
plot( Data$environment + ox , Data$population ,
      col=ifelse(Data$species>med.r,rangi2,"black") ,
      xlim=c(-0.25,1.25) , xaxt="n" , ylab="Population" ,
      xlab="Environment" )
axis( 1 , at=c(0,1) , labels=c("wet","dry") )
lines( 0:1 , mu.specieslo.mean , lty=2 )
shade( mu.specieslo.PI , 0:1 )
lines( 0:1 , mu.specieshi.mean , col=rangi2 )
shade( mu.specieshi.PI , 0:1 , col=col.alpha(rangi2,0.25) )



