#Terrestrial Field Ecology
#Snail Analysis

library(here)
here::i_am("Snails.R")

################################################################################
#Part 1: Reading in Data, Performing Calculations
################################################################################

Snail.Data <- read.table("Snail.Data.txt", sep="\t", header=TRUE)
Snail.Groups <- read.table("Snail.Groups.txt", sep="\t", header=TRUE)
d.detailed <- Snail.Data
native1 <- read.table("native1.txt", sep="\t", header=TRUE)
native2 <- read.table("native2.txt", sep="\t", header=TRUE)
nonnative1 <- read.table("nonnative1.txt", sep="\t", header=TRUE)
nonnative2 <- read.table("nonnative2.txt", sep="\t", header=TRUE)
mixed1 <- read.table("mixed1.txt", sep="\t", header=TRUE)
mixed2 <- read.table("mixed2.txt", sep="\t", header=TRUE)

#Set up summary matrix
first.rows.by.snail <- seq(from=1, to=807, by=62)
first.rows.by.trial <- seq(from=1, to=838, by=31)
summary <- Snail.Data[first.rows.by.snail, 1:4]

#Compute detailed info on movement
#distance from origin
d.detailed$d.frm.origin <- ((d.detailed$x)^2+(d.detailed$y)^2)^0.5
#compute distance this step
d.detailed$d.this.step <- numeric(868)
for (i in 2:868){
  d.detailed$d.this.step[i] <- ((d.detailed$x[i]-d.detailed$x[(i-1)])^2 + 
    (d.detailed$y[i]-d.detailed$y[(i-1)])^2)^0.5
}
d.detailed$d.this.step[first.rows.by.trial] <- 0
#compute total distance this trial
d.detailed$total <- numeric(868)
for (i in first.rows.by.trial){
  for (j in 0:30){
    d.detailed$total[i+j] <- sum(d.detailed$d.this.step[i:(i+j)])
  }
}

#Update the summary data with total distance covered on wet and dry environments
summary$total.dry <- d.detailed$total[(first.rows.by.snail+30)]
summary$total.wet <- d.detailed$total[(first.rows.by.snail+61)]




#COLOR CODING

#Native snails used in the group trials
natives <- c("8","5","7","B")
colnatives <- c("cyan","blue","green","purple")

#Nonnative snails used in the group trials
nonnatives <- c("6","1","11","2")
colnonnatives <- c("red","orange","brown","black")

#All snails used in group trials
group <- c(natives, nonnatives)
colmixed <- c(colnatives, colnonnatives)

#CALCULATIONS FOR MIXED

for (i in 0:7){ #8 snails
  for (j in 1:31){ #31 time steps, including time t=0
    #distance from origin <- (this point x - start point x)^2 ...
    mixed1$d.origin[31*i+j] <- ((mixed1$y[31*i+j]-mixed1$y[31*i+1])^2 + (mixed1$x[31*i+j]-mixed1$x[31*i+1])^2)^0.5
    #distance this step
    if (j == 1){
      mixed1$d.this.step[31*i+j] <- 0
    } else {
      mixed1$d.this.step[31*i+j] <- ((mixed1$y[31*i+j]-mixed1$y[31*i+j-1])^2 + (mixed1$x[31*i+j]-mixed1$x[31*i+j-1])^2)^0.5
    }
    # average distance to conspecifics:
    
    dist.marker <- 1
    if (is.element(mixed1$Snail.ID[31*i+j], natives)){
      snail.set.con <- natives[!natives==mixed1$Snail.ID[31*i+j]]
      snail.set.hetero <- nonnatives
    } else {
      snail.set.con <- nonnatives[!nonnatives==mixed1$Snail.ID[31*i+j]]
      snail.set.hetero <- natives
    }
    distances.con <- numeric(3)
    distances.hetero <- numeric(4)
    for (ID in snail.set.con){
      sub <- subset(mixed1, Snail.ID == ID)
      distances.con[dist.marker] <- ((mixed1$y[31*i+j]-sub$y[j])^2 + (mixed1$x[31*i+j]-sub$x[j])^2)^0.5
      dist.marker = dist.marker + 1
    } 
    mixed1$avg.dist.con[31*i+j] <- mean(distances.con)
    # distance to nearest conspecific:
    mixed1$min.dist.con[31*i+j] <- min(distances.con)
    
    # average distance to heterospecifics:
    dist.marker <- 1
    for (ID in snail.set.hetero){
      sub <- subset(mixed1, Snail.ID == ID)
      distances.hetero[dist.marker] <- ((mixed1$y[31*i+j]-sub$y[j])^2 + (mixed1$x[31*i+j]-sub$x[j])^2)^0.5
      dist.marker = dist.marker + 1
    }
    mixed1$avg.dist.hetero[31*i+j] <- mean(distances.hetero)
    #distance to nearest heterospecific:
    mixed1$min.dist.hetero[31*i+j] <- min(distances.hetero)
  }
}




for (i in 0:7){ #8 snails
  for (j in 1:31){ #31 time steps, including time t=0
    #distance from origin <- (this point x - start point x)^2 ...
    mixed2$d.origin[31*i+j] <- ((mixed2$y[31*i+j]-mixed2$y[31*i+1])^2 + (mixed2$x[31*i+j]-mixed2$x[31*i+1])^2)^0.5
    #distance this step
    if (j == 1){
      mixed2$d.this.step[31*i+j] <- 0
    } else {
      mixed2$d.this.step[31*i+j] <- ((mixed2$y[31*i+j]-mixed2$y[31*i+j-1])^2 + (mixed2$x[31*i+j]-mixed2$x[31*i+j-1])^2)^0.5
    }
    # average distance to conspecifics:
    
    dist.marker <- 1
    if (is.element(mixed2$Snail.ID[31*i+j], natives)){
      snail.set.con <- natives[!natives==mixed2$Snail.ID[31*i+j]]
      snail.set.hetero <- nonnatives
    } else {
      snail.set.con <- nonnatives[!nonnatives==mixed2$Snail.ID[31*i+j]]
      snail.set.hetero <- natives
    }
    distances.con <- numeric(3)
    distances.hetero <- numeric(4)
    for (ID in snail.set.con){
      sub <- subset(mixed2, Snail.ID == ID)
      distances.con[dist.marker] <- ((mixed2$y[31*i+j]-sub$y[j])^2 + (mixed2$x[31*i+j]-sub$x[j])^2)^0.5
      dist.marker = dist.marker + 1
    } 
    mixed2$avg.dist.con[31*i+j] <- mean(distances.con)
    # distance to nearest conspecific:
    mixed2$min.dist.con[31*i+j] <- min(distances.con)
    
    # average distance to heterospecifics:
    dist.marker <- 1
    for (ID in snail.set.hetero){
      sub <- subset(mixed2, Snail.ID == ID)
      distances.hetero[dist.marker] <- ((mixed2$y[31*i+j]-sub$y[j])^2 + (mixed2$x[31*i+j]-sub$x[j])^2)^0.5
      dist.marker = dist.marker + 1
    }
    mixed2$avg.dist.hetero[31*i+j] <- mean(distances.hetero)
    #distance to nearest heterospecific:
    mixed2$min.dist.hetero[31*i+j] <- min(distances.hetero)
  }
}


#SIZE AND DISTANCE

plot(summary$total.dry ~ summary$size)
plot(summary$total.wet ~ summary$size)

native.distances <- summary$total.dry[summary$species=="native"] 
native.distances[7:12] <- summary$total.wet[summary$species=="native"]

non.distances <- summary$total.dry[summary$species=="non"] 
non.distances[9:16] <- summary$total.wet[summary$species=="non"]


t.test(native.distances, non.distances)
summary$mixed1.tot.d <- 0
summary$mixed1.tot.d[14] <-sum(mixed1$d.this.step[156:186]) #snail 2
summary$mixed1.tot.d[9] <- sum(mixed1$d.this.step[32:62]) #snail 6
summary$mixed1.tot.d[2] <- sum(mixed1$d.this.step[94:124]) #snail 1
summary$mixed1.tot.d[6] <- sum(mixed1$d.this.step[218:248]) #snail 11

summary$mixed1.tot.d[13] <- sum(mixed1$d.this.step[125:155]) #snail 8
summary$mixed1.tot.d[12] <- sum(mixed1$d.this.step[63:93]) #snail 5
summary$mixed1.tot.d[5] <- sum(mixed1$d.this.step[187:217]) #snail 7
summary$mixed1.tot.d[3] <- sum(mixed1$d.this.step[1:31]) #snail B

summary$mixed2.tot.d[14] <-sum(mixed2$d.this.step[156:186]) #snail 2
summary$mixed2.tot.d[9] <- sum(mixed2$d.this.step[32:62]) #snail 6
summary$mixed2.tot.d[2] <- sum(mixed2$d.this.step[94:124]) #snail 1
summary$mixed2.tot.d[6] <- sum(mixed2$d.this.step[218:248]) #snail 11

summary$mixed2.tot.d[13] <- sum(mixed2$d.this.step[125:155]) #snail 8
summary$mixed2.tot.d[12] <- sum(mixed2$d.this.step[63:93]) #snail 5
summary$mixed2.tot.d[5] <- sum(mixed2$d.this.step[187:217]) #snail 7
summary$mixed2.tot.d[3] <- sum(mixed2$d.this.step[1:31]) #snail B

personality <- subset(summary, Snail.ID == c(8,5,7,6,1,11,2))

personality <- subset(summary, Snail.ID == c(8,5,7,6,1,11,2))



################################################################################
#Part 2: Plots for Presentation
################################################################################



################################################################################
#Part 3: Metacommunity Simulations
################################################################################

# See "Metacommunity_Simulations.R"

################################################################################
#Part 4: Analysis of Simulations
################################################################################

# See "Analysis_of_Simulations.R"