#===========================================
# SEAGRASS METABOLISM LITERATURE SYNTHESIS
#===========================================
# DATA SOURCE: Johnson et al. (2017) Scientific Reports
# Manuscript & Data: https://www.nature.com/articles/s41598-017-13142-4
#---------------------------
# Script written by G.M. Wilkinson (March 2022) for ZOO 955 as an example of meta-analysis approaches. Analysis intended for demonstration purposes only
# DO NOT USE FOR STATISTICAL INFERENCE - consult original manuscript and analyses therein for scientific conclusions from this literature synthesis
#---------------------------

# Packages =================================
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)
if (!require(rjags)) install.packages('rjags')
library(rjags)
if (!require(coda)) install.packages('coda')
library(coda)

# NOTE!!!!!!!!!!!
# Must have JAGS (just another gibbs sampler) installed on your machine 
# https://sourceforge.net/projects/mcmc-jags/files/

#===========================================
# DATA
data = read.csv("Johnson_2021_Seagrass_LiteratureSynthesis.csv") %>%
  # Alter the method categories to create 3 categories of similar methods
  mutate(method = replace(method, 
                          method=="in situ bell jar (incubation chamber)", 
                          "incubation chamber"),
         method = replace(method, 
                          method=="in situ benthic incubation chamber", 
                          "incubation chamber"),
         method = replace(method, 
                          method=="in situ diel curve of O2 concentration", 
                          "diel O2"),
         method = replace(method, 
                          method=="eddy correlation", 
                          "eddy flux"),
         method = replace(method, 
                          method=="eddy covariance", 
                          "eddy flux")) %>%
  #rename net community production as net ecosystem production
  rename(NEP_mmolCm2d = NCP_mmolCm2d,
         NEP_sd = NCP_sd) %>%
  # filter(!(speciesCat=="Cn" | speciesCat=="Ea" | speciesCat=="Hw")) %>%
  mutate(speciesCat = replace(speciesCat, speciesCat=="Zn", "Zostera"),
         speciesCat = replace(speciesCat, speciesCat=="Zm", "Zostera"),
         speciesCat = replace(speciesCat, speciesCat=="Po", "Posidonia"),
         speciesCat = replace(speciesCat, speciesCat=="Tt", "Thalassia"),
         speciesCat = replace(speciesCat, speciesCat=="mix", "Mixed"))

#Let's learn something about the data ======
unique(data$Source) #29 studies
nrow(data) #58 observations
summary(data$NEP_sd) #10 missing observations of NEP SD
unique(data$speciesCat) #8 species of seagrass in the data set
unique(data$region) #5 regions in the data set

#==============================================================
# Plot the distribution of Net Ecosystem Production (Response Variable)
plot(density(data$NEP_mmolCm2d), 
     col = rgb(34, 124, 157, max = 255), pch = 20, cex = 2,
     xlab = expression(NEP~"("*mmol~m^-2~d^1*")"), main = "",
     cex.lab = 1.4, cex.axis = 1.2)
polygon(density(data$NEP_mmolCm2d), 
        col = rgb(34, 124, 157, max = 255), border = "white", lwd = 2)
lines(c(0,0), c(-1,1))
text(200,0.01, "Sink"); text(-70,0.01, "Source")

#==============================================================
# BAYESIAN MODEL - ESTIMATING MISSING VARIANCE
# As the mean increases, so does SD. We can use this relationship (which is the coefficient of variation) to estimate missing SDs in the data set

#Relationship between mean and sd for studies with that information
plot(data$NEP_sd ~ data$NEP_mmolCm2d)
summary(lm(data$NEP_sd ~ data$NEP_mmolCm2d))
abline(lm(data$NEP_sd ~ data$NEP_mmolCm2d))


set.seed(333)
sink("SDest_JAGS.txt") #creates a file with this model
cat("model{

    # Priors
    int ~ dnorm(0, .001) #intercept
    slope ~ dnorm(0, .001) #slope
    tau <- 1/(sigma * sigma) #precision
    sigma ~ dunif(0,10) #standard deviation

    # Model structure
    for(i in 1:R){
      Y[i] ~ dnorm(m[i],tau) 
      m[i] <- int + slope * X[i]
    }
    }", fill=TRUE)
sink()

#Input observations for the model
jags.data <- list(R = nrow(data), 
                  X = data$NEP_mmolCm2d, 
                  Y = data$NEP_sd)

#Initialize the 3 MCMC chains
inits <- function(){list(int = rnorm(1, 0, 5), 
                         slope = rnorm(1,0,5),
                         sigma = runif(1,0,10))}

#Info to run the model =================================================
params <- c("Y") #Parameter of interest is "y' in y = mx + b in the model
nc <- 3 #number of MCMC chains
n.adapt <-10000 #number of iterations to choose and optimize the sampler
n.burn <- 10000 #number of iterations for the burn-in period
n.iter <- 30000 #number of iterations to keep for the posterior
thin <- 10 #thin the number of iterations kept for the posterior

#initialize and optimize the chains
sd.model <- jags.model('SDest_JAGS.txt',
                       data = jags.data, 
                       inits = inits, 
                       n.chains = nc, 
                       n.adapt = n.adapt)
#burn in
update(sd.model, 
       n.burn) 
#iterations to keep for posterior
sd.model_samples <- coda.samples(sd.model, 
                                 params,
                                 n.iter = n.iter, 
                                 thin = thin)

#Create a new column with the absolute value of the estimated SD
data$NEP_sd_bayes <- abs(summary(sd.model_samples)$statistics[,1])

#=================================================================
# BAYESIAN MODEL - Estimate the mean NEP by seagrass species

data <- data[order(data$speciesCat),]
n.y = nrow(data)
n.spp = length(unique(data$speciesCat))
data = data %>%
  mutate(num_spp = case_when(speciesCat=="Zostera" ~ 1,
                             speciesCat=="Thalassia" ~ 2,
                             speciesCat=="Posidonia" ~ 3,
                             speciesCat=="Mixed" ~ 4))

obs_in = list(n = n.y, 
              n.group = n.spp,
              prior.min = -100, 
              prior.max = 300,
              y = as.double(data$NEP_mmolCm2d), 
              sigma = as.double(data$NEP_sd_bayes),
              group = as.numeric(data$num_spp))

# Initialize chains
inits = list(
  list(mu.alpha = -1.5, precision.alpha = 2, precision.delta = 1.5), #chain 1
  list(mu.alpha = 2, precision.alpha = 5, precision.delta = 3), #chain 2
  list(mu.alpha = 8, precision.alpha = 6, precision.delta = 5) #chain 3
)

n.adapt=50000 #number of iterations to choose and optimize the sampler
n.update=100000 #number of iterations for the burn-in period
n.iter=300000 #number of iterations to keep for the posterior

jm = jags.model("NEP_estimation_JAGS.R", 
                data = obs_in, 
                inits, 
                n.chains = length(inits), 
                n.adapt= n.adapt)

# Burn in the chain
update(jm, n.iter = n.update)
load.module("dic")
zc.spp = coda.samples(jm, 
                      variable.names = c("alpha"),
                      n.iter = n.iter,
                      thin = 10)

# Output - collect the samples from the chains for the posteroir distrubtions
spp.df <- as.data.frame(rbind(zc.spp[[1]], #the three MCMC chains
                              zc.spp[[2]], 
                              zc.spp[[3]]))
## NOTE: Change column headers if adapting for different analysis
colnames(spp.df) <- c("Zostera", "Thalassia", "Posidonia", "Mixed") 
spp.df <- as.data.frame(spp.df)

#===========================================
# Code checking Bayesian model

# 1) Check trace plots with 
plot(zc.spp)

# 2) Check Gelman and Rubin diagnostic: gelman.diag(coda.object)
# Approximate convergence is diagnosed when the upper limit is close to 1
gelman.diag(zc.spp)

# 3) Look at coda summaries and output tablular results
zc.summary<-summary(zc.spp)
spp.results<-cbind(zc.summary$stat[,1:2], zc.summary$quantile[,c(3,1,5)])
row.names(spp.results)<-c("Zostera", "Thalassia", "Posidonia", "Mixed")
spp.results

#=================================================
# Plot the posterior distributions
graphics.off()
plot(density(spp.df$Zostera), col = "gray30", lwd = 3, 
     ylim = c(0,0.07), main = "NEP Seagrass Species Estimate")
lines(density(spp.df$Thalassia), col = "turquoise3", lwd = 3, lty = 2)
lines(density(spp.df$Posidonia), col = "mediumseagreen", lwd = 3, lty = 3)
lines(density(spp.df$Mixed), col = "dodgerblue3", lwd = 3)
legend("topright", 
       legend = c("Zostera", "Thalassia", "Posidonia", "Mixed"),
       lty = c(1,2,3,1), 
       lwd = 3, 
       col = c("gray30", "turquoise3", "mediumseagreen", "dodgerblue3"))
