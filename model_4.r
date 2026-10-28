############################################################
#                                                          #
# MODEL 4. Density dependent growth with MVN's             #
#          vector of means shared by all states            #
#                                                          #
############################################################
rm(list=ls())

root_dir <- readRDS("root_dir.rds")
outputs_dir <- file.path(root_dir,"JAGS_Outputs") #directory for coda and jags objects output from model
figures_dir <- file.path(root_dir,"Figures")

setwd(root_dir)

library("readxl")
library("boot")
library("rjags")
library("coda")
library("scales")
library("lubridate")

bugs.data <- readRDS("bugs_data_DD.rds")

##############################################################################
# BUGS model
##############################################################################

sink("model.bug")
cat("

model
    {

################################
# Priors
################################

  for (j in 1:n.pops){
    tauProc[j] ~ dgamma(0.001, 0.001) #Prior for precision of true population process
    N[first.row[j],j] ~ dunif(0.5*N1[j],2*N1[j]) #The prior for the true population size at timestep 1, from half to twice the observed population
    N.sim[first.row[j],j] ~ dunif(0.5*N1[j],2*N1[j]) # For posterior predictive checks

    beta0[j] <- eps[1,j]
    beta1[j] <- eps[2,j]
    beta2[j] <- eps[3,j]
    
    eps[1:3,j] ~ dmnorm.vcov(mu[],sigma.vcov[,])
  }

    mu[1] ~ dunif(0.5,2)
    mu[2] ~ dnorm(0, 0.000001)
    mu[3] ~ dnorm(0, 0.000001)
    
    # Reparameterize sigma.vcov in terms of SD and rho
    for (i in 1:3){
      sigma.vcov[i,i] <- sigma[i] * sigma[i]
    }

    for (i in 1:2){
      for (j in (i+1):3){
        sigma.vcov[i,j] <- sigma[i] * sigma[j] * rho[i+j-2]
        sigma.vcov[j,i] <- sigma.vcov[i,j]
      } #j
    } #i

    # Specify priors for SD and rho
    for (i in 1:3){
      tau_betas[i] ~ dgamma(0.001, 0.001)
      sigma.squared[i] <- 1/tau_betas[i]
      sigma[i] <- sigma.squared[i]^0.5
      rho[i] ~ dunif(-1,1)
    }

    
    StateErrorObs[1] ~ dunif(0, 0.25) # Idaho observation error 1
    StateErrorObs[2] ~ dunif(0, 0.25) # Idaho observation error 2
    StateErrorObs[3] ~ dunif(0, 1) # Montana observation error 1
    StateErrorObs[4] ~ dunif(0, 1) # Montana observation error 2
    StateErrorObs[5] ~ dunif(0, 0.25) # Wyoming observation error
    StateErrorObs[6] ~ dunif(0, 0.25) # Washington Northern Rocky Mountains observation error
    StateErrorObs[7] ~ dunif(0, 0.25) # Oregon Northern Rocky Mountains observation error

    StateErrorObs[10] ~ dunif(0, 0.25) # Wisconsin observation error 1
    StateErrorObs[11] ~ dunif(0, 0.25) # Wisconsin observation error 2
    StateErrorObs[12] ~ dunif(0, 0.25) # Wisconsin observation error 3
    StateErrorObs[13] ~ dunif(0, 0.25) # Wisconsin observation error 4
    StateErrorObs[14] ~ dunif(0, 0.25) # Wisconsin observation error 5
    StateErrorObs[15] ~ dunif(0, 0.25) # Michigan observation error 1
    StateErrorObs[16] ~ dunif(0, 0.25) # Michigan observation error 2
    StateErrorObs[17] ~ dunif(0, 0.25) # Minnesota observation error
    

################################
# Process model
################################
   for (j in 1:n.pops){ 
    for (t in (first.row[j]):(last.row[j]-1)){
  
      Nproc[t+1,j] <- log(max(1,lambda[t,j]*N[t,j]))
      N[t+1,j] ~ dlnorm(Nproc[t+1,j], tauProc[j])
    
      lambda[t,j] <- max(0, (beta0[j] + beta1[j]*Prop.delisted[t+1,j] + beta2[j]*N[t,j])) 
      
      lambda.obs[t,j] <- N[t+1,j]/N[t,j]
    }
    
    lambda.draws[j] <- beta0[j] + beta1[j]
   }

#################################
# Observation model
#################################

 for (j in 1:n.pops){
   for (t in (first.row[j]:last.row[j])){
   
   sigmaObs[t,j] <- StateErrorObs[ObsIndexMatrix[t,j]]*N[t,j]
   
   #Gamma-Poisson
    shapeObs[t,j] <- N[t,j]*N[t,j]/(sigmaObs[t,j]*sigmaObs[t,j])
    rateObs[t,j] <- N[t,j]/(sigmaObs[t,j]*sigmaObs[t,j])

    lambdaNobs[t,j] ~ dgamma(shapeObs[t,j], rateObs[t,j])
    Nobs[t,j] ~ dpois(lambdaNobs[t,j])

    #for model validation 
     Nobs.sim[t,j] ~ dpois(lambdaNobs[t,j])

   }
 }
 
#################################
#Model predicitons
#################################

#Mean prediction of annual growth rate
for (k in 1:length(x.pred)){
  output[k] <- mu[1] + mu[2]*x.pred[k]
}


######################################
# Model validation
######################################

  for (j in 1:n.pops){
    mean.y[j] <- mean(Nobs[c((first.row[j]):(last.row[j]-1)),j])
    mean.y.sim[j] <- mean(Nobs.sim[c((first.row[j]):(last.row[j]-1)),j])
    p.mean[j] <- step(mean.y.sim[j] - mean.y[j])

    cv.y[j] <- sd(Nobs[c((first.row[j]):(last.row[j]-1)),j])/mean(Nobs[c((first.row[j]):(last.row[j]-1)),j])
    cv.y.sim[j] <- sd(Nobs.sim[c((first.row[j]):(last.row[j]-1)),j])/max(0.001,mean(Nobs.sim[c((first.row[j]):(last.row[j]-1)),j]))
    p.cv[j] <- step(cv.y.sim[j]-cv.y[j])

    for (t in (first.row[j]):(last.row[j]-1)){
      sq[t,j] <- (Nobs[t,j]-lambdaNobs[t,j])^2     #inside parentheses is Nobs - its expectation
      sq.sim[t,j] <- (Nobs.sim[t,j]-lambdaNobs[t,j])^2
    }
    fit[j] <- sum(sq[c((first.row[j]):(last.row[j]-1)),j])
    fit.sim[j] <- sum(sq.sim[c((first.row[j]):(last.row[j]-1)),j])
    p.fit[j] <- step(fit.sim[j]-fit[j])

  }

}
    ", fill=TRUE)
sink()

#############################################################################
# MCMC
#############################################################################

model = "model.bug"

n.chains=6

jm = jags.model(model,
                data = bugs.data,
                n.adapt = 100000, 
                n.chains = n.chains)

burn.in = 100000

update(jm, n.iter=burn.in)

samples = 1000000
n.thin = 50

zc = coda.samples(jm, 
                  variable.names = c("mu", "beta0", "beta1","beta2","StateErrorObs","sigma.vcov","rho"),
                  n.iter=samples,
                  thin=n.thin)


zj = jags.samples(jm, 
                  variable.names=c( "beta0", "beta1","beta2","N","mu","rho",
                                    "lambda.obs","lambda.draws",
                                    "lambda", "output",
                                    "p.mean" ,"mean.y","mean.y.sim", "cv.y", "cv.y.sim", 
                                    "p.cv", "fit", "fit.sim", "p.fit"), 
                  n.iter=samples, 
                  thin=n.thin)

save(zc, file = paste0(outputs_dir,'/4.Model4_coda.RData'))
save(zj, file = paste0(outputs_dir,'/4.Model4_jags.RData'))

# names to remove
drop_vars <- c("sigma.vcov[1,2]","sigma.vcov[1,3]","sigma.vcov[2,3]")

# keep all except those
keep_vars <- setdiff(varnames(zc), drop_vars)

# subset mcmc.list
zc_filtered <- zc[, keep_vars, drop = FALSE]

# now run Gelman diagnostics
gelman.diag(zc_filtered, multivariate = TRUE)

# gelman.diag(zc)

############################################
# Model validation results
############################################

# Bayesian p-values

zj$p.mean
zj$p.cv
zj$p.fit

series_names <- c("ID","MT","WY","WA_NRM","OR_NRM",
                  # "WA_W","OR_W",
                  "WI","MI","MN")

#Make a function that rearranges the zj mc arrays into matrices
#where rows correspond to populations
zj.as.matrix <- function(metric){
  metric.matrix <- matrix(as.vector(zj[[metric]]), nrow = length(series_names), byrow=FALSE)
  return(metric.matrix)
}

list.of.metrics <- c("mean.y","mean.y.sim", "cv.y", "cv.y.sim",
                     "fit", "fit.sim")

validation.metrics <- list()
for (i in 1:length(list.of.metrics)){
  validation.metrics[[i]] <- zj.as.matrix(list.of.metrics[i])
}

validation.plots <- function(name){
  position <- which(series_names[]==name)
  par(mfrow=c(2,2))

  p.mean <- summary(zj$p.mean,mean)$stat[position]
  plot(validation.metrics[[1]][position,],validation.metrics[[2]][position,], xlab="mean data", ylab="mean simulated",cex=.05)
  mtext(paste("P =",round(as.numeric(p.mean),4)), side = 3)
  abline(0,1)

  p.cv <- summary(zj$p.cv,mean)$stat[position]
  plot(validation.metrics[[3]][position,],validation.metrics[[4]][position,], xlab="Coeff of variation real data", ylab="Coeff of variation simulated data",cex=.05)
  mtext(paste("P =",round(as.numeric(p.cv),4)), side = 3)
  abline(0,1)

  p.fit <- summary(zj$p.fit,mean)$stat[position]
  plot(validation.metrics[[5]][position,],validation.metrics[[6]][position,], xlab="Sum of squares real data", ylab="Sum of squares simulated data",cex=.05) #xlim=c(0,10), ylim = c(0,10))
  mtext(paste("P =",round(as.numeric(p.fit),4)), side = 3)
  abline(0,1)
}

#plot model fit by inserting element from series_names into validation.plots function.
#series_names = "ID","MT","WY","WA_NRM","OR_NRM","WI","MI","MN"
validation.plots("ID") #Validation plots for Idaho.


####################
# Figures
####################

#Create data frames
beta0 <- as.vector(zj$mu[1,,])
beta1 <- as.vector(zj$mu[2,,])
b0b1 <- as.vector(zj$mu[1,,]) + as.vector(zj$mu[2,,])

d0 <- density(beta0)
x_fill0 <- d0$x[d0$x <= 1]
y_fill0 <- d0$y[d0$x <= 1]

d1 <- density(b0b1)
x_fill1 <- d1$x[d1$x <= 1]
y_fill1 <- d1$y[d1$x <= 1]

DL.p <- round(ecdf(b0b1)(1),3)
L.p <- round(ecdf(beta0)(1),4)
DLoL.R <- round((ecdf(b0b1)(1))/(ecdf(beta0)(1)),3)

d1 <- density(beta1)
x_fill1 <- d1$x[d1$x <= 0]
y_fill1 <- d1$y[d1$x <= 0]


pdf(file = paste0(figures_dir,"/M4_Beta1.pdf"))

plot(density(beta1), 
     ylim= c(0,max(density(beta1)$y)),
     main = "",xlab = "Effect of delisting on growth rate", ylab = "",
     xlim = c(-0.5,0.5), lty = 1, lwd = 2, col = "grey",
     cex.axis = 1.5, 
     cex.lab = 1.8)


polygon(
  c(x_fill1, rev(x_fill1)),
  c(rep(0, length(x_fill1)), rev(y_fill1)),
  density = 10,           # number of lines per inch
  angle = -45,             # diagonal direction
  col = "black",
  border = NA
)

abline(v = 0, lty = 2)
abline(h = 0)

title(ylab="Probability density", mgp = c(2.5, 1.0, 0),cex.lab = 1.8)

dev.off()



