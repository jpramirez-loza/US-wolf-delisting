############################################################
#                                                          #
# MODEL 7.3 Density dependent growth with state-specific
#          beta2's and harvest data. 
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

sheets_names<-excel_sheets("all_timeseries.xlsx")

all_data <- list()

for (sheet in sheets_names){
  all_data[[sheet]]<-read_excel("all_timeseries.xlsx", sheet = sheet,
                                col_names = TRUE,col_types =c(rep("numeric",2),rep("date",2),"numeric","date","date",rep("numeric", 3)))
}
#warnings are for NAs in columns classified as numeric


#Function to get the dates of listing and delisting for each time series
get.dates.delisted <- function(name){
  c(all_data[[name]]['Date_Delisted'][c(which(is.na(all_data[[name]]['Date_Delisted'])==FALSE)),])
}

get.dates.listed <- function(name){
  c(all_data[[name]]['Date_Listed'][c(which(is.na(all_data[[name]]['Date_Listed'])==FALSE)),])
}

#Create vectors with all dates of listing and delisting per state
listings <- list()
delistings <- list()

for (name in sheets_names){
  listings[[name]] <- rev(get.dates.listed(name))
  delistings[[name]] <- rev(get.dates.delisted(name))
}

#Function to create dataframe for each time series

create_timeseries_dataframe <- function(name){
  data.frame(all_data[[name]][c('Year','Popsize','MP_Total','Days_Listed','Days_Delisted','N_Killed')])
}

for(name in sheets_names){
  assign(name, create_timeseries_dataframe(name))
}

#This has created data frames named: 
# ID, MT, WY, WA_NRM, OR_NRM, WA_W, OR_W, WI, MI, MN
# WA_NRM and OR_NRM are the Northern Rocky Mountain sections of Washington and Oregon, respectively
# WA_W and OR_W are the Western sections of Washington and Oregon, respectively. 


#Names of the 8 time series
series_names <- c("ID","MT","WY","WA_NRM","OR_NRM",
                  #"WA_W","OR_W", # time series with exponential growth (West Washington and West Oregon) won't be included in this model 
                  "WI","MI","MN"
)

#Right now, row 1 in each data frame corresponds to the year 2023, 
#the code below rearranges the data frames so that row 1 of each time series
#is the year in which the series begins
rearrange_chronologically <- function(t.series){
  t.series<-t.series[rev(row.names(t.series)),]
  rownames(t.series) <- NULL
  return(t.series)
}

for (t.series_name in series_names){
  t.series <- get(t.series_name)
  t.series <- rearrange_chronologically(t.series)
  assign(t.series_name, t.series)
}

first_year <- function(name) {
  first.year <- min(name$Year)
  return(first.year)
}

last_year <- function(name) {
  last.year <- name[max(which(is.na(name['Popsize'])==0)),'Year']
  return(last.year)
}

#Make vectors with first and last years of each time series
series.start<-c()
series.end<-c()

for(i in 1:length(series_names)){
  timeseries<- get(series_names[i])
  series.start[i]<-first_year(timeseries)
  series.end[i]<-last_year(timeseries)
}

pop.data<-matrix(nrow=(length(c(1978:2023))), ncol= length(series_names))
days.listed.data <- matrix(nrow=(length(c(1978:2023))), ncol= length(series_names))
days.delisted.data <- matrix(nrow=(length(c(1978:2023))), ncol= length(series_names))
prop.listed.data <- matrix(nrow=(length(c(1978:2023))), ncol= length(series_names))
kill.data<-matrix(nrow=(length(c(1978:2023))), ncol= length(series_names))

for (j in 1:dim(pop.data)[2]){
  current_timeseries <- get(series_names[j])
  first <- series.start[j]
  last <- series.end[j]
  for (i in (first-1978+1):(last-1978+1)){
    pop.data[i,j]<-current_timeseries[(i-(first-1978+1)+1),2]
    days.listed.data[i,j]<-current_timeseries[(i-(first-1978+1)+1),4]
    days.delisted.data[i,j]<-current_timeseries[(i-(first-1978+1)+1),5]
    prop.listed.data[i,j]<-current_timeseries[(i-(first-1978+1)+1),4]/current_timeseries[(i-(first-1978+1)+1),3]
    kill.data[i,j]<-current_timeseries[(i-(first-1978+1)+1),6]
  }
}

# make vectors translating first and last years of time series to row numbers in matrix
first.row <- series.start-1978+1
last.row <- series.end-1978+1

#"Metpop" is a vector with the biological population each state's population belongs to
# 1 is Northern Rockies and 2 is Great Lakes 
Metpop<-c(rep(1,5),rep(2,3))

#Vector of observed population sizes at t = 1, for use in prior of true population size at first timestep
N1 <-c()
for (i in 1:length(series_names)){
  N1[i]<-pop.data[first.row[i],i]
}


ID_ObsErrors <- c(rep(NA,17),rep(1,11),rep(2,10),rep(NA,8))               #1 & 2 
MT_ObsErrors <- c(NA,rep(3,28),rep(4,16),NA)                                #3 & 4
WY_ObsErrors <- c(rep(NA,17),rep(5,28),NA)                                #5
WA_NRM_ObsErrors <- c(rep(NA,31),rep(6,14),NA)                            #6
OR_NRM_ObsErrors <- c(rep(NA,31),rep(7,14),NA)                            #7
# WA_W_ObsErrors <- c(rep(NA,31),rep(8,14),NA)                              #8
# OR_W_ObsErrors <- c(rep(NA,35),rep(9,10),NA)                              #9
WI_ObsErrors <- c(NA,rep(10,15),rep(11,9),rep(12,9),rep(13,8),rep(14,3),NA)   #10-14 #Wisconsin needs 5 different observation errors
MI_ObsErrors <- c(rep(NA,10),rep(15,18), rep(16,16),NA,NA)                   #15-16 #Michigan needs 2 different observation errors
MN_ObsErrors <- c(rep(NA,34),rep(17,10),NA,NA)                            #17

#Combine all ObsErrors vectors into a matrix 

ObsIndexMatrix <- cbind(ID_ObsErrors, MT_ObsErrors, WY_ObsErrors, WA_NRM_ObsErrors,OR_NRM_ObsErrors, 
                        # WA_W_ObsErrors, OR_W_ObsErrors,
                        WI_ObsErrors, MI_ObsErrors, MN_ObsErrors)


##############################################################################
# BUGS input
##############################################################################

bugs.data <- list(
  Nobs = pop.data,
  n.pops = length(series_names),
  first.row = first.row,
  last.row = last.row,
  Prop.delisted = 1-prop.listed.data,
  Metpop = Metpop,
  ObsIndexMatrix = ObsIndexMatrix,
  nObsErrors = max(ObsIndexMatrix, na.rm = TRUE),
  N1 = N1,
  x.pred = seq(0,1,0.01),
  N_killed = kill.data
)

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

    beta0[j] ~ dunif(0.5, 2)
    beta2[j] ~ dnorm(0, 0.000001)
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
  
      Nproc[t+1,j] <- log(max(1,lambda[t,j]*N[t,j]-N_killed[t+1,j]))
      N[t+1,j] ~ dlnorm(Nproc[t+1,j], tauProc[j])
    
      lambda[t,j] <- max(0, (beta0[j] + beta2[j]*N[t,j]))
      
      lambda.obs[t,j] <- N[t+1,j]/N[t,j]
    }

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
                  variable.names = c("beta0", "beta2","StateErrorObs"),
                  n.iter=samples,
                  thin=n.thin)


zj = jags.samples(jm, 
                  variable.names=c( "beta0","beta2","N",
                                    "lambda.obs",
                                    "lambda",
                                    "p.mean" ,"mean.y","mean.y.sim", "cv.y", "cv.y.sim", 
                                    "p.cv", "fit", "fit.sim", "p.fit"), 
                  n.iter=samples, 
                  thin=n.thin)

save(zc, file = paste0(outputs_dir,'/7.3Model7_coda.RData'))
save(zj, file = paste0(outputs_dir,'/7.3Model7_jags.RData'))


 summary(zc)

 # names to remove
 drop_vars <- c("sigma.vcov[1,2]", "sigma.vcov[1,3]", "sigma.vcov[2,3]", "sigma.vcov[3,3]")

 # keep all except those
 keep_vars <- setdiff(varnames(zc), drop_vars)

 # subset mcmc.list
 zc_filtered <- zc[, keep_vars, drop = FALSE]

 # now run Gelman diagnostics
 gelman.diag(zc_filtered, multivariate = TRUE)


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
#series_names = "ID","MT","WY","WA_NRM","OR_NRM","WA_W","OR_W","WI","MI","MN"
validation.plots("ID") #Validation plots for Idaho.