############################################################
#                                                          #
# MODEL 3. Exponential growth with independent estimate    #
#          of beta0 and beta1 per state                    #
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
library("dplyr")

bugs.data <- readRDS("bugs_data.rds")

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
    beta1[j] ~ dnorm(0, 0.000001)
    beta0[j] ~ dunif(0.5, 2)
  }

    
    #Priors for observation errors
    StateErrorObs[1] ~ dunif(0, 0.25) # Idaho observation error 1
    StateErrorObs[2] ~ dunif(0, 0.25) # Idaho observation error 2
    StateErrorObs[3] ~ dunif(0, 1) # Montana observation error 1
    StateErrorObs[4] ~ dunif(0, 1) # Montana observation error 2
    StateErrorObs[5] ~ dunif(0, 0.25) # Wyoming observation error
    StateErrorObs[6] ~ dunif(0, 0.25) # Washington Northern Rocky Mountains observation error
    StateErrorObs[7] ~ dunif(0, 0.25) # Oregon Northern Rocky Mountains observation error
    StateErrorObs[8] ~ dunif(0, 0.25) # West Washington observation error
    StateErrorObs[9] ~ dunif(0, 0.25) # West Oregon observation error
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
    
      lambda[t,j] <- max(0, (beta0[j] + beta1[j]*Prop.delisted[t+1,j]))
      
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

#State-specific predictions of annual growth rates
  for (j in 1:n.pops){
    for (k in 1:length(x.pred)){
      output.states[j,k] <- beta0[j] + beta1[j]*x.pred[k]
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
                  variable.names = c("beta0", "beta1","StateErrorObs"),
                  n.iter=samples,
                  thin=n.thin)


zj = jags.samples(jm, 
                  variable.names=c( "beta0", "beta1","N",
                                    "lambda.obs","lambda.draws",
                                    "lambda", "output.states",
                                    "p.mean" ,"mean.y","mean.y.sim", "cv.y", "cv.y.sim", 
                                    "p.cv", "fit", "fit.sim", "p.fit"), 
                  n.iter=samples, 
                  thin=n.thin)

save(zc, file = paste0(outputs_dir,'/3.Model3_coda.RData'))
save(zj, file = paste0(outputs_dir,'/3.Model3_jags.RData'))

summary(zc)
gelman.diag(zc)

############################################
# Model validation results
############################################

# Bayesian p-values

zj$p.mean
zj$p.cv
zj$p.fit

series_names <- c("ID","MT","WY","WA_NRM","OR_NRM",
                  "WA_W","OR_W",
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

####################
# Time series figures (Column 1 of Figure 4)
####################

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

#Names of the 10 time series
series_names <- c("ID","MT","WY","WA_NRM","OR_NRM",
                  "WA_W","OR_W",
                  "WI","MI","MN")

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

for (j in 1:dim(pop.data)[2]){
  current_timeseries <- get(series_names[j])
  first <- series.start[j]
  last <- series.end[j]
  for (i in (first-1978+1):(last-1978+1)){
    pop.data[i,j]<-current_timeseries[(i-(first-1978+1)+1),2]
    days.listed.data[i,j]<-current_timeseries[(i-(first-1978+1)+1),4]
    days.delisted.data[i,j]<-current_timeseries[(i-(first-1978+1)+1),5]
    prop.listed.data[i,j]<-current_timeseries[(i-(first-1978+1)+1),4]/current_timeseries[(i-(first-1978+1)+1),3]
  }
}

# make vectors translating first and last years of time series to row numbers in matrix
first.row <- series.start-1978+1
last.row <- series.end-1978+1


Pop.pred<-summary(zj$N, quantile, c(.025,.5,.975),na.rm=TRUE)$stat
Pop.pred.vector <- as.vector(Pop.pred)
#45 is the maximum number of timesteps for any population, it can be checked as the number of columns
#of each element in Pop.pred

max.timesteps <-45

Percentiles<-matrix(Pop.pred.vector, nrow = 3, ncol=length(series_names)*max.timesteps, byrow=FALSE)

divided_matrices <- list()

# Loop to divide the original matrix into one matrix per population.
for (i in 1:length(series_names)) {
  start_col <- (i - 1) * max.timesteps + 1   # Calculate starting column index for population size
  
  end_col <- i * max.timesteps                # Calculate ending column index
  
  # Extract the subset of columns from the original matrix
  divided_matrix <- Percentiles[, start_col:end_col]
  
  # Store the subset matrix in the list
  divided_matrices[[i]] <- divided_matrix
}

plot.titles <- c("Idaho", "Montana", "Wyoming","Washington NRM", "Oregon NRM",
                 "Washington West", "Oregon West",
                 "Wisconsin","Michigan","Minnesota")

plot.delistings <- function(name){
  position <- which(series_names[]==name)
  state.data <- get(series_names[position])
  pred<-divided_matrices[[position]]
  x=seq((series.start[position]+1),(series.end[position]+1),1)
  y=pred[,c(first.row[position]:last.row[position])]
  
  plot(x,y[2,], col = "white", type = "l", main = plot.titles[position], 
       xlab="Year", ylab= "", 
       cex.axis = 1.8, cex.lab = 2, tck=-0.03, bty="l", 
       xlim = c(series.start[position]+1,series.end[position]+1),
       ylim = c(0, (max(y)+0.1*(max(y)))) )
  
  series.listings<- rev(as.Date(listings[[name]]$Date_Listed[]))
  series.delistings <-rev(as.Date(delistings[[name]]$Date_Delisted[]))
  listings.years <- year(series.listings[])
  delistings.years <- year(series.delistings[])
  listing.day.within.year <- yday(series.listings[])
  delisting.day.within.year <-yday(series.delistings[])
  continuous.listings <- year(series.listings[]) + yday(series.listings[])/365
  continuous.delistings <- year(series.delistings[]) + yday(series.delistings[])/365
  
  n.listings <- length(continuous.listings)
  n.delistings <- length(continuous.delistings)
  
  pred<-divided_matrices[[position]]
  y=pred[,c(first.row[position]:last.row[position])]
  
  
  if (n.delistings > n.listings){
    for(i in 1:(n.delistings-1)){
      polygon(c(seq(continuous.delistings[i],continuous.listings[i],0.01),rev(seq(continuous.delistings[i],continuous.listings[i],0.01))),
              c(rep(max(y)+600,length(seq(continuous.delistings[i],continuous.listings[i],0.01))),rep(-200,length(seq(continuous.delistings[i],continuous.listings[i],0.01)))),
              col= alpha("grey",alpha = 0.4), border = NA)
    }
    
    polygon(c(seq(continuous.delistings[n.delistings],series.end[position]+1,0.01),rev(seq(continuous.delistings[n.delistings],series.end[position]+1,0.01))),
            c(rep(max(y)+600,length(seq(continuous.delistings[n.delistings],series.end[position]+1,0.01))),rep(-200,length(seq(continuous.delistings[n.delistings],series.end[position]+1,0.01)))),
            col= alpha("grey",alpha = 0.4), border = NA)
  } else {
    for(i in 1:n.delistings){
      polygon(c(seq(continuous.delistings[i],continuous.listings[i],0.01),rev(seq(continuous.delistings[i],continuous.listings[i],0.01))),
              c(rep(max(y)+600,length(seq(continuous.delistings[i],continuous.listings[i],0.01))),rep(-200,length(seq(continuous.delistings[i],continuous.listings[i],0.01)))),
              col= alpha("grey",alpha = 0.4), border = NA)
    }
  }
  
  polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), col=alpha("lightblue", alpha=1),border = NA)
  lines(x,y[2,], col = "black")
  points((state.data$Year+1), state.data$Popsize, pch = 16)
  
  title(ylab="Number of wolves", mgp = c(2.65, 1.0, 0),cex.lab = 2)
  
}

for (i in 1:length(series_names)){
  pdf(file = paste0(figures_dir,"/M3_tseries_",series_names[i], ".pdf"))
  plot.delistings(series_names[i])
  dev.off()
}


####################
# Model fit figures (Column 2 of Figure 4)
####################

states.names <- c("Idaho", "Montana", "Wyoming", "Washington NRM", "Oregon NRM",
                  "Washington West", "Oregon West",
                  "Wisconsin", "Michigan","Minnesota")

pred_samples <-  list()

for(i in 1:10){
  pred_samples[[i]] <- as.vector(zj$output.states[i,,,])
}

x.pred <- seq(0,1,0.01) # values of proportion of the year delisted for which lambda was predicted

#Probabilities of lambda <1 after one year of delisting
Prob <- round(ecdf(zj$output.states[3,101,,])(1), 3)

# Initialize outer list
lambdas <- vector("list", 10)

for (i in 1:10) {
  # Initialize inner list for each i
  lambdas[[i]] <- vector("list", 44)
  
  for (j in 1:44) {
    lambdas[[i]][[j]] <- quantile(
      zj$lambda.obs[j, i, , ],
      probs = c(0.025, 0.5, 0.975),
      na.rm = TRUE
    )
  }
}

predictions <- vector("list", 10)  # list of length 10

for (j in 1:10) {
  predictions[[j]] <- vector("list", 101)  # sublist for each j
  
  for (i in 1:101) {
    predictions[[j]][[i]] <- quantile(
      zj$output.states[j, i, , ],
      probs = c(0.025, 0.5, 0.975),
      na.rm = TRUE
    )
  }
}

df_pred <- list()

for(i in 1:10){
  df_pred[[i]] <- data.frame(
    parameter = rep("Pred", length(x.pred)),
    PYDL = x.pred,
    lower_CrI = sapply(predictions[[i]], `[`, 1),
    median     = sapply(predictions[[i]], `[`, 2),
    upper_CrI  = sapply(predictions[[i]], `[`, 3)
  )
}

#Data frames for "observed" lambdas (computed from model)

first.row <- bugs.data$first.row
last.row <- bugs.data$last.row

lambdas.list <- list()
for (i in 1:bugs.data$n.pops){
  lambdas.list[[i]] <- data.frame(
    parameter = rep(series_names[i],length(seq(first.row[i],(last.row[i]-1),1))),
    PYDL = bugs.data$Prop.delisted[(first.row[i]+1):last.row[i],i],
    lower_CrI = sapply(lambdas[[i]][(first.row[i]):(last.row[i]-1)], `[`, 1),
    median = sapply(lambdas[[i]][(first.row[i]):(last.row[i]-1)], `[`, 2),
    upper_CrI = sapply(lambdas[[i]][(first.row[i]):(last.row[i]-1)], `[`, 3)
  )
}

all_points <- bind_rows(lambdas.list)

state_points <- list()

for (i in 1:10){
  state_points[[i]] <- all_points[which(all_points[,'parameter']==series_names[i]),]
}


for (i in 1:10){
  pdf(file = paste0(figures_dir,"/M3_fit_",states.names[i],".pdf"))
  plot(df_pred[[i]]$PYDL,df_pred[[i]]$median, 
       xlim = c(0,1), ylab = "", ylim = c(0.5,2), type = "l",
       lwd = 2, xlab = "Proportion of the year delisted",
       cex.axis = 1.8, 
       cex.lab = 2.0, 
       main = states.names[[i]])
  
  title(ylab="Growth rate", mgp = c(2.65, 1.0, 0),cex.lab = 2.0)
  
  polygon(
    c(df_pred[[i]]$PYDL, rev(df_pred[[i]]$PYDL)),
    c(df_pred[[i]]$lower_CrI, rev(df_pred[[i]]$upper_CrI)),
    col = alpha("grey", 0.5),
    border = NA
  )
  
  points(state_points[[i]]$PYDL, state_points[[i]]$median, pch = 16)
  lines(df_pred[[i]]$PYDL,df_pred[[i]]$median)
  
  abline(h = 1.0, lty = 2)
  
  dev.off()
}

####################
# Predicted growth rates figures (Column 3 of Figure 4)
####################

beta0 <- list()
beta1 <- list()
b0b1 <- list()


for (i in 1:10){
  beta0[[i]]  <- as.vector(zj$beta0[i,,])
  beta1[[i]]  <- as.vector(zj$beta1[i,,])
  b0b1[[i]] <- as.vector(zj$beta0[i,,]) + as.vector(zj$beta1[i,,])
}

round(ecdf(beta1[[10]])(0),3)

b1.low <- c(-0.8,-0.4,-0.5,-1.1,-1.1,-1,-2,-0.4,-0.8,-0.5)
b1.hi <-c(0.1,0.1,0.2,1.0,1,2.5,1,0.2,0.1,0.5)

for (i in 1:10){
  pdf(file = paste0(figures_dir,"/M3_beta1_",states.names[i],".pdf"))
  
  d1 <- density(beta1[[i]])
  x_fill1 <- d1$x[d1$x <= 0]
  y_fill1 <- d1$y[d1$x <= 0]
  
  plot(density(beta1[[i]]), 
       ylim= c(0,max(density(beta1[[i]])$y)),
       main = "", xlab = "Effect of delisting on growth rate",
       ylab = "", yaxt = "n", lty = 1, lwd = 2, col = "grey",
       xlim = c(b1.low[i],b1.hi[i]),
       cex.axis = 2.3, 
       cex.lab = 2.7)
  
  
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
  
  
  dev.off()
}

for (i in 1:10){
  pdf(file = paste0(figures_dir,"/M3_predGR_",states.names[i],".pdf"))
  
  d0 <- density(beta0[[i]])
  x_fill0 <- d0$x[d0$x <= 1]
  y_fill0 <- d0$y[d0$x <= 1]
  
  d1 <- density(b0b1[[i]])
  x_fill1 <- d1$x[d1$x <= 1]
  y_fill1 <- d1$y[d1$x <= 1]
  
  plot(density(beta0[[i]]),
       ylim= c(0,max(max(density(beta0[[i]])$y,max(density(b0b1[[i]])$y)))),
       main = states.names[i], xlab = "Predicted growth rate", ylab = "",
       xlim = c(0,2), lty = 1, lwd = 2, col = "grey",
       cex.axis = 1.8,
       cex.lab = 2.0)
  
  polygon(
    c(x_fill0, rev(x_fill0)),
    c(rep(0, length(x_fill0)), rev(y_fill0)),
    col = alpha("grey", 0.9),
    border = NA
  )
  
  polygon(
    c(x_fill1, rev(x_fill1)),
    c(rep(0, length(x_fill1)), rev(y_fill1)),
    density = 10,           # number of lines per inch
    angle = 45,             # diagonal direction
    col = "black",
    border = NA
  )
  
  lines(density(b0b1[[i]]), col = "black", lwd = 2)
  abline(v = 1.0, lty = 2)
  abline(h = 0)
  
  title(ylab="Probability density", mgp = c(2.65, 1.0, 0),cex.lab = 2.0)
  
  legend("topright", legend = c("Listed", "Delisted"), lty = c(1,1),
         lwd = c(2,2), col = c("grey","black"), bty = "n", cex = c(1.5,1.5))
  
  dev.off()
}

DL.p <- c()
L.p <- c()
DLoL.R <- c()

for (i in 1:10){
  DL.p[i] <- round(ecdf(b0b1[[i]])(1),3)
  L.p[i] <- round(ecdf(beta0[[i]])(1),3)
  DLoL.R[i] <- round((ecdf(b0b1[[i]])(1))/(ecdf(beta0[[i]])(1)),3)
}

