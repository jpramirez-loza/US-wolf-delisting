## Prepare dataset for JAGS model input ###

rm(list=ls())

# Start R in the repository root before running this script.
root_dir <- getwd()
if (!file.exists(file.path(root_dir, "all_timeseries.xlsx"))) {
  stop("Run data_prep.r from the repository root containing all_timeseries.xlsx.")
}
saveRDS(root_dir, "root_dir.rds")

setwd(root_dir)

dir.create(file.path(root_dir, "Figures"), showWarnings = FALSE)
dir.create(file.path(root_dir, "JAGS_Outputs"), showWarnings = FALSE)

library("readxl")

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

#"Metpop" is a vector with the biological population each state's population belongs to
# 1 is Northern Rockies, 2 is West Coast (Western OR and WA), 3 is Great Lakes

Metpop<-c(rep(1,5),rep(2,2),rep(3,3))

#Vector of observed population sizes at t = 1, for use in prior of true population size at first timestep
N1 <-c()
for (i in 1:length(series_names)){
  N1[i]<-pop.data[first.row[i],i]
}

#Make vectors of observation errors for time series of each state. 

ID_ObsErrors <- c(rep(NA,17),rep(1,11),rep(2,10),rep(NA,8))               #1 & 2 
MT_ObsErrors <- c(NA,rep(3,28),rep(4,16),NA)                                #3 & 4
WY_ObsErrors <- c(rep(NA,17),rep(5,28),NA)                                #5
WA_NRM_ObsErrors <- c(rep(NA,31),rep(6,14),NA)                            #6
OR_NRM_ObsErrors <- c(rep(NA,31),rep(7,14),NA)                            #7
WA_W_ObsErrors <- c(rep(NA,31),rep(8,14),NA)                              #8
OR_W_ObsErrors <- c(rep(NA,35),rep(9,10),NA)                              #9
WI_ObsErrors <- c(NA,rep(10,15),rep(11,9),rep(12,9),rep(13,8),rep(14,3),NA)   #10-14 #Wisconsin needs 5 different observation errors
MI_ObsErrors <- c(rep(NA,10),rep(15,18), rep(16,16),NA,NA)                   #15-16 #Michigan needs 2 different observation errors
MN_ObsErrors <- c(rep(NA,34),rep(17,10),NA,NA)                            #17

#Combine all ObsErrors vectors into a matrix 

ObsIndexMatrix <- cbind(ID_ObsErrors, MT_ObsErrors, WY_ObsErrors, WA_NRM_ObsErrors,OR_NRM_ObsErrors, 
                        WA_W_ObsErrors, OR_W_ObsErrors,
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
  nMetpop = max(Metpop),
  ObsIndexMatrix = ObsIndexMatrix,
  N1 = N1,
  x.pred = seq(0,1,0.01)
)

saveRDS(bugs.data, file = "bugs_data.rds")
