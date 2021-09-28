
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# flow disag #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 8 april 2020 

#### Author: David Woodson

#### CLEAR MEMORY
rm(list=ls())
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx20g"))
gc()
#### Clear the R console
cat("\f")

#packages
#library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(data.table)
library(ggrepel)
library(RColorBrewer)
library(readxl)
library(hydroGOF)
library(kknn)
library(plyr)
library(qpcR)
library(easyVerification)
library(verification)
library(qmap)
library(xlsx)
#library(kknn)
library(class)
library(anytime)
library(lubridate)
library(knnstdisagg)
library(xts)
library(XLConnect)
#library(openxlsx)

#disagg k-nn method annual forecasts to monthly, sub-basin
#read in obs monthly flows at MTOM sub-basins
data=readRDS('K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/esp/allFcstLocs_ESP_fcst.rds')

#read in annual flows from research method (knn resampling of observed flows and climate model temperatures)
#blind loocv ------------
fcst.projFlows.3yr = readRDS("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/fcst.sims_k=31_nsims=300_yrs=19822015_meanLength=3yr_fcstMode=kFold.Rds")
fcst.projFlows.5yr = readRDS("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/fcst.sims_k=27_nsims=300_yrs=19822013_meanLength=5yr_fcstMode=kFold.Rds")
fcst.projFlows.7yr = readRDS("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/fcst.sims_k=23_nsims=300_yrs=19822011_meanLength=7yr_fcstMode=kFold.Rds")
fcstMode = "kFold"
#non-blind loocv ------------
# fcst.projFlows.2yr = readRDS("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/fcst.sims_k=34_nsims=300_yrs=19822016_meanLength=2yr_fcstMode=loocv.Rds")
# fcst.projFlows.5yr = readRDS("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/fcst.sims_k=31_nsims=300_yrs=19822013_meanLength=5yr_fcstMode=loocv.Rds")
# fcst.projFlows.8yr = readRDS("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/fcst.sims_k=28_nsims=300_yrs=19822010_meanLength=8yr_fcstMode=loocv.Rds")
# fcst.projFlows.10yr = readRDS("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/fcst.sims_k=26_nsims=300_yrs=19822008_meanLength=10yr_fcstMode=loocv.Rds")
# fcstMode = "nonBlindLOOCV"

#select bootstrapped annual flows from knnBB flow projections 
#and subset to start in water year 1980 (so it matches ESP time domain)
annFlows.3yr = subset(fcst.projFlows.3yr$AnnualFlowBootstrap, startYear >= 1980 & covariate == "CESM-DPLE, RF")
annFlows.3yr$meanLength = 3
annFlows.5yr = subset(fcst.projFlows.5yr$AnnualFlowBootstrap, startYear >= 1980 & covariate == "CESM-DPLE, RF")
annFlows.5yr$meanLength = 5
annFlows.7yr = subset(fcst.projFlows.7yr$AnnualFlowBootstrap, startYear >= 1980 & covariate == "CESM-DPLE, RF")
annFlows.7yr$meanLength = 7


#annFlows = list(annFlows.8yr)
annFlows = list(annFlows.3yr, annFlows.5yr, annFlows.7yr)

nModels = length(annFlows)

#meanLengths = c(8)
meanLengths = c(3,5,7)


processDisagg = function(fcst.annFlows, model.i, meanLength, fcstMode){

  ###########################################################################################################
  #the translation of naturalized to unregulated flows and the knn space-time disagg requires that the fcst flows
  #be in a certain format 
  #so reformat data:
  ###############
  #identify simulation period and ensemble size:
  years = unique(fcst.annFlows$startYear)
  nFcsts = length(years)
  nSim = max(unique(fcst.annFlows$ensembleNo))
  #initialize list to store annual flows from each meanLength-year block
  af = vector("list", nFcsts)
  
  #loop thru each meanLength-year simulation block
  i = 1
  for(i in 1:nFcsts){
  
    #subset to current (i-th) meanLength-year block
    df.i = subset(fcst.annFlows, startYear == years[i])
    #initialize matrix to store reformatted flows
    mat.i = matrix(NA, nrow = nSim, ncol = meanLength)
    
    #loop through years 1-meanLength of the current meanLength-year forecast block and extract flows
    j = 1
    for(j in 1:meanLength){
      
      df.j = subset(df.i, fcstYear == j)
      
      mat.i[,j] = df.j$q.maf
      
    }
    
    colnames(mat.i) = unique(df.i$simYear)
    
    #put i-th year matrix into list:
    af[[i]] = mat.i
    
  }
  
  #reformating finished
  #######################
  
  #read in Lees Ferry naturalized annual flow
  ann.stats.ucrb = read.csv("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/txt_files/UCRB_annual_stats.txt", sep = " ", header = T)
  
  #select simulations using Nov-March minimum temperature as predictor (has best skill for year 2 forecasts)
  #af1 = ann.flows$DJF.Tmin
  #af2 = ann.flows$NovMar.Tmin
  #n = length(af1)
  
  ######################
  #process monthly data
  
  obs = data$locAll_obs
  #remove first and last few months to create water year format (range is now water years 1980-2019)
  obs = obs[-c(1:9,490:492),]
  #reset row names
  rownames(obs) = NULL
  #rearrange sites (columns)
  obs.mod = obs[,c(1,12,3,4,6,5,7,8,13,2,9,10,11)]
  #convert to xts object
  obs.monthly = xts(obs.mod[,-1], order.by = obs.mod$Timestep)
  
  dt = 1980:2019
  n = length(dt)
  props = vector("list", n)
  
  seq = seq(1, 12*n, 12)
  ann.agg.flows = matrix(NA, nrow = n, ncol = 2)
  ann.agg.flows[,1] = dt
  #ann.agg.flows[,2] = 1:n
  colnames(ann.agg.flows) = c("yrs", "LeesFerry")
  
  #index.flows.mon = data.frame(date = obs.mod[,1], flow = rowSums(obs.mod[,-1]))
  
  i = 1
  j = 1
  for(i in seq){
    
    #pull wy i data
    o.i = obs.mod[i:(i+11),]
    #calc annual flow
    #o.b = apply(o.i[,-1], 2, sum)
    o.a = sum(o.i[,2:13])
    
    ann.agg.flows[j,2] = o.a
    
    #p.i = o.i[,2:13]/o.a 
    
    #props[[j]] = p.i   
    
    j = j + 1
    
  }
  
  ###########################################
  #generate linear regression to transform annual naturalized at Lees Ferry to annual unregulated at LF
  ann.flows.lf = subset(ann.stats.ucrb[,c(1,3)], wyears >= min(ann.agg.flows[,1]))
  q.unreg = head(ann.agg.flows, -2)
  mod.data = data.frame(q.nat = ann.flows.lf$q_maf*1000, q.unreg = q.unreg[,2])
  
  model = lm(q.unreg ~ q.nat, data = mod.data)
  
  par(pty="s")
  title = paste0("Annual water year flow at Lees Ferry (1980-2017)\ny = ", round(model$coefficients[1], 0), " + ", 
                 round(model$coefficients[2], 2), "x")
  plot(q.unreg ~ q.nat, data = mod.data, xlab = "Annual naturalized flow (kaf)", ylab = "Annual unregulated flow (kaf)", font.lab = 2,
       main = title, cex.main = 0.9, cex.lab = 0.9)
  abline(a = model$coefficients[1], b = model$coefficients[2], col = "red")
  
  cor(mod.data$q.unreg, mod.data$q.nat)
  
  
  ###########################################
  ##loop thru processing of novel forecasts
  
  
  y = 1
  for(y in 1:nFcsts){
    
    af.y = af[[y]]
    
    #need to put disagg'd data into same format as required for MTOM input
  
    #apply monthly/spatial disaggregation to each block, each ensemble member
    #st1 = Powell Inflow unreg - GLDA3 -
    #st2 = Blue Mesa flow unreg - BMDC2 - 
    #st3 = CrystalInflow unreg - CLSC2 -
    #st4 = Fontenelle inflow - GBRW4 -
    #st5 = Flaming gorge inflow - GRNU1 - 
    #st6 = Morrow Point inflow - MPSC2 -
    #st7 = Navajo inflow - NVRN5 -
    #st8 = Taylor Park - TPIC2 -
    #st9 = Vallecito finlow - VCRC2 -
    #st10 = Yampa River inflow - YDLC2 -
    #st11 = Animas River inflow - DRGC2 - 
    #st12 = GainsCrystaltoGJ - GJLOC -
  
    #########################################################################
    #Create proportion vector for each sub-basin from unregulated flows
    
    #intialize matrix for each forecast location to store disagg'd data
    #rows are months and columns are ensemble members
    drgc2 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    bmdc2 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    clsc2 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    grnu1 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    gbrw4 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    mpsc2 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    nvrn5 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    gjloc = matrix(NA, nrow = meanLength*12, ncol = nSim)
    glda3 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    tpic2 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    vcrc2 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    ydlc2 = matrix(NA, nrow = meanLength*12, ncol = nSim)
    
    #ann.agg.flows = xts(ann.agg.flows, order.by = as.Date(ann.agg.flows[,1], format = '%Y'))
    
    #################################################
    #disaggregate forecasts to each sub-basin
    dat = af.y*10^3
    
    dat = dat*model$coefficients[2] + model$coefficients[1]
    #test = af.y[1,1]*10^3*model$coefficients[2] + model$coefficients[1]
    #test = af.y[2,2]*10^3*model$coefficients[2] + model$coefficients[1]
    # colnames(af.y)
    # ann.agg.flows.y = subset(ann.agg.flows)
    
    i = 1
    for(i in 1:nSim){
      
      q = as.matrix(data.frame(water_years = dt[y]:dt[y+(meanLength-1)], lees_ferry = dat[i,]))
      
      #names(q) = "q"
      
      disagg = knn_space_time_disagg(
        ann_flow = q,
        ann_index_flow = ann.agg.flows,
        mon_flow = obs.monthly,
        start_month = 10,
        nsim = 1,
        scale_sites = 1:12,
        k_weights = knn_params_default(nrow(ann.agg.flows))
      )
      
      q.p = knnst_get_disagg_data(disagg)
      
      drgc2[,i] = q.p[,1] 
      bmdc2[,i] = q.p[,2] 
      clsc2[,i] = q.p[,3] 
      grnu1[,i] = q.p[,4] 
      gbrw4[,i] = q.p[,5] 
      mpsc2[,i] = q.p[,6] 
      nvrn5[,i] = q.p[,7] 
      gjloc[,i] = q.p[,8] 
      glda3[,i] = q.p[,9] 
      tpic2[,i] = q.p[,10] 
      vcrc2[,i] = q.p[,11] 
      ydlc2[,i] = q.p[,12] 
        
      
      
    }
    
    #check
    #sum(drgc2[1:12,1], bmdc2[1:12,1], clsc2[1:12,1], grnu1[1:12,1], gbrw4[1:12,1], mpsc2[1:12,1], nvrn5[1:12,1], gjloc[1:12,1], glda3[1:12,1], tpic2[1:12,1], vcrc2[1:12,1], ydlc2[1:12,1])
    #sum(drgc2[1:12,2], bmdc2[1:12,2], clsc2[1:12,2], grnu1[1:12,2], gbrw4[1:12,2], mpsc2[1:12,2], nvrn5[1:12,2], gjloc[1:12,2], glda3[1:12,2], tpic2[1:12,2], vcrc2[1:12,2], ydlc2[1:12,2])
    #sum(drgc2[13:24,2], bmdc2[13:24,2], clsc2[13:24,2], grnu1[13:24,2], gbrw4[13:24,2], mpsc2[13:24,2], nvrn5[13:24,2], gjloc[13:24,2], glda3[13:24,2], tpic2[13:24,2], vcrc2[13:24,2], ydlc2[13:24,2])
    
    ######################################################################
    #write to excel 
    
    #initialize list of sub basins to store data (should have created such a list structure at the beginning for more concise code)
    #fcst.locs = vector("list", 12)
    locs = c("DRGC2", "BMDC2", "CLSC2", "GRNU1", "GBRW4", "MPSC2", "NVRN5", "GJLOC", "GLDA3", "TPIC2", "VCRC2", "YDLC2")
    #names(fcst.locs) = locs
    
    #loop thru basins to post-process
    dir = "K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/RiverWare/MTOM_v2/Inflow Forecasts/Testbed_Historical_Forecasts/novelHindcasts"
    sub.dir = model.i
    path = paste0(dir, '/', sub.dir)
    dir.create(path)
    setwd(path)
    
    i = 1
    for(i in 1:12){
    
      #fetch var based on list of names
      l.i = get(tolower(locs[i]))
      
      #replace neg w/ row mean (other other quantity, e.g. lower quartile)
      #m <- which(l.i < 0 , arr.ind=TRUE)
      #l.i[m] <- rowMeans(l.i, na.rm=TRUE)[m[,1]]
    
      #start year
      sy = as.numeric(colnames(af.y)[1])-1
      do = as.Date(paste0('01/10/', as.character(sy)), '%d/%m/%Y')
      #do = as.Date(paste0('10/1/', as.character(sy)), '%m/%d/%Y')
      Date0 = seq(from = do, by = '1 month', length.out = meanLength*12)
      #Date2 = gsub("0(\\d)", "\\1", Date0)    
    
      #Date = matrix(NA, nrow = length(Date2), ncol = 1)
      #for(j in 1:length(Date)){Date[j] = paste0(Date2[j], " 0:00:00")}
      #Date = mdy_hms(Date)
      #class(Date) = "POSIXct"
      #options(xlsx.datetime.format="m-d-yyyy H:mm:ss") # change date format
      
      Date = Date0
      
      #combine with dates and rename columns
      l.i = data.frame(Date, l.i)
      for(j in 1:nSim){colnames(l.i)[j+1] = paste("Trace", j, sep = "")}
      
      #as.Date(paste0('10/', as.character(sy)), '%d/%m')
      fname = paste0("RF_", nSim, "trace_fcst-", sy, "-10.xlsx")
      
      #options(xlsx.date.format="MM/dd/yyyy H:mm:ss")
      
      if(i == 1){
        write.xlsx(l.i, fname, sheetName = locs[i], row.names = F)
      } else {
        write.xlsx(l.i, fname, sheetName = locs[i], row.names = F, append = T)
        }
    
      cat(locs[i])
      cat("\n")
        
    }
  
    cat("\n")
    cat(paste0("Finished start year: ", sy))
    
    cat("\n")
    
  
  }

}

###########################################
#loop thru all different forecast models with function to pre-process, disagg, and output excel files

i = 1
for(i in 1:nModels){
  
  fcst.annFlows = annFlows[[i]]
  meanLength = fcst.annFlows$meanLength[1]

  model.i = paste0(gsub(", ", "_", fcst.annFlows$covariate[1]), "_meanLength_", fcst.annFlows$meanLength[1], "yr")
  
  cat("\n")
  cat(paste0("Starting: ", model.i))
  
  cat("\n")
  
  processDisagg(fcst.annFlows, model.i, meanLength, fcstMode)

}



##################################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Archived code
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################################################################


