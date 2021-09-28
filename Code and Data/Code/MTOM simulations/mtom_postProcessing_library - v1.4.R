#packages
library(dplyr)
#library(tidyr)
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
library(hydroGOF)
library(ggplot2)
library(reshape)
library(ggbeeswarm)

options(warnings = -1)
septOnly = T

### Define skill functions
#############################################
#calc NSE on a time series, using climatology
nse = function(mod, obs, clim){
  
  nse = 1 - sum((mod-obs)^2)/sum((obs-clim)^2)
  
  return(nse)
  
}

########################################################
#calc NSE for each ensemble member or for ensemble mean
#
# mod=fcst.mat
# obs=hist.sub
# clim=esp.mat
ens.nse = function(mod, obs, clim, meanLength.N, septOnly, septMonths) {
  
  #list to store annual and 5-year NSEs for function output
  nse.list = list()
  
  #calc forecast ensemble mean
  fcst.EnsMean = apply(mod, 1, mean) 
  #calc clim ensemble mean
  clim.EnsMean = apply(clim, 1, mean) 
  
  #5-year nse---------------------------------------------
  #calc for ens mean----------
  nse.mat = nse(fcst.EnsMean, obs, clim.EnsMean)

  #annual nse---------------------------------------------
  #nse.ann = matrix(NA, nrow = 1, ncol = 5)
  nse.ann = matrix(NA, nrow = 1, ncol = meanLength.N)
  
  k = 1
  for(i in 1:ncol(nse.ann)){
    nse.ann[i] = nse(fcst.EnsMean[k:(k+11)], obs[k:(k+11)], clim.EnsMean[k:(k+11)])
    k = k + 12
  }
    
  #average monthly nse for years 2-4 ----------------------------------
  nse.mean.24 = mean(nse.ann[2:(meanLength.N-1)]) 
  
  #calc NSE simultaneously on entire ensemble for each EOWY -----------------
  if(septOnly == T){
    nseEns = matrix(NA, ncol = 1, nrow = length(meanLength.N))
    
    #subset matrices to only sept months
    fcst.sept = mod[septMonths,]
    clim.sept = clim[septMonths,]
    obs.sept = obs[septMonths]
    
    i = 1
    for(i in 1:meanLength.N){
      
      obs.i = obs.sept[i]
      nseEns[i] = 1 - sum((obs.i-fcst.sept[i,])^2)/sum((obs.i-clim.sept[i,])^2)
        
    }
    
  }
  
  #output
  nse.list$nse.ann = nse.ann
  nse.list$nse = nse.mat
  nse.list$nse.mean.24 = nse.mean.24
  nse.list$nseEns = nseEns

  return(nse.list)
  
}

##################################
#calc RMSE for each ensemble member
#option to calculate on September 31 pool elevation only
#
mod=fcst.mat
obs=hist.sub
clim=esp.mat
crpss.on=T
ens.rmse = function(mod, obs, septOnly, dates, meanLength.N, crpssOn, clim) {

  #list to store annual and N-year RMSEs for function output
  rmse.list = list()
  
  #5-year rmse-----------------------------------------
  rmse = matrix(NA, nrow = ncol(mod), ncol = 1)
  #annual rmse
  # if(meanLength > 5){
  #   ncol = 5 #since ESP duration is 5-years, subset novel forecast lengths longer than 5-years
  # } else{
  #   ncol = meanLength
  # }
  rmse.ann = matrix(NA, nrow = ncol(mod), ncol = meanLength.N)
  
  i = 1
  for(i in 1:ncol(mod)){
    k = 1
    for(j in 1:ncol(rmse.ann)){
      rmse.ann[i,j] = rmse(mod[k:(k+11),i], obs[k:(k+11)])
      k = k + 12
    }
    rmse[i] = rmse(mod[,i], obs)
  }
  
  #average monthly nse for years 2-4-------------------
  if(meanLength.N >=4){
     rmse.mean.24 = apply(rmse.ann[,2:(meanLength.N-1)], 1, mean)  
  } else{
     rmse.mean.24 = mean(rmse.ann[,2:(meanLength.N-1)])  
  }
  
  #ens mean rmse----------------------------------------
  ens.mean = apply(mod, 1, mean)
  ens.mean.rmse = rmse(ens.mean, obs)
  
  #calc RMSE on dec pool elevations only (if desired)--------
  if(septOnly == T){
    
    tmp.mod = cbind(as.character(dates), mod)
    tmp.mod = subset(tmp.mod, format.Date(dates, "%m")=="09")

    mod = apply(tmp.mod[,-1], 2, as.numeric)
    
    tmp.obs = cbind(as.character(dates), obs)
    tmp.obs = subset(tmp.obs, format.Date(dates, "%m")=="09")

    obs = as.numeric(tmp.obs[,-1])
    
    #5-year rmse-----------------------------------------
    rmse.dec = matrix(NA, nrow = ncol(mod), ncol = 1)
    i = 1
    for(i in 1:ncol(mod)){
      rmse.dec[i] = rmse(mod[,i], obs)
    }
    
    #ensemble RMSE (rmse calculated separately by year)
    #now calc Dec rmse by year but calc RMSE on separate ensemble members grouped by year instead of time series
    rmse.eoy.Ens = matrix(NA, nrow = nrow(mod), ncol = 1)
    i = 1
    for(i in 1:nrow(mod)){
      rmse.eoy.Ens[i] = rmse(mod[i,], rep(obs[i], length(mod[i,])))
    }
  
    #eoy crpss 
    if(crpssOn == T){
      
      tmp.clim = cbind(as.character(dates), clim)
      tmp.clim = subset(tmp.clim, format.Date(dates, "%m")=="09")
      
      clim = apply(tmp.clim[,-1], 2, as.numeric)
      
      #calc RPSS on forecast and climatology
      crps.f = EnsCrps(ens = mod, obs = obs)
      crps.c = EnsCrps(ens = clim, obs = obs)
      
      #calc CRPSS and store
      crpss.eoy = 1 - crps.f/crps.c
      
    }
   
  }
  
  #function output---------------------------------
  rmse.list$rmse.ann = rmse.ann
  rmse.list$rmse = rmse
  rmse.list$rmse.mean.24 = rmse.mean.24
  rmse.list$rmse.ens.mean = ens.mean.rmse
  if(septOnly == T){
    rmse.list$rmse.dec = rmse.dec
    rmse.list$rmse.eoy.byEns = rmse.eoy.Ens
    if(crpssOn == T){
      rmse.list$crpss.eoy = crpss.eoy
    }
  }
  
  return(rmse.list)
  
}



### Loop function
###################################################################################################################################

#processing, plotting and skill calcs for selected reservoirs

###################################################################################################################################

# years = c(1981:2012)
# reservoir = "Powell"
# fcstType = "nonBlindLOOCV"
# meanLength = 5

reservoir.analysis = function(years, reservoir, fcstType, meanLength, subdir, nseOn){
  
  #initialize storage lists/dataframes for later aggregate plots
  if(nseOn == T){
    nse.5year.list = vector("list", length(years)) 
    nse.ann.df = NULL
    nse.24.df = NULL
    nse.ensMean.df = NULL
  }
  
  rmse.5year.list = vector("list", length(years))
  rmse.ann.df = NULL
  rmse.24.df = NULL
  rmse.ensMean.df = NULL
  rmse.eoy.df = NULL
  rmse.eoy.byEns.df = NULL
  
  crpss.eoy.df = NULL
  
  dir = "K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/RiverWare/MTOM_v2/Scenario/"
  
  if(meanLength<5){
    meanLength.N = meanLength    
  } else{
    meanLength.N = 5 #set analysis and plotting for mean lengths > 5 years to 5 years (since ESP and historic scenarios are only 5 years long)
  }

  #Read and calculate climatology
  hist.clim = NULL
  for(year in years){
    
    #read in historical simulation - 1 trace
    file = paste0(dir, "Historical/Historical_DMI,Historical MRM,MTOM Model,Rules,HistoricalStreamflow-", year, "-10/", reservoir, "Output.csv")  
    hist.c = read.csv(file)
    #subset to pool elevation, units are ft
    hist.sub.c = subset(hist.c$Slot.Value, hist.c$Object.Slot == paste0(reservoir, ".Pool Elevation"))
    
    hist.clim = cbind(hist.clim, hist.sub.c)
  }
  clim = mean(hist.clim)
  
  y = 1
  year = 1981
  for(year in years){
    
    #choose epoch - e100 test
    #######################################################################
    #read in ESP simulations - 29 traces
    file = paste0(dir, "ESP/ESP29Trace_DMI,MRM,MTOM Model,Rules,ESP_29Trace_fcst-", year, "-10/ReservoirOutput.csv")  
    esp = read.csv(file)
    
    #subset to pool elevation, units are ft
    esp.sub = subset(esp, esp$Object.Slot == paste0(reservoir, ".Pool Elevation"))
    #rearrange so trace numbers are ascending
    esp.sub = arrange(esp.sub, Trace.Number)
    
    #create matrix to store traces and more easily plot
    n.trace.esp = length(unique(esp.sub$Trace.Number))
    dt = length(unique(esp.sub$Timestep))
    esp.mat = matrix(NA, nrow = dt, ncol = n.trace.esp)
    
    i = 1
    for(i in 1:n.trace.esp){
      
      esp.i = subset(esp.sub$Slot.Value, esp.sub$Trace.Number == i)
      
      esp.mat[,i] = esp.i
      
    }
    
    esp.mat = esp.mat[1:(meanLength.N*12),]
    
    #######################################################################
    #read in fcst BB simulations - 30 traces
    # if(fcstType == "blindLOOCV"){
    #   subdir = paste0("RF300_", fcstType, "/RF300_", meanLength, "yr/RF300_", meanLength, "yr_DMI,RF300_", meanLength, "yr MRM ", fcstType, ",MTOM Model,Rules,RF_300trace_fcst-")
    # }
    # if(fcstType == "nonBlindLOOCV"){
    #   subdir = paste0("RF300_", fcstType, "/RF300_", meanLength, "yr/RF300_", meanLength, "yr_DMI_", fcstType, ",RF300_", meanLength, "yr MRM ", fcstType, ",MTOM Model,Rules,RF_300trace_fcst-")
    # }
    file = paste0(dir, subdir, year, "-10/ReservoirOutput.csv")  
    fcst = read.csv(file)
    
    #subset to pool elevation, units are ft
    fcst.sub = subset(fcst, fcst$Object.Slot == paste0(reservoir, ".Pool Elevation"))
    #rearrange so trace numbers are ascending
    fcst.sub = arrange(fcst.sub, Trace.Number)
    dt.full = length(unique(fcst.sub$Timestep))
    if(meanLength < 5){dt = length(unique(fcst.sub$Timestep))} #if meanLength less than 5 yrs, subset fcst dt (else it's the 5-year dt of ESP)
    
    #create matrix to store traces and more easily plot
    n.trace = length(unique(fcst.sub$Trace.Number))
    fcst.mat = matrix(NA, nrow = dt, ncol = n.trace)
    # fcst.mat.full = matrix(NA, nrow = dt.full, ncol = n.trace)
    
    i = 1
    for(i in 1:n.trace){
      
      fcst.i = subset(fcst.sub$Slot.Value, fcst.sub$Trace.Number == i)
      
      if(meanLength > 5){
        fcst.i = fcst.i[1:60]
      }
      
      
      if(length(fcst.i) > 0 & is.na(fcst.i[1]) == F){
        # cat("\n")
        # cat("No missing data")
        # cat("\n")
        
        fcst.mat[,i] = fcst.i
        
        } else if (length(fcst.i) == 0 | is.na(fcst.i[1]) == T){  #if trace is missing, fill in with trace i + 1 data
        #fcst.mat[,i] = NA
          cat("\n")
          cat("!! Missing data !!")
          cat("\n")
          
          if(meanLength > 5){
            fcst.mat[,i] = subset(fcst.sub$Slot.Value, fcst.sub$Trace.Number == i+1)[1:60]
            }else{
            fcst.mat[,i] = subset(fcst.sub$Slot.Value, fcst.sub$Trace.Number == i+1)
            }
       }
      
    }
    
    # #######################################################################
    #read in historical simulation - 1 trace
    file = paste0(dir, "Historical/Historical_DMI,Historical MRM,MTOM Model,Rules,HistoricalStreamflow-", year, "-10/", reservoir, "Output.csv")  
    hist = read.csv(file)
    #subset to pool elevation, units are ft
    hist.sub = subset(hist$Slot.Value, hist$Object.Slot == paste0(reservoir, ".Pool Elevation"))
    hist.sub = hist.sub[1:(meanLength.N*12)]
    
    #t = clim$Timestep[1:60]
    
    #########################################################################
    #plot - novel fcst - all months - spaghetti plots
    sd = as.Date(paste0(year,'-10-01'))
    ed = as.Date(paste0(year+meanLength.N,'-09-01'))
    dates = seq(sd, ed, by = "month")
    # if(meanLength>5){
    #   ed.full = as.Date(paste0(year+meanLength,'-09-01'))
    #   dates.full = seq(sd, ed.full, by = "month")
    # }
    
    path = paste0("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Code/Phase 1/mtom/plots/MTOM_RF300/", fcstType, "_", meanLength, "yr/")
    dir.create(path)
    file_name = paste0(path, reservoir, "_spaghettiPlot_", sd, "_to_", ed, ".png")
    
    png(file=file_name, width = 8, height = 6, units = "in", res = 400)
    
    col1 = "royalblue"
    col2 = rgb(0,0,0, alpha = 0.7)
    
    title = paste0("Lake ", reservoir, " pool elevation for water years ", year+1, "-", year+meanLength.N, "\nRandom Forest (300 traces disaggregated from ", meanLength, "-year means)")
    #title = paste("Initial pool elevation = 1212 ft")
    par(mar = c(5, 5, 3, 2.5))
    # if(meanLength > 5){
    #   plot(x = dates.full, y = fcst.mat.full[,1], type = "l", col = col1, ylim = range(esp.mat, fcst.mat.full, hist.sub, na.rm=T), xlab = "Time", ylab = "Elevation (ft)", main = title, font.lab = 2, cex.main = 1)
    #   for(i in 2:n.trace){lines(dates.full, fcst.mat.full[,i], col = col1, lwd = 2)}
    # } else{
        plot(x = dates, y = fcst.mat[,1], type = "l", col = col1, ylim = range(esp.mat, fcst.mat, hist.sub, na.rm=T), xlab = "Time", ylab = "Elevation (ft)", main = title, font.lab = 2, cex.main = 1)
        for(i in 2:n.trace){lines(dates, fcst.mat[,i], col = col1, lwd = 2)}
    # }
      abline(v = dates[grep("-09-", dates)], lwd = 1, col = rgb(0,0,0, alpha = 0.9), lty = 2)
    
    #for(i in 2:n.trace.esp){lines(dates, esp.mat[,i], col = "grey", lwd = 2)}
    # #plot fcst 
    # col = rgb(0,1,0, alpha = 0.5)
    # for(i in 1:n.trace){lines(dates, fcst.mat[,i], col = col, lwd = 2)}
    # #replot observed
    # lines(dates, hist.sub, col = "red", lwd = 2)
    # #legend
    # legend("bottomleft", legend=c("Observed", "ESP", "Random Forest"), col = c("red", "grey", "green"), lwd = 2, cex = 1, box.lty=0, inset = c(0.05,0))
    #plot fcst 
    #col1 = rgb(0,1,0, alpha = 1)

    for(i in 1:n.trace.esp){lines(dates, esp.mat[,i], col = col2, lwd = 2)}
    #replot observed
    lines(dates, hist.sub, col = "red", lwd = 2)
    
    #legend
    legend("bottomleft", legend=c("Observed", "ESP", "Random Forest"), col = c("red", col2, col1),
           lwd = 2, cex = 1, box.lty=0, inset = c(0.05,0), bg=rgb(1,1,1, alpha = 0.75))
    dev.off()
    
    #########################################################################
    #plot - novel fcst - September only - boxplots
    sd = as.Date(paste0(year,'-10-01'))
    ed = as.Date(paste0(year+meanLength.N,'-09-01'))
    dates = seq(sd, ed, by = "month")
    title = paste0("Lake ", reservoir, " pool elevation for water years ", year+1, "-", year+meanLength.N, "\nRandom Forest")
    #title = paste("Initial pool elevation = 1212 ft")
    
    melt.fcst = melt(fcst.mat)
    melt.fcst$Method = "RF"
    melt.esp = melt(esp.mat)
    melt.esp$Method = "ESP"
    melt.hist = data.frame(X1 = 1:(meanLength.N*12), X2 = rep(NA, (meanLength.N*12)), value = hist.sub, Method = rep("Historic", meanLength.N*12))
    colnames(melt.hist)[1:2] = colnames(melt.fcst)[1:2]
    
    melt.comb = rbind(melt.fcst, melt.esp, melt.hist)
    colnames(melt.comb) = c("SimMonth", "Trace", "Elev", "Method")
    
    #subset so only September months are left
    septMonths = seq(12, meanLength.N*12, 12) 
    
    melt.comb.sub = filter(melt.comb, SimMonth %in% septMonths)
    melt.comb.sub$Year = rep(1:meanLength.N, nrow(melt.comb.sub)/meanLength.N)
    titleEOY = paste0("Lake ", reservoir, " September pool elevation for water years ", year+1, "-", year+meanLength.N)
    
    #boxplots
    plot = ggplot() + 
      geom_boxplot(data = subset(melt.comb.sub, Method != "Historic"), aes(x = factor(SimMonth), y = Elev, group = interaction(Year, Method), fill = Method)) +
      geom_crossbar(linetype = "dashed", data =  subset(melt.comb.sub, Method == "Historic"), aes(x = factor(SimMonth), y = Elev, ymin=Elev, ymax=Elev, color = Method), position=position_dodge(), color="black", fatten = 1.5) + 
      theme_bw() + ylab("Elevation (ft)") + theme(axis.title = element_text(face = "bold")) + xlab("Lead Time (Number of Months)") +
      ggtitle(titleEOY)
    
    print(plot)
    ggsave(paste0(titleEOY, ".png"), plot = print(plot), path = path, device = "png", width = 10, height = 5, dpi = 450, units = "in")
    
    #violin plots
    dodge = position_dodge(width = 0.8)
    plot = ggplot() + 
      geom_violin(position = dodge, alpha = 0.5, data = subset(melt.comb.sub, Method != "Historic"), aes(x = factor(SimMonth), y = Elev, group = interaction(Year, Method), fill = Method)) +
      geom_boxplot(width=0.1, outlier.colour=NA, position = dodge, data = subset(melt.comb.sub, Method != "Historic"), aes(x = factor(SimMonth), y = Elev, group = interaction(Year, Method), fill = Method)) +
      geom_crossbar(linetype = "dashed", data = subset(melt.comb.sub, Method == "Historic"), aes(x = factor(SimMonth), y = Elev, ymin=Elev, ymax=Elev, color = Method), position=position_dodge(), color="black", fatten = 1.5) + 
      theme_bw() + ylab("Elevation (ft)") + theme(axis.title = element_text(face = "bold")) + xlab("Lead Time (Number of Months)") +
      ggtitle(titleEOY)
    
    print(plot)
    ggsave(paste0(titleEOY, "_violinPlot.png"), plot = print(plot), path = path, device = "png", width = 10, height = 5, dpi = 450, units = "in")
    
    #boxplot with dots plots
    dodge = position_dodge(width = 1)
    plot = ggplot() +
      geom_quasirandom(groupOnX = TRUE, width=.2, alpha = 1, dodge.width = 1, data = subset(melt.comb.sub, Method != "Historic"), aes(x = factor(SimMonth), y = Elev, group = interaction(Year, Method), color = Method)) +
      geom_boxplot(width=.4, alpha = 0.85, outlier.shape = NA, position = dodge, data = subset(melt.comb.sub, Method != "Historic"), aes(x = factor(SimMonth), y = Elev, group = interaction(Year, Method), color = Method)) +
      geom_crossbar(linetype = "dashed", data =  subset(melt.comb.sub, Method == "Historic"), aes(x = factor(SimMonth), y = Elev, ymin=Elev, ymax=Elev, color = Method), position=position_dodge(), color="black", fatten = 1.5) +
      theme_bw() + ylab("Elevation (ft)") + theme(axis.title = element_text(face = "bold")) + xlab("Lead Time (Number of Months)") +
      ggtitle(titleEOY)

    print(plot)
    ggsave(paste0(titleEOY, "_boxplotWithDots.png"), plot = print(plot), path = path, device = "png", width = 10, height = 5, dpi = 450, units = "in")

    
    # print(
    #   ggplot() + 
    #     geom_boxplot(data = subset(melt.comb.sub, Method != "Historic"), aes(x = factor(SimMonth), y = Elev, group = interaction(Year, Method), color = Method)) +
    #     geom_point(data = subset(melt.comb.sub, Method == "Historic"), aes(x = factor(SimMonth), y = Elev, color = Method)) + 
    #     theme_bw() + ylab("Elevation (ft)") + xlab("Lead Time (Number of Months)") +
    #     ggtitle(titleEOY)
    # )
    
    
    #######################################################################
    #NSE calculations
    #######################################################################
    if(nseOn == T){
      #ESP
      #nse.esp = ens.nse(esp.mat, hist.sub, clim, meanLength.N)
      
      #fcst-BB with ESP ensemble mean as climatology
      nse.fcst = ens.nse(fcst.mat, hist.sub, esp.mat, meanLength.N, septOnly, septMonths)
    }
    
    #######################################################################
    #RMSE calculations
    #######################################################################
    
    #ESP
    rmse.esp = ens.rmse(esp.mat, hist.sub, septOnly = T, dates, meanLength.N, crpssOn = F, clim = NA)
    
    #fcst-BB
    rmse.fcst = ens.rmse(fcst.mat, hist.sub, septOnly = T, dates, meanLength.N, crpssOn = T, clim = esp.mat)
    
    ###########################################################################
    #merge and plot - NSE
    ###########################################################################
    if(nseOn == T){
      #NSE for all 5 years
      nse.list = list(nse.fcst[[2]])
      names(nse.list) = c("RF")
      title = paste0(meanLength.N, "-year NSE by forcing for water years ", year+1, "-", year+meanLength.N)
      
      par(font.lab = 2)
      boxplot(nse.list, ylab = "NSE", xlab = "Hydrologic Forcing Method",  boxwex=0.5, main = title)
      abline(h = 0)
    
      #####################
      #NSE for each year
      nse.ann.list = list(nse.fcst[[1]])
      names(nse.ann.list) = c("RF")
      title = paste0("Lake ", reservoir, " Elev. Annual NSE by forcing for water years ", year+1, "-", year+meanLength.N)
      
      nse.ann.melt = melt.list(nse.ann.list)
      colnames(nse.ann.melt) = c("Trace", "Year", "NSE", "Method")
      nse.ann.melt$Year = as.factor(nse.ann.melt$Year)
      
      
      print(ggplot(data = nse.ann.melt, aes(x = Year, y = NSE)) + 
              geom_point(aes(fill=Method)) + 
              coord_cartesian(ylim = c(-5, 1)) +
              geom_hline(yintercept = 0) +
              ggtitle(title))
      
      nse.ann.melt$sim.start.year = year
      
      #####################
      #year 2-4 avg nse
      nse.24.list = list(nse.fcst[[3]])
      names(nse.24.list) = c("RF")
      title = paste0("Year 2-", meanLength.N, " monthly average NSE by forcing and trace\nWater years ", year+2, "-", year+meanLength.N-1)
      nse.24.melt = melt.list(nse.24.list)
      
      par(font.lab = 2)
      boxplot(nse.24.list, ylab = "NSE", xlab = "Hydrologic Forcing Method",  boxwex=0.5, main = title, cex.main = 0.9, ylim = c(-1,1))
      abline(h = 0)
      
      #ensemble NSE
      title = paste0("EOWY NSE\nWater years ", year+1, "-", year+meanLength.N)
      plot(septMonths, nse.fcst[[4]], ylab = "NSE", main = title)
      
      #ensemble mean NSE
    }
    
    ###########################################################################
    #merge and plot - RMSE
    ###########################################################################
    #RMSE for all 5 years
    rmse.list = list(rmse.esp[[2]], rmse.fcst[[2]])
    names(rmse.list) = c("ESP", "RF")
    title = paste0(meanLength.N, "-year RMSE by forcing for water years ", year+1, "-", year+meanLength.N)
    
    par(font.lab = 2)
    boxplot(rmse.list, ylab = "RMSE (ft)", xlab = "Hydrologic Forcing Method",  boxwex=0.5, main = title)
    abline(h = 0)
    
    #####################
    #rmse for each year
    rmse.ann.list = list(rmse.esp[[1]], rmse.fcst[[1]])
    names(rmse.ann.list) = c("ESP", "RF")
    title = paste0("Lake ", reservoir, " Elev. Annual RMSE (ft) by forcing for water years ", year+1, "-", year+meanLength.N)
    
    rmse.ann.melt = melt.list(rmse.ann.list)
    colnames(rmse.ann.melt) = c("Trace", "Year", "RMSE", "Method")
    rmse.ann.melt$Year = as.factor(rmse.ann.melt$Year)
    
    
    print(ggplot(data = rmse.ann.melt, aes(x = Year, y = RMSE)) + 
            geom_boxplot(aes(fill=Method)) + 
            #coord_cartesian(ylim = c(-5, 1)) +
            geom_hline(yintercept = 0) +
            ggtitle(title))
    
    rmse.ann.melt$sim.start.year = year
    
    #####################
    #year 2-4 avg rmse
    rmse.24.list = list(rmse.esp[[3]], rmse.fcst[[3]])
    names(rmse.24.list) = c("ESP", "RF")
    title = paste0("Year 2-", meanLength.N," monthly average RMSE by forcing and trace\nWater years ", year+2, "-", year+meanLength.N-1)
    rmse.24.melt = melt.list(rmse.24.list)
    
    par(font.lab = 2)
    boxplot(rmse.24.list, ylab = "RMSE (ft)", xlab = "Hydrologic Forcing Method",  boxwex=0.5, main = title, cex.main = 0.9)
    abline(h = 0)
    
    #####################
    #ens mean rmse
    rmse.ens.mean.list = list(rmse.esp[[4]], rmse.fcst[[4]])
    names(rmse.ens.mean.list) = c("ESP", "RF")
    #title = paste0("Ensemble mean RMSE by forcing and trace\nWater years ", year+1, "-", year+5)
    rmse.ensMean.melt = melt.list(rmse.ens.mean.list)
    
    #####################
    #September RMSE
    if(septOnly == T){
      rmse.eoy.list = list(rmse.esp[[5]], rmse.fcst[[5]])
      names(rmse.eoy.list) = c("ESP", "RF")

      title = paste0(meanLength.N, "-year RMSE (September only) by forcing for water years ", year+1, "-", year+meanLength.N)

      rmse.eoy.melt = melt.list(rmse.eoy.list)
      
      par(font.lab = 2)
      boxplot(rmse.eoy.list, ylab = "RMSE (ft)", xlab = "Hydrologic Forcing Method",  boxwex=0.5, main = title, cex.main = 0.9)
      abline(h = 0)
      
      #ensemble rmse
      #by ensemble year
      rmse.eoy.byEns.list = list(rmse.esp[[6]], rmse.fcst[[6]])
      names(rmse.eoy.byEns.list) = c("ESP", "RF")
      rmse.eoy.byEns.melt = melt.list(rmse.eoy.byEns.list)
      
      #eoy crpss
      title = paste0("EOWY CRPSS (ESP as climatology)\nWater years ", year+1, "-", year+meanLength.N)
      plot(1:meanLength.N, rmse.fcst[[7]], ylab = "CRPSS", main = title)
      abline(h = 0)
      #axis(1, at=septMonths)
      
      crpss.eoy = melt.list(rmse.fcst[[7]])
      
    }
    
    ############################################
    #insert into aggregated dfs
    #NSE
    if(nseOn == T){
      nse.5year.list[[y]] = nse.list
      nse.ann.df = rbind(nse.ann.df, nse.ann.melt)
      nse.24.df = rbind(nse.24.df, nse.24.melt)
    }
    #RMSE
    rmse.5year.list[[y]] = rmse.list
    rmse.ann.df = rbind(rmse.ann.df, rmse.ann.melt)
    rmse.24.df = rbind(rmse.24.df, rmse.24.melt)
    rmse.ensMean.df = rbind(rmse.ensMean.df, rmse.ensMean.melt)
    if(septOnly == T){
      rmse.eoy.df = rbind(rmse.eoy.df, rmse.eoy.melt)
      rmse.eoy.byEns.df = rbind(rmse.eoy.byEns.df, rmse.eoy.byEns.melt)
      crpss.eoy.df = rbind(crpss.eoy.df, crpss.eoy)
    }
    
    y = y + 1
    
  }
  
  if(nseOn == T){
    names(nse.5year.list) = years
    names(nse.24.df) = c("NSE", "Method")
  }
  names(rmse.24.df) = c("RMSE", "Method")
  names(rmse.ensMean.df) = c("RMSE", "Method")
  
  if(septOnly == T){
    rmse.eoy.df = rmse.eoy.df[,-c(1:2)]  
    names(rmse.eoy.df) = c("RMSE", "Method")
    
    rmse.eoy.byEns.df = rmse.eoy.byEns.df[,-2]
    names(rmse.eoy.byEns.df) = c("Year", "RMSE", "Method")
    rmse.eoy.byEns.df$Year = as.factor(rmse.eoy.byEns.df$Year)
    
    crpss.eoy.df$Method = "RF"
    names(crpss.eoy.df)[1:2] = c("CRPSS", "Year")
    crpss.eoy.df$Year = as.factor(crpss.eoy.df$Year)
  }
  
  #assign function outputs
  out = list()
  #NSE
  if(nseOn == T){
    out$nse.5year.list = nse.5year.list
    out$nse.ann.df = nse.ann.df
    out$nse.24.df = nse.24.df
  }
  #RMSE
  out$rmse.5year.list = rmse.5year.list
  out$rmse.ann.df = rmse.ann.df
  out$rmse.24.df = rmse.24.df
  out$rmse.ensMean.df = rmse.ensMean.df
  if(septOnly == T){
    out$rmse.eoy.df = rmse.eoy.df
    out$rmse.eoy.byEns.df = rmse.eoy.byEns.df
    out$crpss.eoy.df = crpss.eoy.df
  }
  
  ####################
  ##NSE
  #####
  if(nseOn == T){
  #all 5 year blocks
    title = paste0("Lake ", reservoir, " Elev. Annual NSE by forcing for ", length(years), " different ", meanLength, "-year blocks\n", head(years, 1), "-", tail(years, 1))
    
    print(ggplot(data = nse.ann.df, aes(x = Year, y = NSE)) + 
            geom_boxplot(aes(fill=Method)) + 
            coord_cartesian(ylim = c(-40, 1)) +
            geom_hline(yintercept = 0) +
            ggtitle(title))
    
    p0 = ggplot(data = nse.ann.df, aes(x = Year, y = NSE)) + 
            geom_boxplot(aes(fill=Method)) + theme_bw() + 
            coord_cartesian(ylim = c(-5, 1)) +
            geom_hline(yintercept = 0) +
            ggtitle(title) +
            theme(axis.title = element_text(face = "bold"))
      

    print(p0)
    
    ggsave(paste0(reservoir, ".Years1-", meanLength, "_AnnualNSE.png"), p0, "png", path = path, dpi = 400)
    
    
    #years 2-4
    title = paste0("Lake ", reservoir, " Elev. year 2-4 average monthly NSE by forcing for ", length(years), " different ", meanLength, "-year blocks\n", head(years, 1), "-", tail(years, 1))
    print(ggplot(data = nse.24.df, aes(x = Method, y = NSE)) + 
            geom_boxplot(aes(fill=Method)) + 
            #coord_cartesian(ylim = c(-5, 1)) +
            geom_hline(yintercept = 0) +
            ggtitle(title) + 
            theme(plot.title = element_text(size = 10, face = "bold")))
  }
  
  ####################
  ##RMSE
  #####
  #all 5 year blocks
  #title = paste0("Lake ", reservoir, " Elevation Annual RMSE by forcing for ", length(years), " different 5-year blocks\n", head(years, 1), "-", tail(years, 1))
  title = paste0("Lake ", reservoir, " Elevation Annual RMSE by forcing entire hindcast", " (", fcstType, "- ", meanLength, "yr)")
  p1 = ggplot(data = rmse.ann.df, aes(x = Year, y = RMSE)) + theme_bw() +
    geom_boxplot(aes(fill=Method)) + 
    #coord_cartesian(ylim = c(-40, 1)) +
    geom_hline(yintercept = 0) + ylab('RMSE (ft)') +
    theme(axis.title = element_text(face = "bold")) +
    ggtitle(paste0("Lake ", reservoir))
  
  print(p1)
  ggsave(paste0(title, ".png"), plot = print(p1), path = path, device = "png", width = 10, height = 5, dpi = 450, units = "in")
  
  
  #years 2-4
  title = paste0("Lake ", reservoir, " elev. year 2-4 average monthly RMSE by forcing for ", length(years), " different ", meanLength, "-year blocks\n", head(years, 1), "-", tail(years, 1), " (", fcstType, ")")
  print(ggplot(data = rmse.24.df, aes(x = Method, y = RMSE)) + 
          geom_boxplot(aes(fill=Method)) + 
          #coord_cartesian(ylim = c(-5, 1)) +
          geom_hline(yintercept = 0) +
          ggtitle(title) + ylab('RMSE (ft)') +
          theme(plot.title = element_text(size = 10, face = "bold")))
  
  #ens Mean
  title = paste0("Lake ", reservoir, " elev. ensemble mean RMSE by forcing for ", length(years), " different ", meanLength, "-year blocks\n", head(years, 1), "-", tail(years, 1), " (", fcstType, ")")
  p3 = ggplot(data = rmse.ensMean.df, aes(x = Method, y = RMSE)) + 
          geom_boxplot(aes(fill=Method)) + 
          #coord_cartesian(ylim = c(-5, 1)) +
          geom_hline(yintercept = 0) +
          ggtitle(title) + ylab('RMSE (ft)') +
          theme(plot.title = element_text(size = 10, face = "bold"))
  print(p3)
  
  ggsave(paste0(reservoir, ".Years1-", meanLength, "_EnsembleMeanRMSE.png"), p3, "png", path = path, dpi = 400)
  
  #September only
  if(septOnly == T){
    
    rmseMonth = "September"

    title = paste0("Lake ", reservoir, " elevation year 1-", meanLength.N, " ", rmseMonth," RMSE by forcing for ", length(years), " different ", meanLength, "-year blocks\n", head(years, 1), "-", tail(years, 1), " (", fcstType, ")")
    
    p4 = ggplot(data = rmse.eoy.df, aes(x = Method, y = RMSE)) + 
      geom_boxplot(aes(fill=Method)) + 
      #coord_cartesian(ylim = c(-5, 1)) +
      geom_hline(yintercept = 0) +
      ggtitle(title) + ylab('RMSE (ft)') +
      theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(face = "bold"))
    
    print(p4)
    ggsave(paste0(reservoir, ".Years1-", meanLength, "_SeptemberOnlyRMSE.png"), p4, "png", path = path, dpi = 400)
    
    #by ens mems
    title = paste0("Lake ", reservoir, " elev. years 1-", meanLength.N, " ", rmseMonth," ensemble RMSE by forcing for ", length(years), " different ", meanLength, "-year blocks\n", head(years, 1), "-", tail(years, 1), " (", fcstType, ")")
    
    p5 = ggplot(data = rmse.eoy.byEns.df, aes(x = Year, y = RMSE)) + 
      geom_boxplot(aes(fill = Method)) + theme_bw() + 
      coord_cartesian(ylim = c(0, 110)) +
      geom_hline(yintercept = 0) +
      ggtitle(title) + ylab('RMSE (ft)') + xlab("Lead Time (number of months)") +
      theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(face = "bold")) +
      scale_x_discrete(labels = seq(12, (meanLength.N*12), 12))

    
    print(p5)
    
    ggsave(paste0(reservoir, ".Years1-", meanLength, "_SeptemberOnlyRMSEbyEns.png"), p5, "png", path = path, dpi = 400)
    
    d5 = ggplot_build(p5)$data
    print(d5[[1]])
    
    #write.csv(d5, paste0(path, fcstType, "-", reservoir, "-", meanLength, "yr_SeptemberOnlyRMSEbyEns_boxplotSummary.xslx"), sep = ",")
    saveRDS(d5, paste0(path, fcstType, "_", reservoir, "_", meanLength, "yr_SeptemberOnlyRMSEbyEns_boxplotSummary.rds"))
    
    #EOY CRPSS
    title = paste0("Lake ", reservoir, " elev. years 1-", meanLength.N, " ", rmseMonth," CRPSS for ", length(years), " different ", meanLength, "-year blocks\nESP as climatology - ", head(years, 1), "-", tail(years, 1), " (", fcstType, ")")
    
    p6 = ggplot(data = crpss.eoy.df, aes(x = Year, y = CRPSS)) + 
      geom_boxplot(aes(fill = Method)) + theme_bw() + 
      #coord_cartesian(ylim = c(-5, 1)) +
      geom_hline(yintercept = 0) +
      ggtitle(title) + ylab('CRPSS') + xlab("Lead Time (number of months)") +
      theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(face = "bold")) +
      scale_x_discrete(labels = seq(12, (meanLength.N*12), 12))

    print(p6)
    
    p6 = ggplot(data = crpss.eoy.df, aes(x = Year, y = CRPSS)) + 
      geom_boxplot(aes(fill = Method)) + theme_bw() + 
      coord_cartesian(ylim = c(-5, 1)) +
      geom_hline(yintercept = 0) +
      ggtitle(title) + ylab('CRPSS') + xlab("Lead Time (number of months)") +
      theme(plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(face = "bold")) +
      scale_x_discrete(labels = seq(12, (meanLength.N*12), 12))
    
    print(p6)
    
    ggsave(paste0(reservoir, ".Years1-", meanLength, "_SeptemberOnlyCRPSS.png"), p6, "png", path = path, dpi = 400)
    
  }
  
  
  return(out)
  
}

