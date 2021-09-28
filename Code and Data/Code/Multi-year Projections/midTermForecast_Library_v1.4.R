#midterm forecast functions library
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
library(reshape)
library(stringr)
library(ggpubr)
library(randomForest)
library(tidytext)
library(DescTools)
library(viridis)
library(scatterplot3d)
library(xts)
library(ggridges)
#####################################################################################################################
#functions ##########################################################################################################
#####################################################################################################################

#####################################################################################################################
#####################################################################################################################

#define flow averaging function for past- and future flows
flow.avg = function(stats, i, meanLength){
  #add wyear col
  stats.Nyr.corSearch = data.frame(wyears = stats$wyears[(i+1):(nrow(stats)-(meanLength-1))])
  
  #create future 5 year means - remove first 5 years to match with later calculation of past 5 year mean
  stats.Nyr.corSearch$q.f_maf = rollmean(stats$q_maf[-(1:i)], meanLength, align = "left")
  
  #past 5 year mean (mean of scaled anomalies) - remove last 5 years of flow
  #stats.Nyr$q.p = rollmean(head(stats$q_maf, -5), k = 5)
  stats.Nyr.corSearch$q.p = rollmean(head(stats$q_maf, -meanLength), k = i)
  
  return(stats.Nyr.corSearch)
  
}

#####################################################################################################################
#####################################################################################################################

#define flow averaging function for past RE and future flows
flowRE.avg = function(stats, i, meanLength, LE.ann){
  
  stats.sub = subset(stats, wyears >= min(LE.ann$wyears))  
  
  #add wyear col
  stats.Nyr.corSearch = data.frame(wyears = stats.sub$wyears[(i+1):(nrow(stats.sub)-(meanLength-1))])
  
  #create future 5 year means - remove first i years to match with later calculation of past i year mean
  stats.Nyr.corSearch$q.f_maf = rollmean(stats.sub$q_maf[-(1:i)], meanLength, align = "left")
  
  #past i year mean (mean of scaled anomalies) - remove last 5 years of flow
  #stats.Nyr$q.p = rollmean(head(stats$q_maf, -5), k = 5)
  #stats.Nyr.corSearch$re.p = rollmean(head(stats.sub$re.model, - 5), k = i)
  stats.Nyr.corSearch$re.p = rollmean(head(stats.sub$re, - meanLength), k = i)
  
  return(stats.Nyr.corSearch)
  
}

#####################################################################################################################
#####################################################################################################################

#df = dp
mean.anoms = function(df, meanLength, stats.Nyr, window.length.re){
  
  #optional, subset DPLE.ann to remove pre-19XX years (appears to be less skill in those years)
  df = subset(df, wyears <= max(stats.Nyr$wyears))
  
  #subset obs 5 year mean flows to start in same year as the precip and temp covariates or vice versa
  if(nrow(stats.Nyr) >= nrow(df)){
    obs.q.Nyr = subset(stats.Nyr, wyears >= df$wyears[1])
  }
  
  if(nrow(stats.Nyr) <= nrow(df)){
    obs.q.Nyr = stats.Nyr
    df = subset(df, wyears >= stats.Nyr$wyears[1])
  }
  
  #or if LE record is longer than flow record (e.g. long flow averaging times) subset LE to start same year as flow
  #df = subset(df, wyears >= obs.q.Nyr$wyears[1])
  
  #initialize data frame to store 5-year mean anomalies
  df.meanAnoms = NULL
  df.means = NULL
  
  i = 2
  for(i in 2:ncol(df)){
    
    #or calc scale anomalies (for use in KNN, should use scaled covariates when calculating distances)
    #anoms.i = scale(df[,i])
    
    #calc rolling mean on scaled anomalies
    #meanAnoms.i = rollmean(anoms.i, k = 5)
    
    #calc rolling mean on unscaled data for use in RF
    means.i = df[,i]
    
    #calc scaled anomalies of rolling means
    meanAnoms.i = scale(means.i)
    
    #aggregate
    df.meanAnoms = cbind(df.meanAnoms, meanAnoms.i)
    df.means = cbind(df.means, means.i)
  }
  
  #add water years and colnames
  meanAnoms = data.frame(df$wyears[1:nrow(df.meanAnoms)], df.meanAnoms)
  colnames(meanAnoms) = colnames(df)
  
  means = data.frame(df$wyears[1:nrow(df.means)], df.means)
  colnames(means) = colnames(df)
  
  #add past (p) i-year mean flows and runoff efficiencies
  #note that the past flows are scaled anomalies (so that covariates are all scaled anomalies)
  #meanAnoms$q.p = as.vector(scale(obs.q.Nyr$q.p))
  meanAnoms$re.p = as.vector(scale(obs.q.Nyr$re.p))
  
  #means$q.p = as.vector(obs.q.Nyr$q.p)
  means$re.p = as.vector(obs.q.Nyr$re.p)
  
  #calc correlation between each covariate and 5-year mean flow 
  mutInf = matrix(NA, nrow = ncol(meanAnoms)-1, ncol = 1)
  for(i in 2:ncol(meanAnoms)){
    mutInf[i-1,1] = MutInf(meanAnoms[,i], obs.q.Nyr$q.f_maf)
  }
  
  meanAnomsFlowCor = data.frame(cor = mutInf, var = as.character(c(colnames(df[-1]), "re.p")))
  #meanAnomsFlowCor = data.frame(cor = apply(meanAnoms[,-1], 2, cor, obs.q.Nyr$q.f_maf), var = as.character(c(colnames(df[-1]), "re.p")))
  
  #calc weights from correlations
  meanAnomsFlowCor$wts = abs(meanAnomsFlowCor$cor)/sum(abs(meanAnomsFlowCor$cor))
  
  #add future (f) 5-year mean flows (MAF, not anomalies)
  meanAnoms$q.f_maf = as.vector(obs.q.Nyr$q.f_maf)
  means$q.f_maf = as.vector(obs.q.Nyr$q.f_maf)
  
  #plot of flow vs each covariate
  #par(pty = "s")
  #plot(q.f_maf ~ ., data = meanAnoms)
  
  n = ncol(meanAnoms)-2
  plot.list = vector("list", n)
  
  i = 2
  for(i in 2:(ncol(meanAnoms)-1)){
    
    if(i == 2){ylab = paste0(meanLength, "-year mean observed flow (MAF)")}  
    else(ylab = "")
    
    if(str_detect(colnames(meanAnoms)[i], "Tmax") == T){xlab = paste0(meanLength, "-year mean model Tmax\n(scaled anomaly)")}  
    if(str_detect(colnames(meanAnoms)[i], "Tmin") == T){xlab = paste0(meanLength, "-year mean model Tmin\n(scaled anomaly)")}  
    if(str_detect(colnames(meanAnoms)[i], "pcp") == T){xlab = paste0(meanLength, "-year mean model Pcp\n(scaled anomaly)")}  
    #if(str_detect(colnames(meanAnoms)[i], "q.p") == T){xlab = paste0("Past ", window.length.qp, "-year mean flow\n(scaled anomaly)")}  
    if(str_detect(colnames(meanAnoms)[i], "re.p") == T){xlab = paste0("Past ", window.length.re, "-year mean runoff eff.\n(scaled anomaly)")}  
    
    plot = ggplot(meanAnoms, aes_string(x=meanAnoms[,i], y = meanAnoms$q.f_maf, color = meanAnoms$wyears)) +
      geom_point() + xlab(xlab) + ylab(ylab) + theme_bw() + theme(aspect.ratio=1) +
      # scale_color_gradient(low="blue", high="red") +
      scale_colour_gradientn(colours = hcl.colors(40, palette = "viridis")) +
      labs(colour = "Year") +
      geom_smooth(method='loess', se=TRUE, fullrange=T, level=0.95, col = "black") +
      #geom_smooth(method="auto", se=TRUE, fullrange=T, level=0.95, col = "black") +
      #ggtitle(paste0("Mutual Information = ", round(meanAnomsFlowCor[i-1,1], digits = 2))) +
      theme(axis.title = element_text(face = "bold"), text=element_text(size=16))
    
    plot.list[[i-1]] = plot
    
  }
  
  #print(ggarrange(plot.list, ncol = n, nrow = 1, common.legend = T, legend = "right"))
  
  out = list(meanAnoms, meanAnomsFlowCor, plot.list, n, means)
  names(out) = c("meanAnoms", "meanAnomsFlowCor", "plots", "NumCovars", "means")
  
  return(out)
}

#####################################################################################################################
#####################################################################################################################

### Define KNN function
knn = function(fit_updated, pred.j, wts){
  
  k = nrow(fit_updated)
  
  np = ncol(fit_updated) - 2
  dist.mat = matrix(NA, nrow = nrow(fit_updated), ncol = np)
  
  l = 1
  for(l in 1:ncol(dist.mat)){
    dist.mat[,l] = abs(wts[l])*(pred.j[1,l+1] - fit_updated[,l+1])^2
  }
  
  #calc distance for each row
  fit_updated$dist = sqrt(apply(dist.mat, 1, sum))
  
  #select the k smallest distances
  neighbs0 = top_n(fit_updated, -1*k, wt = dist)
  #re-order by distance
  neighbs = neighbs0[order(neighbs0$dist, decreasing = F),]
  
  return(neighbs)
  
}    

#####################################################################################################################
#####################################################################################################################

#calc flow terciles for bootstrapped flows

flowTercileCalc = function(df, stats.Nyr){

  for(r in 1:nrow(df)){
    if(df$q.f_maf[r] >= as.numeric(quantile(stats.Nyr$q.f_maf, 2/3))){
      df$simFlowTercile[r] = "High"
      df$simFlowTercileNo[r] = 3
    } else if(df$q.f_maf[r] <= as.numeric(quantile(stats.Nyr$q.f_maf, 1/3))){
      df$simFlowTercile[r] = "Low"
      df$simFlowTercileNo[r] = 1
    }  else{
      df$simFlowTercile[r] = "Average"
      df$simFlowTercileNo[r] = 2
    }
  }

  return(df)
  
}


#####################################################################################################################
#####################################################################################################################

##############################################
### Skill score calcs
##########

###crpss and rpss on 5-year means
bootstrap = T
terciles = T
tempTerciles = F
plotByYear = T
flowTerciles = T

ss.calc.Nyr = function(sim.results, bootstrap, flowTerciles, plotByYear, tempTerciles, meanLength, vl, vs, ve, stats.Nyr, nsim, fcstMode){
 
  ## Probability of value in distribution
  quantInv <- function(distr, value) ecdf(distr)(value)
  ## Thresholds for Reliability Diagram
  thresholds = seq(0,1,0.2) #0.1 #c(0, .1, .9, 1)
  
  #initialize list to store crpss from different models, currently only 1 model in use
  models = unique(sim.results$`NyrMeanFlowBootstrap`$covariate)
  crpss.list = vector("list", length(models))
  df.reliability = NULL
  df.reliability.score = NULL
  df_rel.dia = NULL
  df_sub.agg = NULL
  
  #loop thru different models
  i = 1
  for(i in 1:length(models)){
    
    #select bootstrapped forecasts vs non-bootstrapped (neighbors only)
    if(bootstrap == T){fcsts.i = sim.results$`NyrMeanFlowBootstrap`}
    if(bootstrap == F){fcsts.i = sim.results$`NyrMeanFlow`}
    
    #subset forecasts to i-th model of loop
    fcsts.i = subset(fcsts.i, covariate == models[i])
    
    #initialize matrix to store CRPSS
    crpss.cov = matrix(NA, nrow = vl, ncol = 4)
    rh.mat.fcst = NULL
    rh.mat.obs = NULL
    rel.mat.fcst = NULL
    rel.mat.obs = NULL

    #loop through each year in validation (hindcast) period
    j = 1
    for(j in 1:vl){
      
      #subset forecasts to j-th start year of validation period
      fcsts.j = subset(fcsts.i, startYear == vs+j-1)
      
      #fetch 5-year mean climatology from observed record
      if(fcstMode == 'loocv'){
        clim = t(matrix(subset(stats.Nyr, wyears != fcsts.j$startYear[1])[,2]))
      } else if (fcstMode == 'kFold'){
        #clim = t(matrix(subset(stats.Nyr, wyears <= vs+j-1-meanLength | wyears >= vs+j-1+meanLength)[,2]))
        clim = t(matrix(tail(subset(stats.Nyr, wyears <= vs+j-1-meanLength | wyears >= vs+j-1+meanLength)[,2], 30)))
      } else if (fcstMode == 'retroBlind'){
        clim = t(matrix(subset(stats.Nyr, wyears <= vs+j-1-meanLength)[,2]))
      }
      
      #re-format ensemble data
      ens = t(matrix(fcsts.j$q.f_maf))
      
      #pull observed 5-year mean
      obs = t(matrix(fcsts.j$obsFlow[1]))
      
      #calc RPSS on forecast and climatology
      crps.f = EnsCrps(ens = ens, obs = obs)
      crps.c = EnsCrps(ens = clim, obs = obs)
      
      #calc CRPSS and store
      crpss = 1 - crps.f/crps.c
      
      #crpss decomposition for reliability
      # cd = crpsDecomposition(obs = obs, eps = ens)
      # reliability = cd$Reli
      
      crpss.cov[j,1] = crpss
      crpss.cov[j,2] = fcsts.j$flowTercile[1]    
      crpss.cov[j,3] = fcsts.j$tempTercile[1] 
      crpss.cov[j,4] = fcsts.j$startYear[1]
      # crpss.cov[j,5] = reliability
      
      #aggregate for rank histogram
      rh.mat.fcst = rbind(rh.mat.fcst, ens)
      rh.mat.obs = rbind(rh.mat.obs, obs)
      
      #aggregate for reliability
      fcst_prob = quantInv(ens, obs) # where the obs value lies in the forecast CDF
      obs_freq = quantInv(clim, obs) # quantile of obs flow in history
      
      rel.mat.fcst = rbind(rel.mat.fcst, fcst_prob)
      rel.mat.obs = rbind(rel.mat.obs, obs_freq)
     

    }
    
    #format for plotting
    df.crpss.cov = data.frame(crpss.cov)
    colnames(df.crpss.cov) = c("CRPSS", "flowTercile", "tempTercile", "startYear")
    #df.crpss.cov$covariate = paste0(covs[i])
    df.crpss.cov$CRPSS = as.numeric(as.character(df.crpss.cov$CRPSS))
    df.crpss.cov$startYear = df.crpss.cov$startYear
    # df.crpss.cov$reliability = as.numeric(as.character(df.crpss.cov$reliability))
    
    #crpss decomposition for reliability
    cd.hindcast = crpsDecomposition(obs = rh.mat.obs, eps = rh.mat.fcst)
    reliability.hindcast = cd.hindcast$Reli
    
    #calc rank histogram
    rh.fcst = Rankhist(ens = rh.mat.fcst, obs = rh.mat.obs, reduce.bins = 7)
    rh.fcst.plot = PlotRankhist(rh.fcst)
    #here I computed the ensemble percentile instead of the ensemble rank and each bin has 5% width
    Nint=100
    nt=seq(5/Nint, 1,5/Nint)
    Nint1=length(nt)
    #simul corresponds to a matrix of M ensembles (rows) for N days (columns)
    #data corresponds a vector of observation for N days
    Frequency = quantile(cumsum(rh.fcst), probs = nt)
    Frequency[2:Nint1] = Frequency[2:Nint1] - Frequency[1:(Nint1-1)]
    rankpl=Frequency#this has the frecuency for each bin
    #Di correspond to the discrepancy index see Delle Monache et al. (2006)
    DI.fcst = sum( sapply( 1:Nint1, function(irank) 100*abs(Frequency[irank]/sum(Frequency)-nt[1]) ) )
    cat("\n")
    cat(paste0(models[i], " Forecast DI = ", DI.fcst))
    cat("\n")
    
    #reliability plot -------------------------------------------------------------------------
    # put forecast probability into bins
    XX <- probcont2disc(rel.mat.fcst, bins = thresholds) # put fcst probability into threshold bins
    fcst_bin <- XX$new # forecast probability bins
    
    ## calculate forecast probability vs observed frequency
    N.pred <- aggregate(fcst_bin, by = list(fcst_bin), length) # count forecasts per bin
    N.obs <- aggregate(rel.mat.obs, by = list(fcst_bin), sum) # sum of obs forecast frequency per bin
    obar.i = N.obs/N.pred # calc average obs frequency in each bin - should correspond to bin median
    obar.i$Group.1 = N.obs$Group.1 # add thresholds to matrix
    obar.i$N.pred = N.pred
    
    #aggregate reliablilty 
    df.reli = data.frame(pred = rel.mat.fcst, obs_freq = rel.mat.obs)
    df.reli$Model = models[i]
      
    df.reliability = rbind(df.reliability, df.reli)
    
    df.reliability.score = rbind(df.reliability.score, reliability.hindcast)
    
    #aggregate reliablity for later plotting

    df_rel.dia = rbind.data.frame(df_rel.dia,
                                  cbind.data.frame(fcst_type = models[i],
                                                   threshold = obar.i$Group.1, # thresholds
                                                   obar = obar.i$V1, # observed relative frequency
                                                   n.pred = N.pred$x,
                                                   label = models[i])) # number of forecasts each bin

    df_sub = data.frame(x = rep(obar.i$Group.1, N.pred$x))
    df_sub$label = models[i]
    
    df_sub.agg = rbind(df_sub.agg, df_sub)
    
    
    
    ###############################################################################
    #calc RPSS on terciles
    
    #transform data into correct format for EnsRps function
    years = unique(fcsts.i$startYear)
    nYrs = length(years)
    mat.Nyr.fcst = matrix(NA, nrow = nYrs, ncol = nsim)
    mat.Nyr.obs = matrix(NA, nrow = nYrs, ncol = 1)
    
    j = 1
    for(j in 1:nYrs){
      
      df.j = subset(fcsts.i, startYear == years[j])
      
      mat.Nyr.fcst[j,] = df.j$simFlowTercileNo
      
      mat.Nyr.obs[j] = df.j$obsFlowTercile[1]

    }
    
    
    #calc terciles on entire historic flow series (N-year means)
    stats.Nyr.tercs = flowTercileCalc(stats.Nyr, stats.Nyr)
    nClim = nrow(stats.Nyr.tercs)
    
    #calc RPS time series (one score per year) for forecast and clim, then calc RPSS from the two scores
    if(fcstMode == 'loocv'){
      
      rpss.f = c()

      j = 1
      for(j in 1:nYrs){
      mat.fcst = t(matrix(mat.Nyr.fcst[j,]))  
      mat.clim = t(matrix(subset(stats.Nyr.tercs, wyears != vs+j-1)[,6]))
      
      rps.f = EnsRps(mat.fcst, mat.Nyr.obs[j])
      rps.c = EnsRps(mat.clim, mat.Nyr.obs[j])
      rpss.f[j] = 1 - rps.f/rps.c
      
      }
      
      
    } else if (fcstMode == 'kFold'){
      
      rpss.f = c()
      
      j = 1
      for(j in 1:nYrs){
        mat.fcst = t(matrix(mat.Nyr.fcst[j,]))  
        mat.clim = t(matrix(subset(stats.Nyr.tercs, wyears <= vs+j-1-meanLength | wyears >= vs+j-1+meanLength)[,6]))
        
        rps.f = EnsRps(mat.fcst, mat.Nyr.obs[j])
        rps.c = EnsRps(mat.clim, mat.Nyr.obs[j])
        rpss.f[j] = 1 - rps.f/rps.c
        
      }
      
    } else if (fcstMode == 'retroBlind'){
      
      rpss.f = c()
      
      j = 1
      for(j in 1:nYrs){
        mat.fcst = t(matrix(mat.Nyr.fcst[j,]))  
        mat.clim = t(matrix(subset(stats.Nyr.tercs, wyears <= vs+j-1-meanLength)[,6]))
        
        rps.f = EnsRps(mat.fcst, mat.Nyr.obs[j])
        rps.c = EnsRps(mat.clim, mat.Nyr.obs[j])
        rpss.f[j] = 1 - rps.f/rps.c
        
      }
      
    }
      
    df.crpss.cov$RPSS = rpss.f
    
    #insert results for i-th covariate into overall list
    crpss.list[[i]] = df.crpss.cov
    
  }
  
  ###################################################################################
  #post processing  
  names(crpss.list) = models
  skillScore.melt = melt.list(crpss.list)
  colnames(skillScore.melt)[ncol(skillScore.melt)] = "model"
  #colnames(crpss.melt)[(ncol(crpss.melt)-1):ncol(crpss.melt)] = c("CRPSS", "covariate")
  #crpss.melt$All.Years = 'All.Years'
  
  ####################################################################################################
  #########################################################################################
  #import and calc skill for ESP 
  

  if(meanLength <= 5){
    
    esp = readRDS("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/esp/allFcstLocs_ESP_fcst.rds")
    
    #
    #first generate linear model for later conversion from unreg to nat flows (obs unreg and obs nat used to train LM)
    
    #####################
    #process monthly data
    
    obs.unreg = esp$locAll_obs
    #remove first and last few months to create water year format (range is now water years 1980-2019)
    obs.unreg = obs.unreg[-c(1:9,490:492),]
    #reset row names
    rownames(obs) = NULL
    #rearrange sites (columns)
    obs.mod = obs.unreg[,c(1,12,3,4,6,5,7,8,13,2,9,10,11)]
    #convert to xts object
    obs.monthly = xts(obs.mod[,-1], order.by = obs.mod$Timestep)
    
    dt = 1980:2019
    n = length(dt)
    props = vector("list", n)
    
    seq = seq(1, 12*n, 12)
    ann.agg.flows = matrix(NA, nrow = n, ncol = 2)
    ann.agg.flows[,1] = dt
    colnames(ann.agg.flows) = c("yrs", "LeesFerry")
    
    i = 1
    j = 1
    for(i in seq){
      
      #pull wy i data
      o.i = obs.mod[i:(i+11),]
      #calc annual flow
      o.a = sum(o.i[,2:13])
      #convert to MAF
      ann.agg.flows[j,2] = o.a/1000
      
      j = j + 1
      
    }
    
    #calc n year mean unreg flows
    obs.unreg.NyrMean = rollmean(ann.agg.flows[,2], k = meanLength)
    df.unreg.NyrMean = data.frame(wyears = dt[1:length(obs.unreg.NyrMean)], q = obs.unreg.NyrMean)
    
    #subset dfs to match year range
    df.unreg.NyrMean = subset(df.unreg.NyrMean, wyears <= max(stats.Nyr$wyears))
    stats.Nyr.sub = subset(stats.Nyr, wyears >= min(df.unreg.NyrMean$wyears))
    
    ###########################################
    #generate linear regression to transform unregulated flows at LF to naturalized
    mod.data = data.frame(q.nat = stats.Nyr.sub$q.f_maf, q.unreg = df.unreg.NyrMean$q)
    
    model = lm(q.nat ~ q.unreg, data = mod.data)
    
    par(pty="s")
    title = paste0(meanLength, "-year mean flow at Lees Ferry (1980-2017)\ny = ", round(model$coefficients[1], 1), " + ", 
                   round(model$coefficients[2], 2), "x")
    plot(q.nat ~ q.unreg, data = mod.data, xlab = "Unregulated flow (maf)", ylab = "Naturalized flow (maf)", font.lab = 2,
         main = title, cex.main = 0.9, cex.lab = 0.9)
    abline(a = model$coefficients[1], b = model$coefficients[2], col = "red")
    
    cor(mod.data$q.unreg, mod.data$q.nat)
    
    ###################################################################################################################
    #process ESP forecasts 
    
    df.esp = NULL
    df.esp.skill = NULL
    
    #subset ESP hindcasts for Oct start months
    espFcsts = subset(esp$locAll_fcst, format.Date(esp$locAll_fcst$Run.Date, "%m") == 10)
    
    #extract the selected starting month (Run Date) for each N-year forecast in 1982-2017 hindcast period
    runDates = unique(espFcsts$Run.Date)
    runDates = runDates[year(runDates) %in% vs:ve]
    
    leadTime = meanLength*12
    
    #aggregate for rank histogram
    rel.mat.fcst = NULL
    rh.mat.fcst = matrix(NA, nrow = length(runDates), ncol = 40)
    rel.mat.obs = NULL
    rh.mat.obs = NULL
    
    #loop thru runDates
    i = 1
    for(i in 1:length(runDates)){
      
      runDate = runDates[i]
      startYear = as.numeric(format(runDate, '%Y'))
      
      dat.i = NULL
      mat.esp = NULL

      #select i-th RunDate (Run Dates mark the beginning of the simulation period, with ~30 traces per RunDate and a simulation length of 5 years)
      df.i = subset(espFcsts, Run.Date == runDate)
      
      #extract TraceYears - each trace year is a different ensemble member and the Trace Year marks the begininning of the 5-year period from which precip and temp were taken to generate the ESP trace
      traceYears = unlist(data.frame(yr = as.numeric(as.character(unique(df.i$TraceYear)))))
      
      #loop through each trace year (ensemble member)
      j = as.numeric(traceYears[1])
      for(j in traceYears){
        
        #subset for the j-th trace
        df.j = subset(df.i, TraceYear == j)
        
        #based on the selected N-year mean length, extract the monthly flow for this forecast period
        #convert to MAF and calc N-year mean flow from monthly values
        q.j = sum(df.j[1:(leadTime), 3:14])/1000/meanLength
        
        #convert to naturalized
        q.j.nat = as.numeric(q.j*model$coefficients[2]+model$coefficients[1])
        
        #pull in data template 
        tmp = data.frame(flow = q.j.nat)
        tmp$EnsMem = as.numeric(j)
        tmp$startYear = startYear
        tmp$tag = "ESP"
        tmp$meanLength = meanLength
        
        #calc flow terciles and assign label
        if(tmp$flow >= as.numeric(quantile(stats.Nyr$q.f_maf, 2/3))){
          tmp$flowTercile = 3
          tmp$flowCategory = "High"
        } else if(tmp$flow <= as.numeric(quantile(stats.Nyr$q.f_maf, 1/3))){
          tmp$flowTercile = 1
          tmp$flowCategory = "Low"
        }  else{
          tmp$flowTercile = 2
          tmp$flowCategory = "Average"
        }
        
        dat.i = rbind(dat.i, tmp)
        
      }
      
      #make ESP blind, then calc skill
      if(fcstMode == 'loocv'){
        esp.blind = subset(dat.i, EnsMem != startYear)
      } else if(fcstMode == 'kFold'){
        esp.blind = subset(dat.i, EnsMem <= startYear-meanLength | EnsMem >= startYear+meanLength)
      } else if(fcstMode == 'retroBlind'){
        esp.blind = subset(dat.i, EnsMem <= startYear-meanLength)
      }
      
      df.esp = rbind(df.esp, esp.blind)
      mat.esp = t(esp.blind$flow)
      mat.esp.tercs = t(esp.blind$flowTercile)
      
      #calc RPS time series (one score per year) for forecast and clim, then calc RPSS from the two scores
      if(fcstMode == 'loocv'){
          mat.clim = t(subset(stats.Nyr.tercs, wyears != startYear)[,6])
      } else if (fcstMode == 'kFold'){
          mat.clim = t(subset(stats.Nyr.tercs, wyears <= startYear-meanLength | wyears >= startYear+meanLength)[,6])
          mat.clim.cont = t(subset(stats.Nyr.tercs, wyears <= startYear-meanLength | wyears >= startYear+meanLength)[,2])
      } else if (fcstMode == 'retroBlind'){
          mat.clim = t(subset(stats.Nyr.tercs, wyears <= startYear-meanLength)[,6])
      }
      
      obs = as.numeric(subset(stats.Nyr.tercs, wyears == startYear)[2])
      obs.terc = as.numeric(subset(stats.Nyr.tercs, wyears == startYear)[6])
      
      rps.f = EnsRps(mat.esp.tercs, obs.terc)
      rps.c = EnsRps(mat.clim, obs.terc)
      rpss.f = 1 - rps.f/rps.c
      
      #calc reliability
      cd.esp = crpsDecomposition(obs, mat.esp)
      reliablity.esp = cd.esp$Reli
      
      #aggregate for rank histogram
      rh.mat.fcst[i,(1:length(mat.esp))] = mat.esp
      rh.mat.obs = rbind(rh.mat.obs, obs)
      
      #aggregate for reliability
      fcst_prob = quantInv(mat.esp, obs) # where the obs value lies in the forecast CDF
      obs_freq = quantInv(mat.clim.cont, obs) # quantile of obs flow in history
      
      rel.mat.fcst = rbind(rel.mat.fcst, fcst_prob)
      rel.mat.obs = rbind(rel.mat.obs, obs_freq)
      
      #post-processing --------------------------------------------------------
      
      terciles.i = skillScore.melt[skillScore.melt$startYear %in% startYear, 1:2]
      
      df.esp.skill.i = data.frame(flowTercile = terciles.i[1,1], stringsAsFactors=F)
      df.esp.skill.i$tempTercile = terciles.i[1,2]
      df.esp.skill.i$startYear = as.factor(startYear)
      df.esp.skill.i$variable = "RPSS"
      df.esp.skill.i$value = rpss.f
      df.esp.skill.i$model = "ESP"
      
      # df.esp.rel.i = data.frame(flowTercile = terciles.i[1,1], stringsAsFactors=F)
      # df.esp.rel.i$tempTercile = terciles.i[1,2]
      # df.esp.rel.i$startYear = as.factor(startYear)
      # df.esp.rel.i$variable = "reliability"
      # df.esp.rel.i$value = reliablity.esp
      # df.esp.rel.i$model = "ESP"

      df.esp.skill = rbind(df.esp.skill, df.esp.skill.i)

    }
  
    # #calc rank histogram
    rh.mat.fcst = rh.mat.fcst[ , apply(rh.mat.fcst, 2, function(x) !any(is.na(x)))]
    
    # rel.mat.fcst.2 = rel.mat.fcst.2[,1:which.min()]
    rh.fcst = Rankhist(ens = rh.mat.fcst, obs = rh.mat.obs)
    rh.fcst.plot = PlotRankhist(rh.fcst)
    #here I computed the ensemble percentile instead of the ensemble rank and each bin has 5% width
    Nint=100
    nt=seq(5/Nint, 1,5/Nint)
    Nint1=length(nt)
    #simul corresponds to a matrix of M ensembles (rows) for N days (columns)
    #data corresponds a vector of observation for N days
    Frequency = quantile(cumsum(rh.fcst), probs = nt)
    Frequency[2:Nint1] = Frequency[2:Nint1] - Frequency[1:(Nint1-1)]
    rankpl=Frequency#this has the frecuency for each bin
    #Di correspond to the discrepancy index see Delle Monache et al. (2006)
    DI.esp = sum( sapply( 1:Nint1, function(irank) 100*abs(Frequency[irank]/sum(Frequency)-nt[1]) ) )
    cat("\n")
    cat(paste0("ESP DI = ", DI.esp))
    cat("\n")

    #calc reliability
    cd.esp = crpsDecomposition(obs = rh.mat.obs,  eps = rh.mat.fcst)
    reliability.esp = cd.esp$Reli
    
    #reliability plot -------------------------------------------------------------------------
    # put forecast probability into bins
    XX <- probcont2disc(rel.mat.fcst, bins = thresholds) # put fcst probability into threshold bins
    fcst_bin <- XX$new # forecast probability bins
    
    ## calculate forecast probability vs observed frequency
    N.pred <- aggregate(fcst_bin, by = list(fcst_bin), length) # count forecasts per bin
    N.obs <- aggregate(rel.mat.obs, by = list(fcst_bin), sum) # sum of obs forecast frequency per bin
    obar.i = N.obs/N.pred # calc average obs frequency in each bin - should correspond to bin median
    obar.i$Group.1 = N.obs$Group.1 # add thresholds to matrix
    
    df.reli = data.frame(pred = rel.mat.fcst, obs_freq = rel.mat.obs)
    df.reli$Model = "ESP"
    
    df.reliability = rbind(df.reliability, df.reli)
    
    df.reliability.score = rbind(df.reliability.score, reliability.esp)
    
    #aggregate reliablity for later plotting
    if(meanLength <= 5){
    
      df_rel.dia = rbind.data.frame(df_rel.dia,
                                    cbind.data.frame(fcst_type = "ESP",
                                                     threshold = obar.i$Group.1, # thresholds
                                                     obar = obar.i$V1, # observed relative frequency
                                                     n.pred = N.pred$x,
                                                     label = "ESP")) # number of forecasts each bin
      
      df_sub = data.frame(x = rep(obar.i$Group.1, N.pred$x))
      df_sub$label = "ESP"
      df_sub.agg = rbind(df_sub.agg, df_sub)
    
    }
    
    skillScore.melt = rbind(skillScore.melt, df.esp.skill)
    
  }
  
  df.reliability$meanLength = meanLength
  
  df.reliability.score = as.data.frame(df.reliability.score)
  if(meanLength <= 5){
    df.reliability.score$model = c(models, "ESP")
    } else{
    df.reliability.score$model = models
    }
  df.reliability.score$meanLength = meanLength
  
  ####################################################################################################
  ####################################################################################################
  ### plotting options:
  
  #plotting 
  if(flowTerciles == T){
    
    #plot CRPSS and group by obs flow terciles
    
    df = subset(skillScore.melt, variable == "CRPSS")
    df$flowTercile = "All Years"
    df = rbind(df, subset(skillScore.melt, variable == "CRPSS" & flowTercile != "Average"))
    
    plot = ggplot(data = df, aes(x = flowTercile, y = value, fill = model)) + 
      geom_boxplot() +  theme_bw() +
      #coord_cartesian(ylim = c(min(crpss.melt$CRPSS), 1)) +
      geom_hline(yintercept = 0) + xlab(paste0("Flow terciles of observed ", meanLength, "-year means")) +
      #ggtitle(paste0("CRPSS for average-, high-, and low- flow years\n", vl, " different 5-year blocks")) +
      theme(axis.title = element_text(face = "bold")) +
      ylab("CRPSS")
    
    print(plot)
    
    plotPath = "K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Code/plots/midtermForecasts_multiYearMeans/"
    ggsave(paste0("crpss_flowTerciles_", meanLength, "yr_fcstMode=", fcstMode, ".png"),  plot = print(plot), path = plotPath, device = "png", width = 10, height = 5, dpi = 450, units = "in")
    
    
    #plot RPSS and group by obs flow terciles
    
    df = subset(skillScore.melt, variable == "RPSS")
    df$flowTercile = "All Years"
    df = rbind(df, subset(skillScore.melt, variable == "RPSS" & flowTercile != "Average"))
    
    plot = ggplot(data = df, aes(x = flowTercile, y = value, fill = model)) + 
      geom_boxplot() +  theme_bw() +
      #coord_cartesian(ylim = c(min(crpss.melt$CRPSS), 1)) +
      geom_hline(yintercept = 0) + xlab(paste0("Flow terciles of observed ", meanLength, "-year means")) +
      #ggtitle(paste0("CRPSS for average-, high-, and low- flow years\n", vl, " different 5-year blocks")) +
      theme(axis.title = element_text(face = "bold")) +
      ylab("RPSS")
    
    print(plot)
    
    ggsave(paste0("rpss_flowTerciles_", meanLength, "yr_fcstMode=", fcstMode, ".png"), plot = print(plot), path = plotPath, device = "png", width = 10, height = 5, dpi = 450, units = "in")
    
  }
  
  # if(flowTerciles == F & tempTerciles == F){
  #   print(ggplot(data = crpss.melt, aes(x = covariate, y = CRPSS)) + 
  #   geom_boxplot(width = 0.3) + 
  #   coord_cartesian(ylim = c(min(crpss.melt$CRPSS), 1)) +
  #   geom_hline(yintercept = 0)+
  #   theme(axis.title.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank()) +
  #   #ggtitle(paste0("CRPSS for all 5-year blocks\n", vl, " different 5-year blocks")) +
  #   theme(axis.title = element_text(face = "bold"))
  #   )
  # }
  
  if(tempTerciles == T){
    
    #plot CRPSS and group by obs flow terciles
    
    df = subset(skillScore.melt, variable == "CRPSS")
    df$tempTercile = "All Years"
    df = rbind(df, subset(skillScore.melt, variable == "CRPSS" & tempTercile != "Average"))
    
    plot = ggplot(data = df, aes(x = tempTercile, y = value, fill = model)) + 
      geom_boxplot() +  theme_bw() +
      #coord_cartesian(ylim = c(min(crpss.melt$CRPSS), 1)) +
      geom_hline(yintercept = 0) + xlab(paste0("Observed ", meanLength, "-year mean temperature terciles")) +
      #ggtitle(paste0("CRPSS for average-, high-, and low- flow years\n", vl, " different 5-year blocks")) +
      theme(axis.title = element_text(face = "bold")) +
      ylab("CRPSS")
    
    print(plot)
    
    ggsave(paste0("crpss_tempTerciles_", meanLength, "yr_fcstMode=", fcstMode, ".png"), plot = print(plot), path = plotPath, device = "png", width = 10, height = 5, dpi = 450, units = "in")
    
    
    #plot RPSS and group by obs flow terciles
    
    df = subset(skillScore.melt, variable == "RPSS")
    df$tempTercile = "All Years"
    df = rbind(df, subset(skillScore.melt, variable == "RPSS" & tempTercile != "Average"))
    
    plot = ggplot(data = df, aes(x = tempTercile, y = value, fill = model)) + 
      geom_boxplot() +  theme_bw() +
      #coord_cartesian(ylim = c(min(crpss.melt$CRPSS), 1)) +
      geom_hline(yintercept = 0) + xlab(paste0("Observed ", meanLength, "-year mean temperature terciles")) +
      #ggtitle(paste0("CRPSS for average-, high-, and low- flow years\n", vl, " different 5-year blocks")) +
      theme(axis.title = element_text(face = "bold")) +
      ylab("RPSS")
    
    print(plot)
    
    ggsave(paste0("rpss_tempTerciles_", meanLength, "yr_fcstMode=", fcstMode, ".png"), plot = print(plot), path = plotPath, device = "png", width = 10, height = 5, dpi = 450, units = "in")
    
    
  }
  
  
  if(plotByYear == T){
    
    skillScore.melt$startYear = as.numeric(as.character(skillScore.melt$startYear))
  
    p1a = ggplot(data = subset(skillScore.melt, variable == "CRPSS")) + 
      geom_line(mapping = aes(x = startYear, y = value, linetype = model)) +
      geom_point(mapping = aes(x = startYear, y = value, color = flowTercile)) +
      theme_bw() +
      #geom_boxplot(aes(fill = covariate)) + 
      geom_hline(yintercept = 0) + xlab(paste0("Starting year of projected ", meanLength, "-year mean flow")) +
      theme(axis.title = element_text(face = "bold"), axis.text.x = element_text()) +
      ylab("CRPSS") 
    #scale_x_discrete(breaks=seq(min(as.numeric(as.character(crpss.melt$startYear))), max(as.numeric(as.character(crpss.melt$startYear))), 10))
    
    print(p1a)
    
    ggsave(paste0("crpss_skillTS_", meanLength, "yr_fcstMode=", fcstMode, ".png"), plot = print(p1a), path = plotPath, device = "png", width = 10, height = 5, dpi = 450, units = "in")
    
    #plot RPSS by year
    p1b = ggplot(data = subset(skillScore.melt, variable == "RPSS")) + 
      geom_line(mapping = aes(x = startYear, y = value, linetype = model)) +
      geom_point(mapping = aes(x = startYear, y = value, color = flowTercile)) +
      theme_bw() +
      #geom_boxplot(aes(fill = covariate)) + 
      geom_hline(yintercept = 0) + xlab(paste0("Starting year of projected ", meanLength, "-year mean flow")) +
      theme(axis.title = element_text(face = "bold"), axis.text.x = element_text()) +
      ylab("RPSS")
    #scale_x_discrete(breaks=seq(min(as.numeric(as.character(crpss.melt$startYear))), max(as.numeric(as.character(crpss.melt$startYear))), 10))
    
    print(p1b)
    
    ggsave(paste0("rpss_skillTS_", meanLength, "yr_fcstMode=", fcstMode, ".png"), plot = print(p1b), path = plotPath, device = "png", width = 10, height = 5, dpi = 450, units = "in")
    
  }
  
  #plot reliablity for all
  if(meanLength <= 5){
    subtitle = paste0("ESP reliability score = ", round(subset(df.reliability.score, model == "ESP")[1],2), 
                      "\nRF reliability score = ", round(subset(df.reliability.score, model == "CESM-DPLE, RF")[1],2))
  } else{
    subtitle = paste0("RF reliability score = ", round(subset(df.reliability.score, model == "CESM-DPLE, RF")[1],2))
  }
  
  p1 <- ggplot(subset(df.reliability, Model != "CESM-DPLE, KNN"), aes(pred, obs_freq, color = Model)) +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey20', size = 0.5) +
    geom_vline(xintercept = 0.2, linetype = 'dashed', colour = 'grey20', size = 0.5) +
    geom_vline(xintercept = 0.4, linetype = 'dashed', colour = 'grey20', size = 0.5) +
    geom_vline(xintercept = 0.6, linetype = 'dashed', colour = 'grey20', size = 0.5) +
    geom_vline(xintercept = 0.8, linetype = 'dashed', colour = 'grey20', size = 0.5) +
    geom_vline(xintercept = 1, linetype = 'dashed', colour = 'grey20', size = 0.5) +
    geom_point() +
    theme_bw() +
    xlab('Forecast Probability') +
    ylab('Observed Frequency') +
    scale_x_continuous(expand = c(0.0,0.0), breaks = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(-0.02,1.02))+
    scale_y_continuous(expand = c(0.0,0.0), breaks = seq(0,1,0.2), limits = c(-0.02,1.02)) +
    theme(panel.grid.minor = element_line(size = 0.02), 
          panel.grid.major = element_line(color = 'grey65', size = 0.5),
          axis.title = element_text(face = "bold")) +
    coord_fixed() + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggtitle(paste0("Mean Length = ", meanLength, "-years"),
            subtitle = subtitle) +
    labs(color = "Model")
  # gg_circle(r=0.02, xc=fcst_prob, yc=obs_freq, color = 'red')
  
  print(p1)
  
  plotPath = "K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Code/plots/midtermForecasts_multiYearMeans/"
  title = paste0("reliability.", meanLength, "yr_fcstMode=", fcstMode, ".png")
  ggsave(title, plot = print(p1), device = "png", path = plotPath, width = 6, height = 6, dpi = 450, units = "in")
  
  
  ## small histogram for plots (subplot)
  sub = ggplot(subset(df_sub.agg, label != "CESM-DPLE, KNN"), aes(x = x,  fill = label)) + 
    geom_histogram(col = 'black', binwidth = 0.1, position = 'dodge') +
    theme_classic() +
    theme(legend.position="none", 
          axis.title.x=element_blank(), axis.text.x=element_blank(), 
          axis.title.y=element_blank(), axis.text.y=element_blank()) +
    scale_x_continuous(breaks = seq(0,10,2), expand = c(0.0,0.0)) +
    scale_y_continuous(expand = c(0.0,0.0)) + 
    theme(panel.grid = element_blank(), panel.border = element_blank()) 
  
  sub
  
  ## reliability plots containing subplots
  rel.plot = ggplot(subset(df_rel.dia, label != "CESM-DPLE, KNN"), aes(x = threshold, y = obar, color = label)) +
    geom_abline(slope = 1, intercept = 0, col = 'black', linetype = "dashed")+
    geom_line() +
    geom_point() +
    scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0,1), expand = c(0,0)) +
    scale_x_continuous(breaks=seq(0, 1, 0.2), limits = c(0,1), expand = c(0,0)) +
    theme_bw() +
    # theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank()) +
    annotation_custom(ggplotGrob(sub), xmin = 0.01, xmax = 0.27, ymin = 0.63, ymax = 0.99) +
    xlab('Forecast Probability') +
    ylab('Observed Frequency') +
    ggtitle(paste0("Mean Length = ", meanLength, "-years"),
            subtitle = subtitle) +
    labs(color = "Model")+
    theme(axis.title = element_text(face = "bold")) 
    # facet_wrap(~label, nrow = 1, scale = 'fixed') #, scales = "free_y"
  
  rel.plot
  
  title = paste0("reliability_Traditional.", meanLength, "yr_fcstMode=", fcstMode, ".png")
  ggsave(title, plot = print(rel.plot), device = "png", path = plotPath, dpi = 450, units = "in")
  
  
  #################################################################################
  skillScore.melt$startYear = as.numeric(as.character(skillScore.melt$startYear))
  
  
  
  return(list(skillScore.melt, df.reliability, df.reliability.score, df_rel.dia, df_sub.agg))
  
}

#####################################################################################################################
#####################################################################################################################

############################
### crpss on annual flows

bootstrap = T
flowTerciles = T

crpss.calc.ann = function(sim.results, bootstrap, flowTerciles, vl, vs, stats, meanLength){
  
  #initialize list to store crpss from different models, currently only 1 model in use
  models = unique(sim.results$AnnualFlowBootstrap$covariate)
  
  #crpss.list = vector("list", length(models))
  crpss.list = NULL
  
  #loop thru different models
  i = 1
  for(i in 1:length(models)){
    
    #select bootstrapped forecasts vs non-bootstrapped (neighbors only)
    if(bootstrap == T){fcsts.i = sim.results$AnnualFlowBootstrap}
    if(bootstrap == F){fcsts.i = sim.results$AnnualFlow}
    
    #subset forCecasts to i-th model of loop
    fcsts.i = subset(fcsts.i, covariate == models[i])
    
    #initialize matrix to store CRPSS and obs flow terciles
    crpss.cov = matrix(NA, nrow = vl*meanLength, ncol = meanLength)
    
    #loop through each year in validation (hindcast) period
    l = 1
    j = 1
    for(j in 1:vl){
      
      #loop thru fcst years 1:5
      k = 1
      for(k in 1:meanLength){
        
        #subset forecasts to j-th start year of validation period
        fcsts.j = subset(fcsts.i, startYear == vs+j-1 & fcstYear == k)
        
        #fetch annual flow climatology from observed record
        clim = t(matrix(subset(stats, wyears != vs+j-1)[,3]))
        #clim = t(matrix(subset(stats, wyears < vs+j-1)[,3]))
        
        #re-format ensemble data
        ens = t(matrix(fcsts.j$q.maf))
        
        #pull observed annual flow
        obs = t(matrix(fcsts.j$obsFlow[1]))
        
        #calc RPSS on forecast and climatology
        rpss.f = EnsCrps(ens = ens, obs = obs)
        rpss.c = EnsCrps(ens = clim, obs = obs)
        
        #calc CRPSS and store
        crpss = 1 - rpss.f/rpss.c
        
        crpss.cov[l+k-1, 1] = crpss
        crpss.cov[l+k-1, 2] = k 
        
        #calc flow terciles and assign label
        if(fcsts.j$obsFlow[1] >= as.numeric(quantile(stats$q_maf, 2/3))){
          crpss.cov[l+k-1,3] = "High"
        } else if(fcsts.j$obsFlow[1] <= as.numeric(quantile(stats$q_maf, 1/3))){
          crpss.cov[l+k-1,3] = "Low"
        }  else{
          crpss.cov[l+k-1,3] = "Average"
        }
        
        crpss.cov[l+k-1,4] = fcsts.j$simYear[1]
        
        crpss.cov[l+k-1,5] = models[i]
        
      }
      
      l = l + meanLength
      
    }
    
    colnames(crpss.cov) = c("CRPSS", "Year", "FlowTercile", "SimYear", "Model")
    
    crpss.list = rbind(crpss.list, crpss.cov)
    crpss.list = as.data.frame(crpss.list)
    
    crpss.list$CRPSS = as.numeric(as.character(crpss.list$CRPSS))
  }
  
  
  #plot
  print(
    ggplot(data = crpss.list, aes(x = Year, y = CRPSS, fill = Model)) + 
      geom_boxplot() + theme_bw() + 
      #geom_boxplot(aes(fill = Covariate)) + 
      coord_cartesian(ylim = c(min(crpss.list$CRPSS), 1)) +
      geom_hline(yintercept = 0) + 
      #ggtitle(paste0("CRPSS by forecast year\n", vl, " different 5-year blocks")) +
      theme(axis.title = element_text(face = "bold"))
    
  )
  
  
  #plotting 
  if(flowTerciles == T){
    
    #plot CRPSS and group by obs flow terciles
    
    print(
      ggplot(data = crpss.list, aes(x = Year, y = CRPSS, fill = FlowTercile)) + 
        geom_boxplot() + theme_bw() + 
        #geom_boxplot(aes(fill = Covariate)) + 
        coord_cartesian(ylim = c(min(crpss.list$CRPSS), 1)) +
        geom_hline(yintercept = 0) + 
        #ggtitle(paste0("CRPSS by forecast year\n", vl, " different 5-year blocks")) +
        theme(axis.title = element_text(face = "bold")) +
        facet_wrap(~Model, scales = "free", shrink = F)
      
      
    )
    
    
  }
  
  return(crpss.list)
  
}

#####################################################################################################################
#####################################################################################################################

############################
### rpss on annual flows

bootstrap = T
flowTerciles = T


rpss.calc.ann = function(sim.results, bootstrap, flowTerciles, vl, vs, stats){
  
  #initialize list to store crpss from different models, currently only 1 model in use
  models = unique(sim.results$AnnualFlowBootstrap$covariate)
  
  #crpss.list = vector("list", length(models))
  rpss.list = NULL
  
  #loop thru different models
  i = 1
  for(i in 1:length(models)){
    
    #select bootstrapped forecasts vs non-bootstrapped (neighbors only)
    if(bootstrap == T){fcsts.i = sim.results$AnnualFlowBootstrap}
    if(bootstrap == F){fcsts.i = sim.results$AnnualFlow}
    
    #subset forCecasts to i-th model of loop
    fcsts.i = subset(fcsts.i, covariate == models[i])
    
    #initialize matrix to store CRPSS and obs flow terciles
    rpss.cov = matrix(NA, nrow = vl*meanLength, ncol = meanLength)
    
    #loop through each year in validation (hindcast) period
    l = 1
    j = 1
    for(j in 1:vl){
      
      #loop thru fcst years 1:5
      k = 1
      for(k in 1:meanLength){
        
        #subset forecasts to j-th start year of validation period
        fcsts.j = subset(fcsts.i, startYear == vs+j-1 & fcstYear == k)
        
        #calc simulated flow terciles and assign label
        m = 1
        for(m in 1:nrow(fcsts.j)){
          
          if(fcsts.j$q.maf[m] >= as.numeric(quantile(stats$q_maf, 2/3))){
            fcsts.j$simFlowTercile[m] = "High"
            fcsts.j$simFlowTercileNo[m] = 3
          } else if(fcsts.j$q.maf[m] <= as.numeric(quantile(stats$q_maf, 1/3))){
            fcsts.j$simFlowTercile[m] = "Low"
            fcsts.j$simFlowTercileNo[m] = 1
          }  else{
            fcsts.j$simFlowTercile[m] = "Average"
            fcsts.j$simFlowTercileNo[m] = 2
          }
          
          
        }
        
        
        #calc obs flow terciles and assign label
        if(fcsts.j$obsFlow[1] >= as.numeric(quantile(stats$q_maf, 2/3))){
          fcsts.j$obsFlowTercile = "High"
          fcsts.j$obsFlowTercileNo = 3
        } else if(fcsts.j$obsFlow[1] <= as.numeric(quantile(stats$q_maf, 1/3))){
          fcsts.j$obsFlowTercile = "Low"
          fcsts.j$obsFlowTercileNo = 1
        }  else{
          fcsts.j$obsFlowTercile = "Average"
          fcsts.j$obsFlowTercileNo = 2
        }
        
        #fetch annual flow climatology from observed record
        clim = t(matrix(subset(stats, wyears != vs+j-1)[,13]))
        
        #re-format ensemble data
        ens = t(matrix(fcsts.j$simFlowTercileNo))
        
        #pull observed annual flow
        obs = t(matrix(fcsts.j$obsFlowTercileNo[1]))
        
        #calc RPSS on forecast and climatology
        rpss.f = EnsRps(ens = ens, obs = obs)
        rpss.c = EnsRps(ens = clim, obs = obs)
        
        #calc CRPSS and store
        rpss = 1 - rpss.f/rpss.c
        
        rpss.cov[l+k-1,1] = rpss
        rpss.cov[l+k-1,2] = k 
        rpss.cov[l+k-1,3] = fcsts.j$obsFlowTercile[1]
        rpss.cov[l+k-1,4] = fcsts.j$simYear[1]
        rpss.cov[l+k-1,5] = models[i]
        
      }
      
      l = l + meanLength
      
    }
    
    colnames(rpss.cov) = c("RPSS", "Year", "FlowTercile", "SimYear", "Model")
    
    rpss.list = rbind(rpss.list, rpss.cov)
    rpss.list = as.data.frame(rpss.list)
    
    rpss.list$RPSS = as.numeric(as.character(rpss.list$RPSS))
  }
  
  
  #plot
  print(
    ggplot(data = rpss.list, aes(x = Year, y = RPSS, color = Model)) + 
      geom_boxplot() + theme_bw() + 
      #geom_boxplot(aes(fill = Covariate)) + 
      coord_cartesian(ylim = c(min(rpss.list$RPSS), 1)) +
      geom_hline(yintercept = 0) + 
      #ggtitle(paste0("CRPSS by forecast year\n", vl, " different 5-year blocks")) +
      theme(axis.title = element_text(face = "bold"))
    
  )
  
  
  #plotting 
  if(flowTerciles == T){
    
    #plot CRPSS and group by obs flow terciles
    
    print(
      ggplot(data = rpss.list, aes(x = Year, y = RPSS, color = FlowTercile)) + 
        geom_boxplot() + theme_bw() + 
        #geom_boxplot(aes(fill = Covariate)) + 
        coord_cartesian(ylim = c(min(rpss.list$RPSS), 1)) +
        geom_hline(yintercept = 0) + 
        #ggtitle(paste0("CRPSS by forecast year\n", vl, " different 5-year blocks")) +
        theme(axis.title = element_text(face = "bold")) +
        facet_wrap(~Model, scales = "free", shrink = F)
      
      
    )
    
    
  }
  
  return(rpss.list)
  
}


#####################################################################################################################
#####################################################################################################################

#####################################################################################################################
#####################################################################################################################

### Main Function code (read in obs, climate data, pre-process, run simulations, calculate skill)
#input with varied lead times
# fcstMode = 'kFold'
# meanLength = 2

midTermProjection = function(meanLength, fcstMode){
  
  ###########################################################################################################################
  ### DATA PROCESSING: OBSERVED DATA AND CESM-DPLE DATA (5 YEAR MEANS)
  ###########################################################################################################################

  # setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/txt_files")
  # write.table(fwdNyr_means, file = "fwdNyr_means.txt", row.names = F)
   
  ### Read in and process CESM LE or LME MONTHLY data
   
  setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/txt_files/CESM")
  #setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/CESM/CESM_LME/monthly")
  
  #DP
  dp.dat = paste0("dp.le.", meanLength, "yr_means.txt")
  dp.le = read.csv(file = dp.dat, sep = " ")
  
  #calc ensemble mean
  DPLE.pcp = apply(dp.le[,2:6], 1, mean)
  DPLE.Tmax = apply(dp.le[,7:11], 1, mean)
  DPLE.Tmin = apply(dp.le[,12:16], 1, mean)
  
  #all quantiles
  dp = data.frame(wyears = dp.le$wyears, DPLE.Tmax, DPLE.Tmin, DPLE.pcp)

  ### Read in and process CESM LE or LME MONTHLY data
  setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/CESM/CESM_LE/post1920/monthly")
  #setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/CESM/CESM_LME/monthly")
  
  #Tmax
  Tmax.le = read.csv(file = "TREFHTMX_quantile_UCO_1920-2080.txt", sep = "") 
  #Tmax.le = read.csv(file = "TREFHTMX_quantile_anom_UCO_1920-2080.txt", sep = "") 
  #Tmax.le = read.csv(file = "TREFHTMX_quantile_UCO_LME.txt", sep = ",") 
  yearmo = Tmax.le$yearmoth
  
  #Tmin
  Tmin.le = read.csv(file = "TREFHTMN_quantile_UCO_1920-2080.txt", sep = "") 
  #Tmin.le = read.csv(file = "TREFHTMN_quantile_anom_UCO_1920-2080.txt", sep = "") 
  #Tmin.le = read.csv(file = "TREFHTMN_quantile_UCO_LME.txt", sep = ",") 
  
  #precip
  pcp.le = read.csv(file = "PRECT_quantile_UCO_1920-2080.txt", sep = "")
  #pcp.le = read.csv(file = "PRECT_quantile_anom_UCO_1920-2080.txt", sep = "")
  
  #calc ensemble mean
  ensMean.Tmax.C = apply(Tmax.le[,2:6], 1, mean) - 273.15
  ensMean.Tmin.C = apply(Tmin.le[,2:6], 1, mean) - 273.15
  ensMean.pcp.mm = apply(pcp.le[,2:6], 1, mean)
  
  #all quantiles
  LE.monthly = cbind(yearmo, ensMean.Tmax.C, ensMean.Tmin.C, pcp.le$med)
  #colnames(DPLE.seas)[1] = "yearmo"
  

  ### calc annual (water year) covariates from CESM-DPLE
   
  #remove first 8 rows of DPLE.seas to prepare for water year format and remove all rows after 082017 to match end of record for obs flows
  LE.monthly.sub = subset(LE.monthly, yearmo >= 192010 & yearmo <= 201709)
  #DPLE.seas = subset(DPLE.seas, DPLE.seas$wyears >= 85009 & DPLE.seas$wyears <= 201708)
  wyears = 1921:2017
  
  #seq = seq(1, nrow(DPLE.seas), 3)
  #seq = 1:(nrow(DPLE.seas)-11)
  seq = seq(1, nrow(LE.monthly.sub)-11, 12)
  
  #initialize dataframes to store annual precip, temperature
  LE.Tmax = NULL
  LE.Tmin = NULL
  LE.pcp = NULL
  
  i = 1
  for(i in seq){
    
    ann.i.Tmax =  mean(LE.monthly.sub[i:(i+11),2])
    ann.i.Tmin =  mean(LE.monthly.sub[i:(i+11),3])
    ann.i.pcp =  sum(LE.monthly.sub[i:(i+11),4])
    
    LE.Tmax = rbind(LE.Tmax, ann.i.Tmax)
    LE.Tmin = rbind(LE.Tmin, ann.i.Tmin)
    LE.pcp = rbind(LE.pcp, ann.i.pcp)
    
  }
  
  LE.ann = data.frame(wyears, LE.Tmax, LE.Tmin, LE.pcp)
  rownames(LE.ann) = NULL
  
  ### Read in UCRB stats and naturalized flow at Lees Ferry
  stats = read.csv(file = "K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/txt_files/UCRB_annual_stats.txt", sep = " ")
  #stats = read.csv(file = "meko.comb.flows.txt", sep = ",")
  
  clim.ann = mean(subset(stats, wyears >= 1981 & wyears <= 2010)[,2])

  setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Code/plots/observedMultiYearMeanRelationships")

  ### plot obs annual flows
  plot = ggplot(data = stats, aes(x = wyears, y = q_maf)) +
    geom_line() + geom_point() + theme_bw() + 
    xlab("Water Year") + ylab("Flow (MAF)") +
    geom_smooth(method="loess", se=F, fullrange=T, level=0.95, col = "red") +
    theme(axis.title = element_text(face = "bold"))
  plot 
  
  ggsave("obsFlowTS.png", plot = print(plot), device = "png", width = 10, height = 5, dpi = 450, units = "in")
    
  ### calc an annual Tmean for later use
  # mean.temp.ann = data.frame(wyears = DPLE.ann$wyears, t.C = apply(DPLE.ann[,2:3], 1, mean))
  # mean.temp.Nyr = data.frame(wyears = head(mean.temp.ann$wyears, -4), t.C = rollmean(mean.temp.ann$t.C, 5))
  mean.temp.Nyr = data.frame(wyears = head(stats$wyears, -(meanLength-1)), t.C = rollmean(stats$basin_tmean_C, meanLength))

  ### Correlation search (past-n year mean flow correlation with future 5-year mean flow)
  #loop to find optimal correlations
  corSearch = NULL
  i = 2
  for(i in 2:70){
    
    stats.Nyr.corSearch = flow.avg(stats, i, meanLength)
    
    cor.i = cor(stats.Nyr.corSearch$q.p, stats.Nyr.corSearch$q.f_maf)
    
    df.i = cbind(i, cor.i)
    
    corSearch = rbind(corSearch, df.i)
    
  }
  
  corSearch.qp = data.frame(corSearch)
  colnames(corSearch.qp) = c("Years", "Correlation")
  title = "Streamflow"
  #title = "Future 5-year mean flow correlation with different window lengths for past flow averaging"
  plot(Correlation ~ Years, corSearch.qp, type = 'l', xlab = "Past flow averaging window length (years)", main = title, ylab = paste0("Correlation with future", meanLength, "-year mean flow"), cex.main = 1, font.lab = 2)
  abline(h = 0)
  
  maxCor.qp = top_n(corSearch.qp, 1, Correlation)
  minCor.qp = top_n(corSearch.qp, -1, Correlation)
  
  ################################################################################
  #plot optimal correlations
  title = paste0("Future ", meanLength, "-year mean flow")
  
  #######################
  #max cor
  i = maxCor.qp$Years
  stats.Nyr.corSearch = flow.avg(stats, i, meanLength)
  par(pty = 's')
  plot(stats.Nyr.corSearch$q.p, stats.Nyr.corSearch$q.f_maf, xlab = paste0("Past ", i, "-year mean flow"), ylab = title)
  
  #######################
  #min cor
  i = minCor.qp$Years
  stats.Nyr.corSearch = flow.avg(stats, i, meanLength)
  par(pty = 's')
  plot(stats.Nyr.corSearch$q.p, stats.Nyr.corSearch$q.f_maf, xlab = paste0("Past ", i, "-year mean flow"), ylab = title)
  
  ### Correlation search (past-n year runoff efficiency correlation with future 5-year mean flow)
  #loop to find optimal correlations
  corSearch = NULL
  i = 2
  for(i in 2:70){
    
    stats.Nyr.corSearch = flowRE.avg(stats, i, meanLength, LE.ann)
    
    cor.i = cor(stats.Nyr.corSearch$re.p, stats.Nyr.corSearch$q.f_maf)
    
    df.i = cbind(i, cor.i)
    
    corSearch = rbind(corSearch, df.i)
    
  }
  
  corSearch.re = data.frame(corSearch)
  colnames(corSearch.re) = c("Years", "Correlation")
  title = "Runoff efficiency"
  #title = "Future 5-year mean flow correlation with different window lengths for past runoff efficiency averaging"
  plot(Correlation ~ Years, corSearch.re, type = 'l', xlab = "Past runoff efficiency averaging window length (years)", main = title, ylab = paste0("Correlation with future ", meanLength, "-year mean flow"), cex.main = 1, font.lab = 2)
  abline(h = 0)
  
  maxCor.re = top_n(corSearch.re, 1, Correlation)
  minCor.re = top_n(corSearch.re, -1, Correlation)
  
  ################################################################################
  #plot optimal correlations
  title = paste0("Future ", meanLength, "-year mean flow")
  
  #######################
  #max cor
  i = maxCor.re$Years
  stats.Nyr.corSearch = flowRE.avg(stats, i, meanLength, LE.ann)
  par(pty = 's')
  plot(stats.Nyr.corSearch$re.p, stats.Nyr.corSearch$q.f_maf, xlab = paste0("Past ", i, "-year mean RE"), ylab = title)
  cor(stats.Nyr.corSearch$re.p, stats.Nyr.corSearch$q.f_maf)
  
  #######################
  #min cor
  i = minCor.re$Years
  stats.Nyr.corSearch = flowRE.avg(stats, i, meanLength, LE.ann)
  par(pty = 's')
  plot(stats.Nyr.corSearch$re.p, stats.Nyr.corSearch$q.f_maf, xlab = paste0("Past ", i, "-year mean RE"), ylab = title)
  cor(stats.Nyr.corSearch$re.p, stats.Nyr.corSearch$q.f_maf)
  
  
  #######################
  #custom cor
  i = corSearch.re[which.max(abs(corSearch.re$Correlation[1:14])),1]
  stats.Nyr.corSearch = flowRE.avg(stats, i, meanLength, LE.ann)
  par(pty = 's')
  plot(stats.Nyr.corSearch$re.p, stats.Nyr.corSearch$q.f_maf, xlab = paste0("Past ", i, "-year mean RE"), ylab = title)
  cor(stats.Nyr.corSearch$re.p, stats.Nyr.corSearch$q.f_maf)
  

  ### plot both past flow and past RE on same plot
  #setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Code/plots")
  
  png(paste0("autocorAnalysis_", meanLength, "yr_meanFlow.png"), width = 1920, height = 1080,
      units = "px", bg = "white", pointsize = 24, res = NA,
      restoreConsole = TRUE)
  
  colnames(corSearch.re) = c("Years", "Correlation")
  title = "Runoff efficiency"
  #title = "Future 5-year mean flow correlation with different window lengths for past runoff efficiency averaging"
  plot(Correlation ~ Years, corSearch.re, type = 'l', xlab = "Averaging window length of past data (years)", ylab = paste0("Correlation with future ", meanLength, "-year mean flow"), cex.main = 1,
       font.lab = 2, ylim = range(corSearch.qp$Correlation, corSearch.re$Correlation),  col = "red")
  lines(x = corSearch.qp$Years, y = corSearch.qp$Correlation, lty = 2)
  abline(h = 0)
  legend("bottomleft", c("Past flow", "Past RE"), col = c("black", "red"), lty = c(2,1), bty = "n")
  
  dev.off()
  
  ### Calc 5-year mean obs flow at Lees Ferry, pair with optimal past mean flow and past mean runoff efficiency covariates
   
  ###########################################
  #past flow as covariate
  #use window length with optimal correlation
  
  i = 2
  #i = maxCor.qp$Years
  window.length.qp = i
  #window.length.qp = maxCor.qp$Years
  
  stats.Nyr.qp = flow.avg(stats, i, meanLength)
  
  ###########################################
  #past runoff efficiency as covariate
  #use window length with optimal correlation
  
  i = corSearch.re[which.max(abs(corSearch.re$Correlation[1:14])),1]
  #i = 30
  #i = minCor.re$Years
  window.length.re = i
  #window.length.re = minCor.re$Years
  stats.sub = subset(stats, wyears >= min(LE.ann$wyears))  
  
  stats.Nyr.rep = flowRE.avg(stats.sub, i, meanLength, LE.ann)
  
  ############################################
  #merge stats.Nyr... dataframes
  
  stats.Nyr = inner_join(stats.Nyr.qp, stats.Nyr.rep, by = "wyears")
  stats.Nyr = stats.Nyr[,-4]
  colnames(stats.Nyr)[2] = "q.f_maf"
  
  ###############
  # #use past 5-year mean flow
  # window.length = 5
  #
  # #add wyear col
  # stats.Nyr = data.frame(wyears = stats$wyears[6:(nrow(stats)-4)])
  # 
  # #create future 5 year means - remove first 5 years to match with later calculation of past 5 year mean
  # stats.Nyr$q.f_maf = rollmean(stats$q_maf[-(1:5)], 5, align = "left")
  # 
  # #past 5 year mean (mean of scaled anomalies) - remove last 5 years of flow
  # #stats.Nyr$q.p = rollmean(head(stats$q_maf, -5), k = 5)
  # stats.Nyr$q.p = rollmean(scale(head(stats$q_maf, - 5)), k = 5)
  
  clim.Nyr = mean(subset(stats.Nyr, wyears >= 1951 & wyears <= 2017)[,2])
  

  ### create 5 year mean anomalies to be used as model covariates
  #CESM-DPLE annual data (DPLE.ann) is used here, but seasonal LE data could also be used
  DPLE.Nyr = mean.anoms(dp, meanLength, stats.Nyr, window.length.re)
  
  p = ggarrange(DPLE.Nyr$plots[[1]], DPLE.Nyr$plots[[2]], DPLE.Nyr$plots[[4]], ncol = DPLE.Nyr$NumCovars-1, nrow = 1, common.legend = T, legend = "right")
  p2 = annotate_figure(p, top = text_grob(paste0("CESM-DPLE correlation with obs. flow (", meanLength,"-year means), ", DPLE.Nyr$meanAnoms$wyears[1], "-", tail(DPLE.Nyr$meanAnoms$wyears, 1)), face = "bold", size = 14, vjust = 0.5))
  
  p2
  
  title = paste0("flowVsDPLEmet_post1981_noPCP_", meanLength, "yrMeans.png")
  ggsave(title, plot = print(p2), device = "png", width = 15, height = 5, dpi = 450, units = "in")
  
  dat = DPLE.Nyr$meanAnoms
  
  dat$DPLE.Tmean = 0.4*dat$DPLE.Tmin + 0.6*dat$DPLE.Tmax
  
  plot = ggplot(data = dat, aes(x = DPLE.Tmean, y = DPLE.pcp, color = q.f_maf)) +
    geom_point() + 
    xlab("Mean Temp. Scaled Anomaly") + ylab("Precip. Scaled Anomaly") + theme_bw() + theme(aspect.ratio=1) +
    scale_colour_gradientn(colours = rev(hcl.colors(40, palette = "viridis"))) + labs(colour = paste0(meanLength, "-year Mean Flow\n(MAF)")) +
    geom_smooth(method='loess', se=TRUE, fullrange=T, level=0.95, col = "black") +
    #geom_smooth(method="auto", se=TRUE, fullrange=T, level=0.95, col = "black") +
    #ggtitle(paste0("Mutual Information = ", round(meanAnomsFlowCor[i-1,1], digits = 2))) +
    theme(axis.title = element_text(face = "bold"))
  
  plot
   
  
  # scatter plot
  lim = c(-2,2)
  nclr <- 15 # number of colors
  plotclr <- rev(viridis(nclr)) # get the colors
  colornum <- cut(rank(dat$q.f_maf), nclr, labels=FALSE)
  colcode <- plotclr[colornum] # assign color
  plot.angle <- 45
  scatterplot3d(x = dat$DPLE.Tmean, y = dat$DPLE.pcp, z = dat$q.f_maf, type="h", angle=plot.angle, color=colcode, pch=20, cex.symbols=1,
                box = F, xlab = "Mean Temp. Scaled Anomaly", ylab = "", zlab = paste0(meanLength, "-year Mean Flow (MAF)"), font.axis = 2, main = "CESM-DPLE", xlim = lim, ylim = lim)
  text(x = 3.5, y = 4.55, "Precip. Scaled Anomaly", srt = 22, font = 2)
   
  data.loess <- loess(q.f_maf ~ DPLE.pcp * DPLE.Tmean, data = dat)
  
  # Create a sequence of incrementally increasing (by 0.3 units) values for both wt and hp
  xgrid <-  seq(min(dat$DPLE.pcp)-0.1, max(dat$DPLE.pcp)+0.1, 0.01)
  ygrid <-  seq(min(dat$DPLE.Tmean)-0.1, max(dat$DPLE.Tmean)+0.1, 0.01)
  # Generate a dataframe with every possible combination of wt and hp
  data.fit <-  expand.grid(DPLE.pcp = xgrid, DPLE.Tmean = ygrid)
  # Feed the dataframe into the loess model and receive a matrix output with estimates of
  # acceleration for each combination of wt and hp
  mtrx3d <-  predict(data.loess, newdata = data.fit)
  # Abbreviated display of final matrix
  mtrx3d[1:4, 1:4]
  
  filled.contour(x = xgrid, y = ygrid, z = mtrx3d,
                 plot.axes={points(dat$DPLE.pcp, dat$DPLE.Tmean, pch = 16, cex = 0.7)},
                 color.palette = function(n) hcl.colors(n, "viridis", rev = TRUE)
  )
   
  
   
  # wyrs = subset(dp, wyears <= max(stats.Nyr$wyears))[,1]
  # min = min(wyrs)
  # max = max(wyrs)
  # xlab = "Water Year"
  # 
  # plot(wyrs, subset(dp, wyears <= max)[,2], col = "red", type = "l", ylim = range(dp$DPLE.Tmax, stats$basin_tmax_C), xlab = xlab, ylab = "Tmax (C)")
  # # lines(wyrs, subset(DPLE.ann, wyears >= min)[,2], col = "blue")
  # lines(wyrs, subset(stats, wyears >= min & wyears <= max)[,2])
  # legend("bottomright", legend = c("DP-LE", "Obs"), col = c("red", "black"), lty = 1)
  # 
  # plot(wyrs, subset(dp, wyears <= max)[,3], col = "red", type = "l", ylim = range(dp$DPLE.Tmin, stats$basin_tmin_C), xlab = xlab, ylab = "Tmin (C)")
  # # lines(wyrs, subset(DPLE.ann, wyears >= min)[,3], col = "blue")
  # lines(wyrs, subset(stats.Nyr, wyears >= min & wyears <= max)[,2])
  # legend("bottomright", legend = c("DP-LE", "Obs"), col = c("red", "black"), lty = 1)
  # 
  # plot(wyrs, subset(dp, wyears <= max)[,4], col = "red", type = "l", ylim = range(dp$DPLE.pcp, stats$avg_p_in*25.4), xlab = xlab, ylab = "Precip. (mm)")
  # # lines(wyrs, subset(DPLE.ann, wyears >= min)[,4], col = "blue")
  # lines(wyrs, subset(stats.Nyr, wyears >= min & wyears <= max)[,2]*25.4)
  # legend("bottomright", legend = c("DP-LE", "Obs"), col = c("red", "black"), lty = 1)

  
  ###########################################################################################################################
  ### KNN BLOCK BOOTSTRAPPING
  ###########################################################################################################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### subset annual obs flows to study period
   
  #~~~~~~~~This chunk can be edited to select hindcast period~~~~~~~~~~#
  #set validation start and end years
  #vs = 1980 #forecast starts in year vs and ends in year ve
  #ve = 2013
  
   
  
  ###########################################################################################################################
  ### FORECAST - ANNUAL FLOWS - CESM LE - Years 1-5
  ###########################################################################################################################
  
  
   
  # #~~~~~~~~This chunk can be edited to select experimental groups of different covariates~~~~~~~~~~#
  #
  # # create n data frames for select covariates 
  # # (currently only one in use, previously included CMIP5 and/or seasonal CESM data as separate set(s) of inputs)
  
  mod.1 = DPLE.Nyr
  #remove precip covariate
  mod.1$meanAnoms = mod.1$meanAnoms[,-4]
  mod.1$meanAnomsFlowCor = mod.1$meanAnomsFlowCor[-3,]
  mod.1$NumCovars = 3
  mod.1$means = mod.1$means[,-4]
  
  #
  # #covariates for graph labeling purposes
  c.1 = "CESM-DPLE"
  #covs = cbind(c.1, c.2)
  covs = cbind(c.1)
  
  #specify no. simulations
  #k = ceiling(sqrt(nrow(DPLE.1[[1]])))
  #k = 30
  #nsim = k
  nsim = 300
  
  #n data frames for loop
  a = ncol(covs)
  
  #
  # #~~~~~~~~Do not edit below this chunk~~~~~~~~~~#
  # #~~~~~~~~The following code is designed to run of off the inputs you have selected for hindcast period and covariate selections~~~~~~~~~~#

  ###########################################################################################################################
  ### Hindcast
  ###########################################################################################################################
   
  #create empty dfs to store simulation results
  fcsts.Nyr = NULL
  fcsts.Nyr.bootstrap = NULL
  fcsts.ann = NULL
  fcsts.ann.bootstrap = NULL
  varImp = NULL
  
  i = 1
  for(i in 1:a){
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
    ### Pre-process data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
    #select i-th set of model covariates
    m.list = get(paste0("mod.", i))
    covariates = get(paste0("c.", i))
    
    #extract feature vector and correlations to be used as weights in knn
    mod = data.frame(m.list$meanAnoms)
    wts = matrix(m.list$meanAnomsFlowCor$wts)
    wts[1:length(wts)] = 1 #optionally make weights equal
    mod.RF = data.frame(m.list$means)
    
    #initial validation period length (years)
    if(fcstMode == 'loocv' | fcstMode == 'kFold'){
      vs = min(mod$wyears)
    } else if(fcstMode == 'retroBlind'){
      vs = 2000
    }
    ve = max(mod$wyears)
    vl = ve - vs + 1
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
    ### KNN Block Bootstrap ###
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
    
    j = 20
    for(j in 1:vl){
      
      #select j-th row/year from dataframe for use as predictand
      #subset fitting data to remove future 5-year mean flows or past mean flows/REs that overlap with j-th (current) forecast period 
      if(fcstMode == 'loocv'){
        pred.j = mod[j,]   
        fit_updated = subset(mod, wyears != vs+j-1)
      } else if(fcstMode == 'kFold'){
        pred.j = subset(mod, wyears == vs+j-1)   
        #fit_updated = subset(mod, wyears <= vs+j-1-meanLength | wyears >= vs+j-1+meanLength)
        fit_updated = subset(mod, wyears >= (vs+j-1+max(window.length.re,meanLength)) | wyears <= vs+j-1-meanLength)
      } else if(fcstMode == 'retroBlind'){
        pred.j = subset(mod, wyears == vs+j-1)   
        fit_updated = subset(mod, wyears <= vs+j-1-meanLength)
      }
      
      #1. perform knn resampling for row j (year j) of the validation dataframe----------------------------------
      
      #fitting data with Euclidean distance and weights (wts come from correlation of covariate with future flow)
      #equalWts = matrix(1, nrow = length(wts), ncol = 1)
      #neighbs = knn(fit_updated, pred.j, equalWts)
      neighbs = knn(fit_updated, pred.j, wts)
      #alternatively, fit data with kknn package function
      #knn.pkg = kknn(q.f_maf ~ ., fit_updated[,-1], pred.j, k = nrow(fit_updated))
      #neighbs = data.frame(fit_updated[knn.pkg$C,], dist = knn.pkg$D)
      
      #calc weights for k-neighbors
      k = nrow(fit_updated)
      #knn.wts = (1/(1:k))/sum(1/(1:k))
      #alternatively, weight using Sarah Baker's method
      knn.wts = matrix(NA, k, 1)
      r = 1
      for(r in 1:k){
        knn.wts[r] = (1 - neighbs$dist[r]/tail(neighbs$dist, 1))^2
      }
      #normalize so weights sum to 1
      knn.wts = knn.wts/sum(knn.wts)
      
      #2.  create dataframe of the N years of annual flows corresponding to the time period that generates the selected N year mean neighbor
      #pull starting year of each neighbor
      y = neighbs$wyears

      #create matrix to store annual flow sequences
      #qmtx.ann => matrix of annual flows (N-year sequences based on the N-year means pulled from the knn process)
      qmtx.ann = matrix(NA, nrow = k, ncol = meanLength)
      
      #loop to create N-year sequence for each neighbor (year) selected, pull annual flows for each N year sequence, weight, then store in list for       later summing
      r = 1
      for(r in 1:k){
        q.r = subset(stats, wyears >= y[r] & wyears <= y[r]+(meanLength-1))[,3]
        qmtx.ann[r,] = q.r
      }
      
      #sample nsim times from kxN matrix of N-year annual flow blocks to generate ensemble
      #ens.sims => ensemble simulations
      ens.sims = matrix(NA, nrow = nsim, ncol = meanLength)
      
      for(r in 1:nrow(ens.sims)){
        set.seed(r)
        ind = sample(nrow(qmtx.ann), size = 1, replace = T, prob = knn.wts)
        #ind = sample(nrow(qmtx.ann), size = 1, replace = T)
        ens.sims[r,] = qmtx.ann[ind,]  
      }
      
      #assign simulation years as colnames
      names = (vs+j-1):(vs+j-1+(meanLength-1))
      colnames(qmtx.ann) = names
      colnames(ens.sims) = names
      
      #melt N-year annual flow sequences and further label
      qmtx.ann.melt = melt(qmtx.ann)
      colnames(qmtx.ann.melt) = c("neighborNo", "simYear", "q.maf")
      qmtx.ann.melt$startYear = vs+j-1
      qmtx.ann.melt$fcstYear = rep(1:meanLength, each = k, times = 1)
      qmtx.ann.melt$covariate = paste0(covs[i],", KNN")
      qmtx.ann.melt$obsFlow = rep(subset(stats, wyears >= (vs+j-1) & wyears <= (vs+j-1+(meanLength-1)))[,3], each = k, times = 1)
      
      #melt N-year annual flow sequence bootstraped simulations and further label
      ens.sims.melt = melt(ens.sims)
      colnames(ens.sims.melt) = c("ensembleNo", "simYear", "q.maf")
      ens.sims.melt$startYear = vs+j-1
      ens.sims.melt$fcstYear = rep(1:meanLength, each = nsim, times = 1)
      ens.sims.melt$covariate = paste0(covs[i],", KNN")
      ens.sims.melt$obsFlow = rep(subset(stats, wyears >= (vs+j-1) & wyears <= (vs+j-1+(meanLength-1)))[,3], each = nsim, times = 1)
      
      #reformat the N-year mean flows selected through knn
      means.Nyr = cbind(neighbs[,c(1,ncol(neighbs)-1)], knn.wts)
      means.Nyr$neighborNo = 1:k
      means.Nyr$startYear = vs+j-1
      means.Nyr$simPeriod = paste0((vs+j-1), "-", (vs+j-1+(meanLength-1)))
      means.Nyr$covariate = paste0(covs[i],", KNN")
      means.Nyr$obsFlow = pred.j[1,ncol(pred.j)]
      #calc flow terciles and assign label
      if(means.Nyr$obsFlow[1] >= as.numeric(quantile(stats.Nyr$q.f_maf, 2/3))){
        means.Nyr$flowTercile = "High"
        means.Nyr$obsFlowTercileNo = 3
      } else if(means.Nyr$obsFlow[1] <= as.numeric(quantile(stats.Nyr$q.f_maf, 1/3))){
        means.Nyr$flowTercile = "Low"
        means.Nyr$obsFlowTercileNo = 1
      }  else{
        means.Nyr$flowTercile = "Average"
        means.Nyr$obsFlowTercileNo = 2
      }
      #calc flow terciles for KNN resampled flows
      means.Nyr = flowTercileCalc(means.Nyr, stats.Nyr)

      #calc temp terciles and assign label
      means.Nyr$temp.C = subset(mean.temp.Nyr, wyears == means.Nyr$startYear[1])[1,2]
      if(means.Nyr$temp.C[1] >= as.numeric(quantile(mean.temp.Nyr$t.C, 2/3))){
        means.Nyr$tempTercile = "High"
      } else if(means.Nyr$temp.C[1] <= as.numeric(quantile(mean.temp.Nyr$t.C, 1/3))){
        means.Nyr$tempTercile = "Low"
      }  else{
        means.Nyr$tempTercile = "Average"
      }
      
      #bootstrap to generate an ensemble of N-year mean flows using the neighbors and knn wts
      means.Nyr.sim = data.frame(q.f_maf = sample(means.Nyr$q.f_maf, nsim, replace = T, prob = means.Nyr$knn.wts))
      #means.Nyr.sim = data.frame(q.f_maf = sample(means.Nyr$q.f_maf, nsim, replace = T))
      means.Nyr.sim$ensembleNo = 1:nsim
      means.Nyr.sim$startYear = vs+j-1
      means.Nyr.sim$simPeriod = paste0((vs+j-1), "-", (vs+j-1+(meanLength-1)))
      means.Nyr.sim$covariate = paste0(covs[i],", KNN")
      means.Nyr.sim$obsFlow = pred.j[1,ncol(pred.j)]
      means.Nyr.sim$flowTercile = means.Nyr$flowTercile[1]
      means.Nyr.sim$obsFlowTercile = means.Nyr$obsFlowTercile[1]
      means.Nyr.sim$tempTercile = means.Nyr$tempTercile[1]
      
      #calc flow terciles for bootstrapped flows
      means.Nyr.sim = flowTercileCalc(means.Nyr.sim, stats.Nyr)
      
      #-----------------------------------------------------------------------------------------------------------
      #1. repeat but with random forest instead----------------------------------
      
      #select j-th row/year from dataframe for use as predictand
      #subset fitting data to remove future N-year mean flows or past mean flows/REs that overlap with j-th (current) forecast period 
      if(fcstMode == 'loocv'){
        pred.j = mod.RF[j,] #unscaled covariates  
        fit_updated = subset(mod.RF, wyears != vs+j-1)
      } else if(fcstMode == 'kFold'){
        pred.j = subset(mod.RF, wyears == vs+j-1)   
        #fit_updated = subset(mod.RF, wyears <= vs+j-1-meanLength | wyears >= vs+j-1+meanLength)
        fit_updated = subset(mod.RF, wyears >= (vs+j-1+max(window.length.re,meanLength)) | wyears <= vs+j-1-meanLength)
      } else if(fcstMode == 'retroBlind'){
        pred.j = subset(mod.RF, wyears == vs+j-1)   
        fit_updated = subset(mod.RF, wyears <= vs+j-1-meanLength)
      }
      
      #random forest fitting and prediction
      n.tree = nsim
      
      if(meanLength != 10){
        rf.train = randomForest(q.f_maf ~ ., data = fit_updated[,-1], ntree = n.tree, importance = T)
        #plot variable importance plot for RF
        #varImpPlot(rf.train, main = paste0("RF Variable importance plot for ", pred.j$wyears, " forecast"))
        #varImp.j = varImpPlot(rf.train, main = paste0("RF Variable importance plot for ", pred.j$wyears, " forecast"))
        plot(rf.train)
        #create variable importance for RF
        vip.j = varImpPlot(rf.train, main = paste0("RF - Mean Length = ", meanLength, " years - ", pred.j$wyears))
        vip.jm = melt(vip.j)
        colnames(vip.jm)[1:2] = c("Covariate", "Metric")
        vip.jm$fcstYear = pred.j$wyears
        vip.jm$model = covs[i]
        #store results with other hindcast years
        varImp = rbind(varImp, vip.jm)
      } else{
        rf.train = randomForest(q.f_maf ~ ., data = fit_updated[,-1], ntree = n.tree, importance = F)
      }
      
      #predict ensemble using RF
      fcst.rf = data.frame(q.f_maf = as.numeric(predict(rf.train, newdata = pred.j[-c(1,ncol(pred.j))], predict.all=TRUE)[[2]]))
      rm(rf.train) #remove training model to avoid overloading memory
      
      #2.  create dataframe of the meanLength years of annual flows corresponding to the time period that generates the selected meanLength year mean neighbor
      
      #random forest generates estimates that both match and do not match the training flow data, so for the latter find closet match and use that year for annual flows
      y = matrix(NA, nrow = n.tree, ncol = 1)
      for(r in 1:n.tree){
        y[r] = fit_updated[which.min(abs(fcst.rf$q.f_maf[r] - fit_updated$q.f_maf)), 1]
      }
      
      #create matrix to store annual flow sequences
      #ens.sims => matrix of annual flows (meanLength-year sequences based on the meanLength-year means pulled from the rf process)
      ens.sims = matrix(NA, nrow = n.tree, ncol = meanLength)
      
      #loop to create meanLength-year sequence for each neighbor (year) selected, pull annual flows for each meanLength year sequence, weight, then store in list for       later summing
      r = 1
      for(r in 1:n.tree){
        q.r = subset(stats, wyears >= y[r] & wyears <= y[r]+(meanLength-1))[,3]
        ens.sims[r,] = q.r
      }
      
      #assign simulation years as colnames
      names = (vs+j-1):(vs+j-1+(meanLength-1))
      colnames(ens.sims) = names
      
      #melt meanLength-year annual flow sequence simulations and further label
      ens.sims.melt.RF = melt(ens.sims)
      colnames(ens.sims.melt.RF) = c("ensembleNo", "simYear", "q.maf")
      ens.sims.melt.RF$startYear = vs+j-1
      ens.sims.melt.RF$fcstYear = rep(1:meanLength, each = n.tree, times = 1)
      ens.sims.melt.RF$covariate = paste0(covs[i],", RF")
      ens.sims.melt.RF$obsFlow = rep(subset(stats, wyears >= (vs+j-1) & wyears <= (vs+j-1+(meanLength-1)))[,3], each = n.tree, times = 1)
      
      
      #reformat RF ensemble of meanLength-year mean flows ----------------------------------
      means.Nyr.sim.RF = data.frame(q.f_maf = fcst.rf)
      #means.Nyr.sim = data.frame(q.f_maf = sample(means.Nyr$q.f_maf, nsim, replace = T))
      means.Nyr.sim.RF$ensembleNo = 1:n.tree
      means.Nyr.sim.RF$startYear = vs+j-1
      means.Nyr.sim.RF$simPeriod = paste0((vs+j-1), "-", (vs+j-1+(meanLength-1)))
      means.Nyr.sim.RF$covariate = paste0(covs[i],", RF")
      means.Nyr.sim.RF$obsFlow = pred.j[1,ncol(pred.j)]
      means.Nyr.sim.RF$flowTercile = means.Nyr$flowTercile[1]
      means.Nyr.sim.RF$obsFlowTercile = means.Nyr$obsFlowTercile[1]
      means.Nyr.sim.RF$tempTercile = means.Nyr$tempTercile[1]
      
      #calc flow terciles for bootstrapped flows
      means.Nyr.sim.RF = flowTercileCalc(means.Nyr.sim.RF, stats.Nyr)

      #assign to storage for the current covariate's full validation period--------------------------------
      fcsts.ann = rbind(fcsts.ann, qmtx.ann.melt)
      fcsts.ann.bootstrap = rbind(fcsts.ann.bootstrap, ens.sims.melt, ens.sims.melt.RF)
      fcsts.Nyr = rbind(fcsts.Nyr, means.Nyr)
      fcsts.Nyr.bootstrap = rbind(fcsts.Nyr.bootstrap, means.Nyr.sim, means.Nyr.sim.RF)
      
      
      ##########################################################################
      
    }
    
    #put all results into list
    sim.results = list(fcsts.ann, fcsts.ann.bootstrap, fcsts.Nyr, fcsts.Nyr.bootstrap )
    names(sim.results) = c("AnnualFlow", "AnnualFlowBootstrap", "NyrMeanFlow", "NyrMeanFlowBootstrap")
    
    #################################################################################################
  }
  
  if(meanLength != 10){
      vip = subset(varImp, Metric == "%IncMSE") %>%
      mutate(leadTime = as.factor(model),
             name = reorder_within(Covariate, value, model)) %>%
      ggplot(aes(x = name, y = value)) +
      geom_boxplot() +
      facet_wrap(~model, scales = "free", shrink = F) +
      coord_flip() +
      theme_bw() +
      theme(axis.title = element_text(face = "bold"), text=element_text(size=16)) +
      #geom_hline(yintercept = 0, col = "red") +
      scale_x_reordered() +
      #scale_y_continuous(expand = c(0,0)) +
      labs(y = "% Increase in MSE if variable is randomly permuted",
           x = "Covariate", 
           title = "Hindcast aggregated variable importance",
           subtitle = paste0(meanLength, "-year means for flow and temperature"))
    
    print(vip)
    
    plotPath = "K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Code/plots/midtermForecasts_multiYearMeans/"
    title = paste0("varImportance_CESM.", meanLength, "yr_fcstMode=", fcstMode, ".png")
    ggsave(title, plot = print(vip), device = "png", path = plotPath, width = 7, height = 5, dpi = 450, units = "in")
    
    # vip2 = subset(varImp, Metric == "%IncMSE") %>%
    #   mutate(leadTime = as.factor(model),
    #          name = reorder_within(Covariate, value, model)) %>%
    #   ggplot(aes(x = name, y = value)) +
    #   geom_density_ridges(scale = 5, alpha = 0.7) +
    #   labs(y = "% Increase in MSE if variable is randomly permuted",
    #        x = "Covariate")
    # 
    #   print(vip2)
    #   
    #   plotPath = "K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Code/plots/midtermForecasts_multiYearMeans/"
    #   title = paste0("varImportance.PDF_CESM.", meanLength, "yr_fcstMode=", fcstMode, ".png")
    #   ggsave(title, plot = print(vip2), device = "png", path = plotPath, width = 6, height = 6, dpi = 450, units = "in")
      
  }
  
  #examine esp
  
  #esp.test = esp$locAll_fcst

  ### Plot 5-year mean flow (no bootstrapping)
  # df = sim.results$`NyrMeanFlow`
  # 
  # ggplot(df, aes(x = startYear, y = q.f_maf, group = startYear, color = tempTercile)) + 
  #   geom_boxplot() + 
  #   geom_point(aes(x = startYear, y = obsFlow), color = "black") +
  #   xlab("Year") + ylab("5-year mean flow (MAF)")
  
  ### Plot 5-year mean flow (bootstrapping)
  df = sim.results$`NyrMeanFlowBootstrap`
  n = length(unique(df$startYear))
  
  i = 1
  for(i in 1:n){
    dat.i = subset(stats.Nyr, wyears < min(df$startYear)+i-1)
  }
  
  
  p = ggplot(df, aes(x = startYear, y = q.f_maf, group = interaction(startYear, covariate), fill = covariate)) + 
    geom_boxplot() + labs(colour = "Covariates & Model") +  theme_bw() +
    #geom_point(aes(x = startYear, y = obsFlow), color = "black") +
    geom_crossbar(data =  df, aes(x = startYear, y = obsFlow, ymin=obsFlow, ymax=obsFlow), position=position_dodge(), color="blue", fatten = 1.5) + 
    xlab("Year") + ylab(paste0(meanLength, "-year mean flow (MAF)")) +
    theme(axis.title = element_text(face = "bold")) +
    geom_hline(yintercept = clim.Nyr, col = "black", linetype="dashed") +
    scale_x_continuous(breaks=seq(min(df$startYear), max(df$startYear), 2)) +  
    scale_linetype_manual(name = "Climatology", values = 1, guide = guide_legend(override.aes = list(color = c("purple"))))
  
  p2 = p + geom_boxplot()
  
  ### Plot 5-year mean flow (bootstrapping) - CESM-DPLE only 
  
  df = sim.results$`NyrMeanFlowBootstrap`
  df = df[str_detect(df$covariate, "CESM-DPLE"),]
  
  plot = ggplot(df, aes(x = startYear, y = q.f_maf, group = interaction(startYear, covariate), fill = covariate)) + 
    geom_boxplot() + labs(colour = "Covariates & Model") + theme_bw() +
    #geom_point(aes(x = startYear, y = obsFlow), color = "black", size = 1.5) +
    geom_crossbar(data =  df, aes(x = startYear, y = obsFlow, ymin=obsFlow, ymax=obsFlow), position=position_dodge(), color="blue", fatten = 1.5) + 
    xlab("Year") + ylab(paste0(meanLength, "-year mean flow (MAF)")) +
    theme(axis.title = element_text(face = "bold")) +
    geom_hline(yintercept = clim.Nyr, col = "black", linetype="dashed") +
    scale_x_continuous(breaks=seq(min(df$startYear), max(df$startYear), 2)) +  
    scale_linetype_manual(name = "Climatology", values = 1, guide = guide_legend(override.aes = list(color = c("purple"))))
  
  plot
  print(plot)
  
  plotPath = "K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Code/plots/midtermForecasts_multiYearMeans/"
  title = paste0("hindcast_CESM.", meanLength, "yr_fcstMode=", fcstMode, ".png")
  ggsave(title, plot = print(plot), device = "png", path = plotPath, width = 10, height = 5, dpi = 450, units = "in")
  
  
  ### Plot N-year mean flow (bootstrapping) - CESM-DPLE only and KNN only

  df = sim.results$`NyrMeanFlowBootstrap`
  df = df[str_detect(df$covariate, "CESM-DPLE, KNN"),]
  
  plot = ggplot(df, aes(x = startYear, y = q.f_maf, group = interaction(startYear, covariate), colour = covariate)) + 
    geom_boxplot() + labs(colour = "Covariates & Model") + theme_light() +
    geom_point(aes(x = startYear, y = obsFlow), color = "black", size = 1.5) +
    xlab("Year") + ylab(paste0(meanLength, "-year mean flow (MAF)")) +
    theme(axis.title = element_text(face = "bold")) +
    geom_hline(yintercept = clim.Nyr, col = "black", linetype="dashed") +
    scale_x_continuous(breaks=seq(min(df$startYear), max(df$startYear), 2)) +  
    scale_linetype_manual(name = "Climatology", values = 1, guide = guide_legend(override.aes = list(color = c("purple"))))
  
  plot
  
  title = paste0("hindcast_CESM_KNN.", meanLength, "yr_fcstMode=", fcstMode, ".png")
  ggsave(title, plot = print(plot), device = "png", path = plotPath, width = 15, height = 5, dpi = 450, units = "in")
  
  
  ### Plot annual flow (no bootstrapping)
   
  # i = 1
  # for(i in 1:vl){
  # 
  # df = subset(sim.results$AnnualFlow, startYear == vs+i-1)
  # 
  # p = ggplot(df, aes(x = simYear, y = q.maf, group = fcstYear)) + 
  #   geom_boxplot() + 
  #   geom_point(aes(x = simYear, y = obsFlow), color = "red") +
  #   xlab("Year") + ylab("Flow (MAF)") +
  #  theme(axis.title = element_text(face = "bold"))
  # 
  # print(p)
  # 
  # }
  
  ### Plot annual flow (bootstrapping)
   
  # i = 1
  # for(i in 1:vl){
  # 
  # df = subset(sim.results$AnnualFlowBootstrap, startYear == vs+i-1)
  # 
  # p = ggplot(df, aes(x = simYear, y = q.maf, group = fcstYear)) +
  #   geom_boxplot() +
  #   geom_point(aes(x = simYear, y = obsFlow), color = "red") +
  #   xlab("Year") + ylab("Flow (MAF)") +
  #   theme(axis.title = element_text(face = "bold"))
  # 
  # print(p)
  # 
  # }
  
  
  ### output flow projections to RDS file
  setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/")
  fname = paste("fcst.sims_k=", k, "_nsims=", nsim, "_yrs=", vs, ve, "_meanLength=", meanLength, "yr_fcstMode=", fcstMode,".Rds", sep = "")
  saveRDS(sim.results, file = fname)
  
  ##############################################
  ### Skill score calcs
  ##########
  
  ### Plotting
  ss.Nyr.bootstrap = ss.calc.Nyr(sim.results, T, T, T, T, meanLength, vl, vs, ve, stats.Nyr, nsim, fcstMode)
  
  ss.Nyr.bootstrap.flowTerciles = ss.Nyr.bootstrap[[1]]
  
  boxplot(subset(ss.Nyr.bootstrap.flowTerciles, startYear >= 1980 & variable == "RPSS")[,5])
  abline(h = 0, col = "red")
  
  ### output flow skill calcs to RDS file
  setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/")
  fname = paste("fcst.skill_k=", k, "_nsims=", nsim, "_yrs=", vs, ve, "_meanLength=", meanLength, "yr_fcstMode=", fcstMode,".Rds", sep = "")
  saveRDS(ss.Nyr.bootstrap, file = fname)
  
   
  ###Ensemble mean RMSE (5 year mean projected flows)
  meanFlow.Nyr.sim = sim.results$`NyrMeanFlowBootstrap`
  
  #transform data into correct format for EnsRps function
  years = unique(meanFlow.Nyr.sim$startYear)
  nYrs = length(years)
  mods = unique(meanFlow.Nyr.sim$covariate)
  nMods = length(mods)
  
  ensMean.Nyr.fcst.list = vector("list", nMods)
  
  j = 1
  for(j in 1:nMods)
  {
    
    ensMean.Nyr.fcst = matrix(NA, nrow = nYrs, ncol = 1)
    mat.Nyr.obs = matrix(NA, nrow = nYrs, ncol = 1)
    
    i = 1
    for(i in 1:nYrs){
      
      df.i = subset(meanFlow.Nyr.sim, startYear == years[i] & covariate == mods[j])
      
      ensMean.Nyr.fcst[i] = mean(df.i$q.f_maf)
      
      mat.Nyr.obs[i] = df.i$obsFlow[1]
      
    }
    
    ensMean.Nyr.fcst.list[[j]] = ensMean.Nyr.fcst
    
  }
  
  mat.clim.Nyr = matrix(rep(clim.Nyr, nYrs))
  rmse.Nyr.clim = rmse(mat.clim.Nyr, mat.Nyr.obs)
  
  rmse.Nyr.list = vector("list", nMods)
  
  i = 1
  for(i in 1:nMods){
    rmse.Nyr.list[[i]] = rmse(ensMean.Nyr.fcst.list[[i]], mat.Nyr.obs)
  }
  
  rmse.Nyr.list[[3]] = rmse.Nyr.clim
  names(rmse.Nyr.list) = c(mods, "Climatology")
  
  ############################
  ### crpss on annual flows
  ### Plotting
   
  
  # crpss.ann.bootstrap = crpss.calc.ann(sim.results, T, T, vl, vs, stats, meanLength)
  # #crpss.ann.noBootstrap = crpss.calc.ann(sim.results, F)
  # 
  # ### output annual flow skill calcs to RDS file
  # setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/")
  # fname = paste("crpss.ann.fcst.skill_k=", k, "_nsims=", nsim, "_yrs=", vs, ve, "_meanLength=", meanLength,".Rds", sep = "")
  # saveRDS(crpss.ann.bootstrap, file = fname)
  
  ############################
  ### rpss on annual flows
  
  for(i in 1:nrow(stats)){
    
    if(stats$q_maf[i] >= as.numeric(quantile(stats$q_maf, 2/3))){
      stats$obsFlowTercile[i] = "High"
      stats$obsFlowTercileNo[i] = 3
    } else if(stats$q_maf[i] <= as.numeric(quantile(stats$q_maf, 1/3))){
      stats$obsFlowTercile[i] = "Low"
      stats$obsFlowTercileNo[i] = 1
    }  else{
      stats$obsFlowTercile[i] = "Average"
      stats$obsFlowTercileNo[i] = 2
    }
    
  }
  
  ### Plotting
  # rpss.ann.bootstrap = rpss.calc.ann(sim.results, T, T, vl, vs, stats)
  # #crpss.ann.noBootstrap = crpss.calc.ann(sim.results, F)
  # 
  # ### output annual flow skill calcs to RDS file
  # setwd("K:/My Drive/Phd Research/CRB Midterm Temperature Perturbed Predictions/Data/Phase 1/knn_fcsts/")
  # fname = paste("rpss.ann.fcst.skill_k=", k, "_nsims=", nsim, "_yrs=", vs, ve, "_meanLength=", meanLength,".Rds", sep = "")
  # saveRDS(rpss.ann.bootstrap, file = fname)
  
  return(list(sim.results, ss.Nyr.bootstrap.flowTerciles, rmse.Nyr.list, data.frame(clim.Nyr), ss.Nyr.bootstrap[[2]], ss.Nyr.bootstrap[[3]], ss.Nyr.bootstrap[[4]], ss.Nyr.bootstrap[[5]]))
  
}


