plot_case_test <- function(case, par, pred=NULL, ext=90, cex=1.4, line=3,
                      showall=FALSE, pch=4) {
  
  demo_sub <- case$demo
  meds_sub <- case$meds
  tof_sub <-  case$tof
  
  ## combine doses that occur at the same time
  time_drug <- interaction(meds_sub$DoseTime, meds_sub$DrugName)
  if(length(time_drug) > length(unique(time_drug))) {
    meds_sub <- 
      do.call(rbind, lapply(unique(time_drug), function(td) {
        meds_sub_sub <- subset(meds_sub, time_drug == td)
        ## ensure that DoseUnit is consistent
        if(length(unique(meds_sub_sub$DoseUnit)) > 1)
          stop("different DoseUnit's for same drug")
        total_dose <- sum(meds_sub_sub$DoseAmount)
        meds_sub_sub <- meds_sub_sub[1,]
        meds_sub_sub$DoseAmount <- total_dose
        meds_sub_sub
      }))
  }
  
  ## if there is a test object and predictions are not null, determine which are accurate
  if(!is.null(case$tof_test) & !is.null(pred)) {
    ## determine which predictions are the same as actual records
    observed <- unclass(case$tof_test$TOFValue)-1
    observed.times <- case$tof_test$RelTOFTime
    predicted <- pred
    plotting <- data.frame(observed.times,observed,predicted)
  }
  
  ## get time of first dose
  first_dose_time <- min(meds_sub$DoseTime)
  ## compute dose times relative to first dose
  rel_dose_time <- difftime(meds_sub$DoseTime,
                            first_dose_time, units="mins")
  ## compute TOF times relative to first dose
  rel_tof_time <- difftime(tof_sub$TOFTime,
                           first_dose_time, units="mins")
  
  ## standardize and extend time range
  xrng <- range(c(rel_dose_time, rel_tof_time))
  xrng[2] <- xrng[2] + ext
  ## plot full TOF count range 
  yrng <- c(0,4) 
  
  ## change layout when full model details are plotted
  if(showall) {
    layout(matrix(c(1,1,2,3,4,5),6,1))
  } else {
    layout(matrix(c(1,2),2,1))
  }
  
  ## set margins
  par(omi=c(1,1,1,2), mai=rep(0,4))
  
  ## create TOF count plot
  plot(rel_tof_time, unclass(tof_sub$TOFValue)-1,
       xaxt="n", yaxt="n",
       xlab="", ylab="",
       pch=pch, lwd=2,
       ylim=yrng,
       xlim=xrng,
       cex=1.5)
  
  ## if there is a tof_test object and pred, plot those records and predictions with different symbols
  if(!is.null(case$tof_test) & !is.null(pred)) {
    points(plotting$observed.times,
           plotting$observed,
           pch=pch, cex=1.5, lwd=2)
    
    points(plotting$observed.times,
           plotting$predicted,
           pch=0, cex=1.5, lwd=2)
  }
  
  ## add dashed horizontal lines at each TOF count value
  for(tof_val in 0:4)
    abline(h=tof_val, lty=3, col="lightgray")
  
  ## add y-axis
  axis(2, cex.axis=cex); mtext("TOF Count", 2, line=line, cex=cex)
  
  ## add main title
  mtext(paste0("Patient ID: ", demo_sub$PatientID), 3,
        line=line+0.1, cex=cex)
  
  
  
  ## add dashed vertical lines corresponding to doses
  ## add drug and amount labels to top of plot
  for(i in 1:nrow(meds_sub)) {
    abline(v=rel_dose_time[i],
           lty=2, col=meds_sub$DrugID[i])
    text(x=rel_dose_time[i], y=par("usr")[4], #yrng[1],
         offset=0.2, xpd=NA, adj=c(0,0), #pos=4,
         labels=paste0(#meds_sub$DrugAbbr[i], ". ",
           round(meds_sub$DoseAmount[i]*1000,1),
           "mg"), cex=0.8*cex, srt=40)
  }
  
  ## create ivt_list, for doses (needed for pk model)
  rel_max_time <- max(c(rel_tof_time, rel_dose_time)) + ext
  ivt_list <- list()
  for(i in 1:nrow(meds_sub)) {
    if(meds_sub$DrugAbbr[i] %in% 
       c("Vec", "Roc", "Atr", "Cis")) {
      ivt_list[[length(ivt_list) + 1]] <- 
        list(time=as.numeric(rel_dose_time[i]/60),
             dose=meds_sub$DoseAmount[i]/(1/60))
    }
  }
  
  ## compute PK model time solutions for each compartment
  pk_sol <- 
    pk_solution(k10=exp(par["lk10"]),k12=exp(par["lk12"]),k21=exp(par["lk21"]),
                k13=exp(par["lk13"]),
                ivt=case$ivt,
                tmax=case$tmax)
  ptime <- seq(0, as.numeric(rel_max_time)/60, 1/60)
  psols <- t(pk_sol(ptime))
  
  ## compute probabilities for each TOF count from OLR model
  pprbs <- olr_predict(ptime, par, case)
  
  ## plot probabilities from OLR model
  for(i in 1:ncol(pprbs))
    points(ptime*60, rep(i-1,length(ptime)),
           pch=16, cex=3*pprbs[,i], col=rgb(0,0,0,0.05))
  
  if(showall) {
    plot(ptime, psols[,1], type="l",
         xaxt="n", yaxt="n",
         xlab="", ylab="", xlim=xrng/60)
    axis(4, cex.axis=cex); mtext("Cent. Conc.",4, line=line, cex=cex)
    for(i in 1:nrow(meds_sub)) {
      abline(v=rel_dose_time[i]/60,
             lty=2, col=meds_sub$DrugID[i])
    }
    
    plot(ptime, psols[,2], type="l",
         xaxt="n", yaxt="n",
         xlab="", ylab="", xlim=xrng/60)
    axis(2, cex.axis=cex); mtext("Peri. Conc.",2, line=line, cex=cex)
    for(i in 1:nrow(meds_sub)) {
      abline(v=rel_dose_time[i]/60,
             lty=2, col=meds_sub$DrugID[i])    
    }
  }
  
  plot(ptime, psols[,3], type="l",
       xaxt="n", yaxt="n",
       xlab="", ylab="", xlim=xrng/60)
  abline(h=0, lty=3, col="gray")
  axis(2, at=0, cex.axis=cex)
  mtext("Biophase Conc.",2, line=line, cex=cex)
  
  for(i in 1:nrow(meds_sub)) {
    abline(v=rel_dose_time[i]/60,
           lty=2, col=meds_sub$DrugID[i])
  } 
  
  if(showall) {
    oprbs <- olr_predict(tof_sub$RelTOFTime/60, par, case)
    plot(ptime, pprbs[,1], type="n",
         xaxt="n", yaxt="n",
         xlab="", ylab="",
         xlim=xrng/60, ylim=c(0,1))
    
    #pcols <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15")
    #pcols <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
    pcols <- c("#1A9641", "#A6D96A", "#FFFF30", "#FDAE61", "#D7191C")
    for(i in 1:ncol(pprbs)) {
      lines(ptime, pprbs[,i], col=pcols[i], lwd=2)
    }
    for(i in 1:nrow(meds_sub)) {
      abline(v=rel_dose_time[i]/60,
             lty=2, col=meds_sub$DrugID[i])
      
    } 
    points(tof_sub$RelTOFTime/60, oprbs[cbind(1:nrow(tof_sub), unclass(tof_sub$TOFValue))], pch=4)
    legend("topright", colnames(pprbs), col=pcols, lwd=2, bty="n", cex=cex)
    axis(2, cex.axis=cex); mtext("Pred. Prob.",2, line=line, cex=cex)
  }
  
  axis(1, cex.axis=cex); mtext("Time relative to first dose (h)",1, line=line, cex=cex)
  par(cex=1)
  
}