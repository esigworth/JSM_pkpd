###Simulation study code
source("R/three-compartment-model.R")
library(mvtnorm)


###Pull same prior mean and variance vectors from original model
vec_prior_mean <- c("lk10"=log(0.770)+log(19.01),
                    "lk12"=log(3.156)+log(19.01),
                    "lk21"=log(3.156)+log(19.01),
                    "lk13"=1+log(19.01),
                    "a1"=-13.65,
                    "la2"=log(4),
                    "b"=-0.0959*1e5) #needs to be scaled up if concentration values are scaled down
vec_prior_var  <- diag(rep(c(10,10,10,10,10,4000),c(1,2,1,1,1,1))) 


##Create empty list for filling with formatted original subject profiles and parameters
subject_profile <- list()
par_profile <- list()

##Loop through original subject list and...
for(i in 1:length(par_opt_test_list_last_dose)){
  ##Extract case-by-case
  case <- par_opt_test_list_last_dose[[i]]$case
  ##Put back together the training and testing sets
  tof <- rbind(case$tof, case$tof_test)
  ##Replace tof with all tof counts and drop testing portion
  case$tof <- tof
  case$tof_test <- NULL
  ##Add all case information to subject profile
  subject_profile[[i]] <- case
  ##Fit OLR model to all TOF counts and store parameters
  par_profile[[i]] <- olr_fit(case)$par
}

##Function to predict probabilities of each TOF count at each time
##i = element of subject profile (1:140), incr=frequency of training TOF counts
pred_probs <- function(i,incr){
  ##Pull the parameters and case from subject i
  pars <- par_profile[[i]]
  case <- subject_profile[[i]]
  ##Split the case into training and testing base on timing of last dose
  lastdose <- unclass(case$meds$RelDoseTime[nrow(case$meds)])[1]
  tof_trn <- case$tof[which(unclass(case$tof$RelTOFTime) < lastdose),]
  tof_tst <- case$tof[which(unclass(case$tof$RelTOFTime) >= lastdose),]
  ##Create training times as a sequence from 0 to the last increment of incr before last dose
  tof_trn <- seq(0, ((lastdose)-incr),
                 by=incr)
  ##Convert times to hours
  times <- c(tof_trn, tof_tst$RelTOFTime)/60
  ##Calculate probabilities at training and testing times 
  probs <- olr_predict(times, pars, case)
  ##Return probabilities, time sequence, and parameters
  list(probs,times,pars)
}


##Function to pick simulated TOF counts based on predicted probabilities
##i = subject, incr = time increment
pick_counts <- function(i,incr){
  ##Predict probabilities at time increments
  out <- pred_probs(i,incr)
  ##Pull probabilities, times, and parameters from output
  mat <- out[[1]]
  times <- out[[2]]
  pars <- out[[3]]
  ##Create empty vector to fill with simulated counts
  choices <- rep(NA, nrow(mat))
  ##Cycle through times and...
  for(i in 1:nrow(mat)){
    ##Use multinom to pick a count from 1 to 5 based on probabilities
    v <- rmultinom(n=1, size=1, prob=mat[i,])
    ##Shift selection to be on 0-4 scale
    choices[i] <- which(v==1) - 1
  }
  ##Return simulated counts, times, and parameters
  list(choices, times, pars)
}


##Function to generate n subjects with freq timing of TOF counts
simulate_people <- function(n, freq){
  ##Create empty list
  subj <- list()
  ##Randomly sample n profiles from our total of 140
  profiles <- sample(1:140, n, replace=TRUE)
  ##Initiate subject ID at 1
  b <- 1
  for(i in profiles){
    ##Pull profile i from subject profiles
    pers <- subject_profile[[i]]
    ##Generate counts at given frequency 
    v <- pick_counts(i,freq)
    ##Pull counts, times, and "true" parameters from output
    vals <- v[[1]]
    times <- v[[2]]
    fit_full <- v[[3]]
    
    ##Create formatted tof object 
    TOFValue <- vals
    RelTOFTime <- times
    PatientID <- rep(b, length(vals))
    TOFParameterID <- rep(pers$tof$TOFParameterID, length(vals))
    
    ##Format TOF object and add to case
    pers$tof <- data.frame(PatientID, TOFValue, RelTOFTime=(RelTOFTime*60))
    
    ##Update patient ID number
    pers$demo$PatientID <- b
    pers$meds$PatientID <- b

    ##Separate TOF values into training and testing
    pers$tof$TOFValue[which(pers$tof$RelTOFTime <= 0)] <- 4
    pers$tof$TOFValue   <- factor(pers$tof$TOFValue, levels=c(0:4), ordered=TRUE)
    last_dose <- unclass(pers$meds$RelDoseTime[nrow(pers$meds)])[1]
    tof_trn <- pers$tof[which(unclass(pers$tof$RelTOFTime) < last_dose),]
    tof_tst <- pers$tof[which(unclass(pers$tof$RelTOFTime) >= last_dose),]
    pers$tof <- tof_trn
    pers$tof_test <- tof_tst
    
    ##Fit OLR to training values to get "estimated" parameters
    fit <- try(olr_fit(pers))
    
    ##Format case and store in subjects
    subj[[b]] <- list(id=b, case=pers, fit=fit, fit_full=fit_full)
    
    ##Index patient ID
    b <- b+1
  }
  ##Return list of subjects
  return(subj)
}


###Assess predictive accuracy

assess <- function(patients){
  probs <- c()
  predicted_list <- lapply(patients, function(x) {
    ##Compute predicted values for each TOF test time
    if(inherits(x$fit, "try-error"))
      return(NA)
    make_pred <- olr_classify(x$case$tof_test$RelTOFTime/60,
                              x$fit$par, x$case)
    cpred <- as.numeric(make_pred[[1]]) - 1
    probs <- make_pred[[2]]
    
    val <- matrix(c(cpred,probs),nrow=length(cpred),byrow=FALSE)
    
    return(val)
  })
  
  ##Unlist returned predictions
  names(predicted_list) <- paste("Subj_", 1:length(predicted_list), sep="")
  predicted <- do.call("rbind", predicted_list)[,1]
  
  
  ##Pull list of actual TOF values from patient profile
  
  actual_list <- lapply(patients, function(x) {
    ## compute predicted probabilities for each TOF count
    if(inherits(x$fit, "try-error"))
      return(NA)
    cactl <- as.numeric(as.character(x$case$tof_test$TOFValue))
    return(cactl)
  })
  
  ##Unlist returned true values
  names(actual_list) <- names(predicted_list)
  actual <- unlist(actual_list)
  
  ##Calculate proportion of correct predictions
  (val_last_dose <- mean(actual == predicted, na.rm=TRUE))
  
  
  
  ##Create matrix of types of disparity in predictions
  combos <- rep(NA,length(predicted))
  for(i in 1:length(predicted)) {
    if(is.na(actual[i])){
      combos[i] <- NA
    }
    else if(actual[i]==0) {
      combos[i] <- letters[predicted[i]+1]
    }
    else if(actual[i]==1) {
      combos[i] <- letters[predicted[i]+6]
    }
    else if(actual[i]==2) {
      combos[i] <- letters[predicted[i]+11]
    }
    else if(actual[i]==3) {
      combos[i] <- letters[predicted[i]+16]
    }
    else if(actual[i]==4) {
      combos[i] <- letters[predicted[i]+21]
    }
  }
  
  combos <- combos[!is.na(combos)]
  to.add <- which(!(letters[1:25] %in% combos))
  combos_table <- c(table(combos), rep(0,length(to.add)))
  if(length(to.add)!=0){
    names(combos_table)[(25-length(to.add)+1):25] <- letters[to.add]
  }
  else {names(combos_table) <- letters[1:25]}
  combos_table <- combos_table[sort(names(combos_table))]
  
  
  combo_mat <- matrix(combos_table,nrow=5,byrow = TRUE)
  colnames(combo_mat) <- c("0","1","2","3","4")
  rownames(combo_mat) <- c("0","1","2","3","4")
  
  ## table of percentages to show relative accuracy at each combo
  combo_props <- round(prop.table(combo_mat,margin=2)*100,2)
  
  (prop_diff <- 100*round(prop.table(table((predicted-actual))),4))
  
  ##Find proportion of predictions that are within 1 count of the true value
  (within_1 <- sum(prop_diff[c("-1","0", "1")]))
  
  ##Return list of accuracy measures
  return(list(N_obs=length(actual),Accuracy=val_last_dose, Combos_numeric=combo_mat, 
              Combos_prop=combo_props, Diffs_prop = prop_diff, 
              Within_one=within_1))
}

###Calculate mean absolute percent error between full and test parameter fits
mape <- function(subj){
  ##Parameters from only training data
  red <- subj$fit$par
  ##Parameters from all data
  full <- subj$fit_full
  ##Calculate MAPE
  perc_error <- abs((full-red)/full)*100
}

###Create data frame of all MAPE values, pull mean and sd
comb <- function(patients){
  ##apply MAPE function to all simulated patients
  errors <- try(lapply(patients, mape))
  ##Put all errors into a dataframe
  frame <- do.call("rbind",errors)
  ##Find mean and SD
  means <- apply(frame, 2, mean)
  sds <- apply(frame, 2, sd)
  ##Return list of MAPE, SE, and actual value of errors
  return(list("Summary"=cbind("MAPE"=means, "SE"=sds),"Values"=frame))
}

###Do all steps in one function
gen_assess <- function(n, freq, seed=sample(x = 1:50000,1)){
  ##Potentially set seed to get consistent results
  set.seed(seed)
  #Simulate n patients with freq timing
  patients <- simulate_people(n,freq)
  ##Assess accuracy of all predictions
  valid <- assess(patients)
  ##Calculate MAPE for all predictions
  #errors <- comb(patients)
  return(list(Validate=valid, Subjects=patients))
}



###Generate populations of 200 and find accuracy and errors at 1 min, 10 minutes, 20 minutes

###Every minute

every_min <- gen_assess(200,1)

###Every 10 minutes

every_ten <- gen_assess(200,10)


###Every 20 minutes

every_twenty <- gen_assess(200,20)

# ##Pull accuracy values into table
# comp_acc <- rbind(every_sec$Validate$Accuracy,
#                   every_min$Validate$Accuracy, 
#                   every_ten$Validate$Accuracy, 
#                   every_twenty$Validate$Accuracy)
# 
# rownames(comp_acc) <- c("30 seconds", "1 minute", "10 minutes", "20 minutes")
# colnames(comp_acc) <- "Accuracy"
# 
# ##Pull distribution of predictive errors into table
# comp_diffs <- rbind(every_sec$Validate$Diffs_prop,
#                     every_min$Validate$Diffs_prop,
#                     every_ten$Validate$Diffs_prop,
#                     every_twenty$Validate$Diffs_prop)
# 
# rownames(comp_diffs) <- c("30 seconds", "1 minute", "10 minutes", "20 minutes")
# 
# ##Pull how often each timing results in "within one" predictions
# comp_within <- rbind(every_sec$Validate$Within_one,
#                       every_min$Validate$Within_one, 
#                       every_ten$Validate$Within_one, 
#                       every_twenty$Validate$Within_one)
# 
# rownames(comp_within) <- c("30 seconds", "1 minute", "10 minutes", "20 minutes")
# colnames(comp_within) <- "Within One"
# 
# 
# ##Pull summary of errors into table
# mape <- cbind(every_sec$MSE$Summary[,1], 
#               every_min$MSE$Summary[,1],
#               every_ten$MSE$Summary[,1],
#               every_twenty$MSE$Summary[,1])
# 
# rownames(mape) <- names(every_sec$MSE$Summary[,1])
# colnames(mape) <- c("30 seconds", "1 minute", "10 minutes", "20 minutes")
# 
# tab
# round(mape,2)














