


## Three-compartment Pharmacodynamic ("biophase") Model 
## This model implements the three-compartment model similar to that 
## described in Hull et al. (1978) "A pharmacodynamic model for pancuronium" BJA v50.
## The model is a simple extension of the two-compartment PK model, 
## where an additional compartment is added to reflect the "biophase,"
## or region where the drug is active.a

## dm1/dt = -k12*m1 +k21*m2 -k13*m1 +k31*m3 -k10*m1
## dm2/dt =  k12*m1 -k21*m2
## dm3/dt =                 +k13*m1 -k31*m3

## K = matrix(c(-k12-k13-k10,  k21,  k31,
##               k12,         -k21,  0,
##               k13,               -k31),
##            3,3,byrow=TRUE)


## The CM3 function follows 
## "Appendix: Two-compartment Model Solution" from 
## Shotwell et al. (2014) "Optimal Design of
## Perturbations For Individual Two-Compartment
## Pharmacokinetic Analysis"
CM3 <- function(dose, k10, k12, k21=k12, k13, k31=k13*1e4,
                init=c(0,0,0)) {
  K = matrix(
    c(-k12-k13-k10,  k21,  k31,
       k12,         -k21,    0,
       k13,            0, -k31),
    3,3,byrow=TRUE)
  E <- eigen(K)
  r <- qr.solve(E$vectors, init+c(dose,0,0))
  function(t)
    E$vectors %*% diag(exp(E$values * t)) %*% r
}

## The following code attempts to reproduce figure 6 of Hull et al. (1978).
## See "docs/Hull1978.pdf". However, it looks like Hull et al used umol
## rather than nmol, and it also appears that the biophase concentration
## was rescaled for plotting purposes (equivalent to adjusting v3, the 
## volume of the biophase compartment; I used 19 l to make it look similar)
## pancuronium bromide has molecuar mass 572.861 g/mol
## 1 g = 1/572.861*1e6 umol
# g2umol <- function(x) x/572.861*1e6
# g2nmol <- function(x) x/572.861*1e9
# cm3 <- CM3(4, k10=0.770, k12=3.156, k21=3.241, k13=13.09)
# times <- seq(0.01,4,length.out=200)
# plt <- t(sapply(times, cm3)) %*% diag(1/c(19.01,18.51,19))
# plt <- log(g2nmol(plt))
# plot(times, plt[,1],
#      ylim=range(plt),
#      xlab="time (h)",
#      ylab="Conc. (ln umol/l)",
#      type="l")
# lines(times, plt[,2], lty=2)
# lines(times, plt[,3], lty=3)
# legend("topright", c("blood", "tissue", "biopahse"),
#        lty=1:3, bty='n')

## default PK parameters for our model
## assume that k12=k21 and k13=k31
## don't need volume parameters
k10_d <- 0.770
k12_d <- 3.156
k21_d <- 3.156
k13_d <- 13.090

## list of default parameters
par_d <- c(lk10=log(k10_d),
           lk12=log(k12_d),
           lk21=log(k21_d),
           lk13=log(k13_d))

## default dosing schedule
ivt_d <- list(list(time=0, dose=4),
              list(time=2, dose=4))

## initial compartment concentrations
init_d <- c(0,0,0)
## maximum simulation time (h)
tmax_d <- 16

pk_solution <- function(k10=k10_d, k12=k12_d, k21=k21_d, k13=k13_d,
                        ivt=ivt_d, init=init_d, tmax=tmax_d) {
  rits <- list()
  rits[[1]] <- ivt[[1]]
  rits[[1]]$init <- init
  rits[[1]]$begin <- ivt[[1]]$time
  if(length(ivt) > 1) {
    rits[[1]]$end <- ivt[[2]]$time
  } else {
    rits[[1]]$end <- Inf
  }
  rits[[1]]$sols <- CM3(dose=ivt[[1]]$dose,
                        k10=k10,k12=k12,k21=k21,k13=k13,
                        init=init)
  if(length(ivt) > 1) {
    for(i in 2:length(ivt)) {
      rits[[i]] <- ivt[[i]]
      rits[[i]]$init <- drop(rits[[i-1]]$sols(rits[[i]]$time - rits[[i-1]]$time))
      rits[[i]]$begin <- ivt[[i]]$time
      if(length(ivt) > i) {
        rits[[i]]$end <- ivt[[i+1]]$time
      } else {
        rits[[i]]$end <- Inf
      }
      rits[[i]]$sols <- CM3(dose=ivt[[i]]$dose,
                            k10=k10,k12=k12,k21=k21,k13=k13,
                            init=rits[[i]]$init)
    }
  }
  
  function(tms) {
    sapply(tms, function(t) {
      if(t < 0)
        return(init)
      for(rit in rits) {
        if(t >= rit$begin && t < rit$end) {
          return(rit$sols(t-rit$begin))
        }
      }
      return(NA)
    })
  }
}


## The code in this file implements ordinal logistic regression (OLR)
## in a Bayesian framework. The posterior distribution is summarized
## using the posterior mode and Hessian matrix at the mode. These
## values form the basis of a Laplace approximation to the posterior.
## This approximation is used to characterize the sampling 
## uncertainty about model predictions.

library("psych") ## cohen.kappa

# original parameter prior for vecuronium:
# Normal(vec_prior_mean, vec_prior_var)
vec_prior_mean <- c("lk10"=log(0.770)+log(19.01),
                    "lk12"=log(3.156)+log(19.01),
                    "lk21"=log(3.156)+log(19.01),
                    "lk13"=1+log(19.01),
                    "a1"=-13.65,
                    "la2"=log(4),
                    "b"=-0.0959*1e5) #needs to be scaled up if concentration values are scaled down
vec_prior_var  <- diag(rep(c(40,10,10,40,40,4000),c(1,2,1,1,1,1))) 

## the mean and covariance below come from the first round of fitting
# vec_prior_mean <- 
#   structure(c(-1.627, 1.299, 1.281, 1.094, 0.881, 1.596, -15.823, 
#               0.77, -0.824), .Names = c("lv1", "lk10", "lk12", "lk21", "lk13", 
#                                         "lk31", "a1", "la2", "b"))
# vec_prior_var <- 
#   structure(c(0.109, -0.181, 0.017, 0.009, -0.023, -0.127, 0.078, 
#               -0.012, 0.187, -0.181, 1.493, -0.221, 0.273, 0.119, 0.327, -0.002, 
#               0.057, -0.412, 0.017, -0.221, 0.379, -0.213, 0.076, 0.052, -0.114, 
#               -0.032, 0.017, 0.009, 0.273, -0.213, 0.655, -0.103, -0.12, 0.127, 
#               0.005, 0.02, -0.023, 0.119, 0.076, -0.103, 0.24, 0.098, 0.002, 
#               0.026, -0.03, -0.127, 0.327, 0.052, -0.12, 0.098, 0.445, -0.07, 
#               0.027, -0.27, 0.078, -0.002, -0.114, 0.127, 0.002, -0.07, 8.529, 
#               -1.062, 0.213, -0.012, 0.057, -0.032, 0.005, 0.026, 0.027, -1.062, 
#               0.377, -0.024, 0.187, -0.412, 0.017, 0.02, -0.03, -0.27, 0.213, 
#               -0.024, 0.364), .Dim = c(9L, 9L),
#             .Dimnames = list(c("lv1", "lk10", "lk12", "lk21", "lk13", "lk31",
#                                "a1", "la2", "b"),
#                              c("lv1", "lk10", "lk12", "lk21", "lk13", "lk31",
#                                "a1", "la2", "b")))

## logit and expit functions
logit <- function(p) log(1/(1/p-1))
expit <- function(x) 1/(1/exp(x) + 1)

## compute the probability of each outcome using parameters
## 'par' for the times (h) and case specified
olr_predict <- function(times, par, case) {
  
  ## PK 'biophase' model
  pk_sol <- 
    pk_solution(k10=exp(par["lk10"]),
                k12=exp(par["lk12"]),
                k21=exp(par["lk21"]),
                k13=exp(par["lk13"]),
                ivt=case$ivt, tmax=case$tmax)
  x <- pk_sol(times)[3,]*10000 ## (mg) # or multiply here by a bigger #
  vv <- pk_sol(times)
  
  ## POLR model
  ylev <- levels(case$tof$TOFValue)
  nlev <- length(ylev)
  yalp <- c(-Inf, par["a1"] + exp(par["la2"])*1:(nlev-1), Inf)
  xbet <- matrix(c(rep(0, length(x)),
                   rep(x * par["b"] , nlev-1), #try multiplying x here by a big #
                   rep(0, length(x))), length(x), nlev+1)
  
  prbs <- sapply(1:nlev, function(lev) {
    yunc <- rep(lev, length(times))
    expit(yalp[yunc+1] - xbet[cbind(1:length(times),yunc+1)]) - 
      expit(yalp[yunc]   - xbet[cbind(1:length(times),yunc)])
  })
  
  prbs <- matrix(prbs, ncol=nlev)
  rownames(prbs) <- NULL
  colnames(prbs) <- ylev
  attr(prbs, 'x') <- x
  attr(prbs, 'sol') <- vv
  return(prbs)
}

## classify predictions using a loss spending function
olr_classify <- function(times, par, case){
  prbs <- olr_predict(times, par, case)
  wts <- matrix(c(0:4, .5, 0:3, 1, .5, 0:2, 1.5, 1, .5, 0:1, 2, 1.5, 1, .5, 0), 
                nrow=5, byrow = TRUE)
  loss <- prbs %*% wts
  ylev <- 0:4
  choices <- factor(apply(loss, 1, which.min) - 1, levels = ylev, ordered = TRUE)
  return(list(choices, prbs))
}


## compute the probabilities associated with the observed
## TOF counts, given the model parameters
olr_probs <- function(par, case) {
  prbs <- olr_predict(case$tof$RelTOFTime/60, par, case)
  yunc <- unclass(case$tof$TOFValue)
  prbs[cbind(1:length(yunc),yunc)]
}

## compute the (weighted) log likelihood
olr_llikelihood  <- function(par, case, wts=NULL) {
  prbs <- olr_probs(par, case)
  if(is.null(wts)) 
    wts <- rep(1, length(prbs))
  sum(wts*log(prbs))
}

## compute gradients of the 'biophase' concentrations
## with respect to model parameters
olr_pkgrads <- function(par, case) {
  exr <- quote(attr(olr_predict(
    case$tof$RelTOFTime/60, 
    c(lk10=lk10, lk12=lk12, lk13=lk13,
      a1=a1, la2=la2, b=b),
    case=case), "x"))
  the <- names(par)
  env <- list2env(as.list(par))
  env$case <- case
  del <- numericDeriv(exr, the, env)
  colnames(attr(del, "gradient")) <- the
  return(del)
}

## plot gradients of the 'biophase' concentrations
## with respect to model parameters
plot_olr_pkgrads <- function(del) {
  grd <- attr(del, "gradient")
  par(mfrow=c(ncol(grd)+1, 1), mai=rep(0,4), omi=rep(1,4))
  plot(del, type="l", xaxt="n", xlab="", ylab="")
  for(i in 1:ncol(grd)) {
    plot(grd[,i], type="l",
         xaxt="n", yaxt="n",
         xlab="", ylab="")
    mtext(colnames(grd)[i], 2, 1)
  }
}

## compute the log prior density function
olr_lprior <- 
  function(par, par_mean=vec_prior_mean, par_var=vec_prior_var) {
    par_diff <- par - par_mean
    -1/2*drop(solve(par_var, par_diff) %*% par_diff)
  }

## compute the log posterior density (objective) function
olr_objective <- function(par, case, wts=NULL) {
  llik <- try(olr_llikelihood(par, case, wts), silent=TRUE)
  if(inherits(llik, "try-error"))
    llik <- -Inf
  llik + olr_lprior(par)
}

## fit the POLR model to a case
olr_fit <- function(case, wts=NULL)
  optim(vec_prior_mean,
        olr_objective, case=case, wts=wts,
        control=list(maxit=10000, fnscale=-1),
        method="BFGS", hessian=TRUE)


## compute weights for weighted fitting
## weight TOFValues 1:3 more than 0 or 4
olr_wts <- function(case, ratio=3)
  ifelse(case$tof$TOFValue %in% 1:3, ratio, 1)

## compute weighted kappa
olr_kappa <- function(par, case) {
  
  tof <- case$tof
  if(!is.null(case$tof_test))
    tof <- case$tof_test
  xp <- olr_classify(tof$RelTOFTime/60, par, case)
  xo <- tof$TOFValue
  levs <- as.numeric(union(levels(droplevels(xo)),
                           levels(droplevels(xp))))
  nlev <- length(levs)
  wts <- abs(levs %o% rep(1,nlev) - rep(1,nlev) %o% levs)
  list(n=length(xp), xo=xo, xp=xp,
       agree=mean(xp==xo),
       agree1=mean(abs(unclass(xp) - unclass(xo)) <= 1),
       kappa=cohen.kappa(x=cbind(xo,xp), w=wts))
}

setwd("~/Desktop/Dropbox/TOFPredictClean")


## edit format_case function to subset based on time of last dose
format_case_last_dose <- function(ID) {
  demo_sub <- subset(demo_data, PatientID == ID)
  meds_sub <- subset(meds_data, PatientID == ID)
  tof_sub <- subset(tof_data, PatientID == ID)
  first_dose_time <- min(meds_sub$DoseTime)
  
  meds_sub$RelDoseTime <- 
    difftime(meds_sub$DoseTime,
             first_dose_time, units="mins")
  
  tof_sub$RelTOFTime <- 
    difftime(tof_sub$TOFTime,
             first_dose_time, units="mins")
  
  tof_sub$TOFValue <- 
    factor(tof_sub$TOFValue,
           levels=c(0:4), ordered=TRUE)
  
  demo_sub$RelPACU_ArrivalTime <- 
    difftime(demo_sub$PACU_ArrivalTime,
             first_dose_time, units="mins")
  
  demo_sub$RelExtubationTime <- 
    difftime(demo_sub$ExtubationTime,
             first_dose_time, units="mins")
  
  ## if neostigmine was given, subset to
  ## records that occurred before the 
  ## first administration
  if(any(meds_sub$DrugAbbr == "Neo")) {
    neo_idx <- which(meds_sub$DrugAbbr == "Neo")
    first_neo_time <- min(meds_sub$RelDoseTime[neo_idx])
    meds_sub <- subset(meds_sub, RelDoseTime < first_neo_time)
    tof_sub  <- subset(tof_sub,  RelTOFTime < first_neo_time)
  }
  
  ## if multiple doses of the same drug were given
  ## at the same time, collapse into one dose
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
  
  
  ## assume that TOF count is 4 just before first dose
  tof_tmp <- tof_sub[1,]
  tof_tmp$RelTOFTime <- as.difftime(0, units="mins")
  tof_tmp$TOFValue   <- factor(4, levels=c(0:4), ordered=TRUE)
  tof_tmp$TOFTime    <- first_dose_time
  tof_sub <- rbind(tof_tmp, tof_sub)
  
  ## sort TOF records by RelTOFTime
  tof_ord <- order(tof_sub$RelTOFTime)
  tof_sub <- tof_sub[tof_ord,]
  
  ## subset TOF records into before and after last dose
  last_dose <- unclass(meds_sub$RelDoseTime[nrow(meds_sub)])[1]
  tof_trn <- tof_sub[which(unclass(tof_sub$RelTOFTime) < last_dose),]
  tof_tst <- tof_sub[which(unclass(tof_sub$RelTOFTime) >= last_dose),]
  
  ivt <- list()
  for(i in 1:nrow(meds_sub)) {
    if(meds_sub$DrugAbbr[i] %in% 
       c("Vec", "Roc", "Atr", "Cis")) {
      ivt[[length(ivt) + 1]] <- 
        list(time=as.numeric(meds_sub$RelDoseTime[i]/60),
             dose=meds_sub$DoseAmount[i])
    }
  }
  
  rel_max_time <- max(c(tof_sub$RelTOFTime,
                        meds_sub$RelDoseTime))
  list(demo=demo_sub, meds=meds_sub, tof=tof_trn,
       tof_test=tof_tst, ivt=ivt, tmax=rel_max_time/60)
}

# source("R/clinical-data.R")
# ## subset to patients with at least 2 doses, at least one TOF count < 4, at least one TOF count after final dose, at least 5 training points; remove observation 626 because training data is only 0, 1 and after is 2-4
# NMB_one_vec_tof_ne4_ge2_IDs <- 
#   NMB_one_vec_tof_IDs[sapply(vec_case_list, 
#                              function(x) any(x$tof$TOFValue < 4) 
#                              & (nrow(x$meds) > 1) 
#                              & (unclass(x$meds$RelDoseTime[nrow(x$meds)] < unclass(x$tof$RelTOFTime[nrow(x$tof)]))) 
#                              & (x$tof$PatientID[1] != 626)
#                              & (length(which(x$tof$RelTOFTime < unclass(x$meds$RelDoseTime[length(x$meds$RelDoseTime)])[1])) > 4)
#   )]


## repeat steps of clinical-analysis but using edited format_case function and new patient subset
## fit model to all TOF counts before final dose
## for cases with at least 2 doses on record
par_opt_test_list_last_dose_file <- "data/par_opt_test_list_last_dose.RData"
if(!file.exists(par_opt_test_list_last_dose_file)) {
  par_opt_test_list_last_dose <- list()
  for(id in NMB_one_vec_tof_ne4_ge2_IDs) {
    cat(paste0("ID: ", id, "..."))
    idx <- length(par_opt_test_list_last_dose) + 1
    
    if(!any(tof_data$PatientID == id)) {
      cat("no TOF data\n")
      par_opt_test_list_last_dose[[idx]]      <- list()
      par_opt_test_list_last_dose[[idx]]$id   <- id
      par_opt_test_list_last_dose[[idx]]$case <- format_case_last_dose(id)
      par_opt_test_list_last_dose[[idx]]$fit  <- try(stop("no TOF data\n"), silent=TRUE)
      next
    }
    
    par_opt_test_list_last_dose[[idx]]      <- list()
    par_opt_test_list_last_dose[[idx]]$id   <- id
    par_opt_test_list_last_dose[[idx]]$case <- format_case_last_dose(id)
    par_opt_test_list_last_dose[[idx]]$fit  <- 
      try(olr_fit(par_opt_test_list_last_dose[[idx]]$case))
    cat("done\n")
  }
  save(par_opt_test_list_last_dose, file=par_opt_test_list_last_dose_file,
       compress="xz", compression_level=9)
} else {
  load(par_opt_test_list_last_dose_file)
}

## note: fit does not work for ID 626, only observe 0 and 1 prior to last dose, so not predicting 2, 3, or 4, which are only values observed AFTER last dose


## create list of predicted values

predicted_list_last_dose <- lapply(par_opt_test_list_last_dose, function(x) {
  ## compute predicted probabilities for each TOF count
  if(inherits(x$fit, "try-error"))
    return(NA)
  cpred <- as.numeric(olr_classify(x$case$tof_test$RelTOFTime/60,
                                   x$fit$par, x$case)[[1]]) - 1
  
  return(cpred)
})

names(predicted_list_last_dose) <- as.character(NMB_one_vec_tof_ne4_ge2_IDs)
predicted_last_dose <- unlist(predicted_list_last_dose)

## create list of actual values

actual_list_last_dose <- lapply(par_opt_test_list_last_dose, function(x) {
  ## compute predicted probabilities for each TOF count
  if(inherits(x$fit, "try-error"))
    return(NA)
  cactl <- as.numeric(as.character(x$case$tof_test$TOFValue))
  return(cactl)
  #c(terr = mean(cpred != cactl),
  #terr1 = mean(abs(cpred - cactl) > 1))
})

names(actual_list_last_dose) <- names(predicted_list_last_dose)
actual_last_dose <- unlist(actual_list_last_dose)

## summary statistics

# number of total cases
length(par_opt_test_list_last_dose)

# number of total TOF counts in set
TOFs <- data.frame(train=rep(NA,140),test=rep(NA,140),total=rep(NA,140))
for(i in 1:length(par_opt_test_list_last_dose)){
  TOFs$train[i] <- nrow(par_opt_test_list_last_dose[[i]]$case$tof)
  TOFs$test[i] <- nrow(par_opt_test_list_last_dose[[i]]$case$tof_test)
  TOFs$total[i] <- TOFs$train[i] + TOFs$test[i]
}
summary(TOFs)
apply(TOFs,2,sum)
apply(TOFs,2,mean)
apply(TOFs,2,sd)

# demographics
demos <- data.frame(ID = rep(NA,140), VecDose = rep(NA,140), NeoDose = rep(NA,140), Duration = rep(NA,140), Height = rep(NA,140), Weight = rep(NA,140), BMI = rep(NA,140), Age = rep(NA,140), Sex = rep(NA,140), TOFs.train = rep(NA,140), TOFs.test = rep(NA,140), TOFs.total = rep(NA,140))

for(i in 1:length(par_opt_test_list_last_dose)){
  demos$ID[i] <- par_opt_test_list_last_dose[[i]]$case$demo$PatientID
  demos$VecDose[i] <- as.numeric(as.character(par_opt_test_list_last_dose[[i]]$case$demo$VecDoseTotal))
  demos$NeoDose[i] <- as.numeric(as.character(par_opt_test_list_last_dose[[i]]$case$demo$NeoDoseTotal))
  demos$Duration[i] <- as.numeric(as.character(par_opt_test_list_last_dose[[i]]$case$demo$SurgeryDuration))
  demos$Height[i] <- as.numeric(as.character(par_opt_test_list_last_dose[[i]]$case$demo$Height))
  demos$Weight[i] <- as.numeric(as.character(par_opt_test_list_last_dose[[i]]$case$demo$Weight))
  demos$BMI[i] <- as.numeric(as.character(par_opt_test_list_last_dose[[i]]$case$demo$BMI))
  demos$Age[i] <- as.numeric(as.character(par_opt_test_list_last_dose[[i]]$case$demo$Age))
  demos$Sex[i] <- as.numeric(as.character(par_opt_test_list_last_dose[[i]]$case$demo$Sex))
  demos$TOFs.train[i] <- nrow(par_opt_test_list_last_dose[[i]]$case$tof)
  demos$TOFs.test[i] <- nrow(par_opt_test_list_last_dose[[i]]$case$tof_test)
  demos$TOFs.total[i] <- demos$TOFs.train[i] + demos$TOFs.test[i]
}

dems <- demos[,2:12]


## how many predictions were accurate?

(val_last_dose <- mean(actual_last_dose == predicted_last_dose, na.rm=TRUE))

## prediction matrix

combos_last_dose <- rep(NA,335)
for(i in 1:335) {
  if(is.na(actual_last_dose[i])){
    combos_last_dose[i] <- NA
  }
  else if(actual_last_dose[i]==0) {
    combos_last_dose[i] <- letters[predicted_last_dose[i]+1]
  }
  else if(actual_last_dose[i]==1) {
    combos_last_dose[i] <- letters[predicted_last_dose[i]+6]
  }
  else if(actual_last_dose[i]==2) {
    combos_last_dose[i] <- letters[predicted_last_dose[i]+11]
  }
  else if(actual_last_dose[i]==3) {
    combos_last_dose[i] <- letters[predicted_last_dose[i]+16]
  }
  else if(actual_last_dose[i]==4) {
    combos_last_dose[i] <- letters[predicted_last_dose[i]+21]
  }
}

combos_last_dose <- combos_last_dose[!is.na(combos_last_dose)]
to.add <- which(!(letters[1:25] %in% combos_last_dose))
combos_table <- c(table(combos_last_dose), rep(0,length(to.add)))
#names(combos_table)[(25-length(to.add)+1):25] <- letters[to.add]
combos_table <- combos_table[sort(names(combos_table))]


combo_mat_last_dose <- matrix(combos_table,nrow=5,byrow = TRUE)
colnames(combo_mat_last_dose) <- c("0","1","2","3","4")
rownames(combo_mat_last_dose) <- c("0","1","2","3","4")

## table of percentages to show relative accuracy at each combo
combo_props_last_dose <- round(prop.table(combo_mat_last_dose,margin=2)*100,2)


## frequency of each absolute value disparity

(prop_diff_last_dose <- 100*round(prop.table(table((predicted_last_dose-actual_last_dose))),4))

## percent within 1
(within_1 <- sum(prop_diff_last_dose[4:6]))



