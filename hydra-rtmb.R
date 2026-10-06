library(RTMB)
library(tidyverse)

hydra_data <- readRDS("inputs/hydra_sim_GB_new5bin_1978_3fleets_4surveys.rds")

#data list object
data <- list()

#model dimensions
data$Nsp <- hydra_data$Nspecies
data$Nlen <- hydra_data$Nsizebins
data$Nstep <- 1
data$Nflt <- hydra_data$Nfleets
data$Nsurvey <- 4 #hydra_data$Nsurveys
data$Nyr <- hydra_data$Nyrs
data$Ntime <- hydra_data$Nyr * data$Nstep
data$t_yrs <- rep(1:(data$Nyr/data$Nstep), each = data$Nstep)

Ntime <- data$Ntime
Nsp <- data$Nsp
Nyr <- data$Nyr
Nflt <- data$Nflt
Nlen <- data$Nlen
Nsurvey <- data$Nsurvey
Ntarget <- Nsp*Nflt #placeholder

#recruitment
data$rec_type <- hydra_data$recType. # 4 = Ricker, 9 = avg rec.

#fishery targeting
fishind <- hydra_data$indicatorFisheryq
data$fish_primary <- apply(fishind, 1, function(x) which(x==1)[1])
data$q_map <- map_dfr(apply(fishind, 1, function(x) which(x==1)[-1]), as_tibble,.id = "flt") |> mutate(flt = as.numeric(as.factor(flt))) |>
  as.matrix()

#growth
data$lbin_hi <- t(apply(hydra_data$binwidth, 1, cumsum))
data$lbin_lo <- cbind(rep(0,Nsp), data$lbin_hi[,-Nlen])
data$lbin_mid <- data$lbin_lo + (data$lbin_hi - data$lbin_lo)/2
data$weight_mid <- apply(data$lbin_mid, 2, FUN = function(x) {hydra_data$lenwta*x^hydra_data$lenwtb})

#consumption
data$stomach_wt <- array(rep(hydra_data$intakeStomach,Ntime),dim=c(Nsp,Nlen,Ntime))
#object contains old values for stomach weights - needs updating!
data$intake_a <- hydra_data$intakeAlpha
data$intake_b <- hydra_data$intakeBeta

#covariates
data$temp_consumption <- hydra_data$observedTemperature[2,] #could go in consumption block
data$covar_growth <- hydra_data$observedTemperature[2,]
data$covar_maturity <- hydra_data$observedTemperature[2,]
data$covar_recruitment <- hydra_data$observedTemperature[2,]

#observations for objective function
data$obs_surveyindex <- as.data.frame(hydra_data$observedBiomass[,c(3,1,2,4,5)])
data$obs_surveysize <- as.data.frame(hydra_data$observedSurvSize[,c(3,1,2,5:ncol(hydra_data$observedSurvSize))])
data$obs_surveydietcomp <- as.data.frame(hydra_data$observedSurvDiet[,c(3,1,2,4:ncol(hydra_data$observedSurvDiet))])
data$obs_catch <- as.data.frame(hydra_data$observedCatch[,c(4,1,3,5,6)])
data$obs_catchsize <- as.data.frame(hydra_data$observedCatchSize[,c(4,1,3,6:ncol(hydra_data$observedCatchSize))])



#parameter list object
parameters <- list()

parameters$log_F_flt <- array(rep(log(0.2),Nyr*Nflt),dim = c(Nflt, Nyr)) #probably derived
parameters$log_fishsel_par <- matrix(c(rep(4,Nflt),rep(2,Nflt)), byrow = TRUE, nrow = 2)
parameters$log_fish_q <- rep(0, Ntarget - Nflt)

parameters$logit_vulnerability <- matrix(rep(5,Nsp*Nsp),nrow = Nsp) #vuln pars would be number of positives in isprey matrix
parameters$log_sigma_predpref <- log(hydra_data$M2sizePrefSigma) #rep(1,Nsp)
parameters$log_psi_pref <- log(hydra_data$M2sizePrefMu) #rep(log(0.5), Nsp)
parameters$log_intake_a <- log(hydra_data$intakeAlpha) #rep(log(0.2), Nsp)
parameters$log_intake_b <- log(hydra_data$intakeBeta)  #rep(log(0.2), Nsp)

parameters$log_vb_Linf <- log(as.numeric(hydra_data$growthLinf))
parameters$beta_growth <- rep(0, Nsp)
parameters$log_vb_k  <- rep(log(as.numeric(hydra_data$growthK)), Nsp)

parameters$log_mat_a <- rep(log(20),Nsp)
parameters$log_mat_b <- rep(log(4), Nsp)
parameters$beta_maturity <- rep(0, Nsp)

parameters$log_initN <- matrix(rep(4,Nsp*Nlen),byrow=TRUE,nrow = Nsp)

parameters$log_ricker_alpha <- rep(log(4e-06), Nsp)
parameters$log_ricker_beta <- rep(log(1), Nsp)
parameters$beta_recruits <- rep(0, Nsp)
parameters$log_avg_recruit <- rep(5, Nsp)
parameters$log_sigma_r <- rep(log(0.5), Nsp)
parameters$rec_devs <- matrix(rep(0,Nsp*Nyr),byrow=TRUE,nrow = Nsp)

parameters$log_M1 <- rep(log(0.2),Nsp)
parameters$log_other_food <- rep(14, Nsp)

parameters$log_survsel_par <- matrix(c(rep(4,Nsurvey),rep(2,Nsurvey)), byrow = TRUE, nrow = 2)
parameters$log_surv_q <- matrix(rep(0,Nsp*Nsurvey), byrow=TRUE, nrow = Nsurvey)


baby_hydra_fun <- function(data, parameters) {
  "[<-" <- ADoverload("[<-")

  #make things visible
  getAll(data, parameters, warn=FALSE)
# procedure section structure

#varibales
neglog19 <- -1*log(19)

#transform parameters
## recruitment
sigma_r <- exp(log_sigma_r)
ricker_alpha <- exp(log_ricker_alpha)
ricker_beta <- exp(log_ricker_beta)
avg_recruit <- exp(log_avg_recruit)
## maturity
mat_a <- exp(log_mat_a)
mat_b <- exp(log_mat_b)
## growth
vb_k <- exp(log_vb_k)
vb_Linf <- exp(log_vb_Linf)
## M
M1_ann <- exp(log_M1) #hydra code doing this slightly differently
other_food <- array(rep(exp(log_other_food),Ntime),dim=c(Nsp,Ntime))
vulnerability <- exp(logit_vulnerability)/(1+exp(logit_vulnerability))
sigma_predpref <- exp(log_sigma_predpref)
psi_pref <- exp(log_psi_pref)
intake_a <- exp(log_intake_a)
intake_b <- exp(log_intake_b)
##surveys
survsel_par <- exp(log_survsel_par)
surv_q <- exp(log_surv_q)
##fishing
F_flt <- exp(log_F_flt)
fishsel_par <- exp(log_fishsel_par)
fish_q <- exp(log_fish_q)


#GF calc other things not dependent on states
#fishing
## selectivity
fish_sel <- array(1, dim=c(Nflt, Nsp, Nlen))
for (iflt in 1:Nflt) {
    fish_sel[iflt,,] = 1/(1+exp(neglog19*(lbin_mid-fishsel_par[1,iflt])/fishsel_par[2,iflt]))
    for (isp in 1:Nsp) {
      fish_sel[iflt,isp,] <- fish_sel[iflt,isp,]/max(fish_sel[iflt,isp,])
    }
}
## q
fishery_q <- array(0, dim = c(Nflt, Nsp))
for (iflt in 1:Nflt) fishery_q[iflt,fish_primary[iflt]] <- 1
for (irow in 1:nrow(q_map))
  fishery_q[q_map[irow,1], q_map[irow,2]] <- fish_q[irow]
## F
FF <- array(0,dim=c(Nsp,Nflt,Nyr))
for (isp in 1:Nsp) {
  for (iflt in 1:Nflt) {
    FF[isp,iflt,] <- fishery_q[iflt, isp]*F_flt[iflt,]
  }
}
#sum over fleets. (GF can make this more efficient later)
t <- 0
Frate <- array(0, dim=c(Nsp,Nlen,Ntime))
for (iyr in 1:Nyr) {
  for (istep in 1:Nstep) {
    t <- t+1
    for (isp in 1:Nsp) {
      for (ilen in 1:Nlen) {
        Frate[isp, ilen, t] <- sum(fish_sel[, isp, ilen]*FF[isp, ,iyr]/Nstep)
      }
    }
  }
}


#surveys
survey_sel <- array(1, dim=c(Nsurvey, Nsp, Nlen))
for (isurv in 1:Nsurvey) {
    survey_sel[isurv,,] = 1/(1+exp(neglog19*(lbin_mid-survsel_par[1,isurv])/survsel_par[2,isurv]))
    for (isp in 1:Nsp) {
      survey_sel[isurv, isp, ] <- survey_sel[isurv, isp, ]/max(survey_sel[isurv, isp, ])
    }
}


## M1
M1 <- M1_ann/Nstep

## other food

## suitability etc.
preference <- array(0, dim=c(Nsp,Nlen,Nsp,Nlen))
suitability <- array(0, dim=c(Nsp,Nlen,Nsp,Nlen))
daily_intake <- array(0, dim=c(Nsp,Nlen,Ntime))
for (predsp in 1:Nsp) {
  for (predlen in 1:Nlen) {
    for (preysp in 1:Nsp) {
      for (preylen in 1:Nlen) {
    #size pref (equation 13)
    preference[predsp, predlen, preysp, preylen] <- (1/(
      (weight_mid[preysp, preylen]/weight_mid[predsp, predlen])*sigma_predpref[predsp]*sqrt(2*pi)))*
      exp(-1*(log(weight_mid[preysp, preylen]/weight_mid[predsp, predlen])-psi_pref[predsp])^2/(2*sigma_predpref[predsp]^2))
    #suitability (equation 12)
    suitability[predsp, predlen, preysp, preylen] <- preference[predsp, predlen, preysp, preylen] * vulnerability[predsp, preysp]
    # intake (equation 14)
    daily_intake[predsp, predlen, ] <- 24*(intake_a[predsp]*exp(intake_b[predsp]*temp_consumption)*
                                             stomach_wt[predsp, predlen, ])
      }
    }
  }
}


## growth
#VB growth (equation 7)
phi <- array(0,dim=c(Nsp,Nlen,Ntime))
for (isp in 1:Nsp) {
  for (ilen in 1:Nlen) {
    for (t in 1:Ntime) {
    numerator <- vb_Linf[isp]*exp(beta_growth[isp]*covar_growth[t])-lbin_lo[isp,ilen]
    denominator <- vb_Linf[isp]*exp(beta_growth[isp]*covar_growth[t])-lbin_hi[isp,ilen]

    if (numerator <0 | denominator < 0) {
      phi[isp, ilen, t] <- 0
    }
    else
    {
      phi[isp, ilen, t] <- vb_k[isp]/log(numerator/denominator)
    }
    }
  }
}
## maturity
#maturity (equation 9)
maturity <- array(0, dim = c(Nsp,Nlen,Ntime))
for (isp in 1:Nsp) {
  for (ilen in 1:Nlen) {
    maturity[isp, ilen, t] <- 1/(1 + exp(-1*(mat_a[isp]+mat_b[isp]*lbin_mid[isp, ilen])+beta_maturity[isp]*covar_maturity[t]))
  }
}

# storage variables for quantities that are functions of the states
N <- array(0, dim = c(Nsp,Nlen,Ntime))
S <- array(0, dim = c(Nsp,Nlen,Ntime))
recruits <- array(0, dim = c(Nsp, Nyr))
M2 <- array(0, dim = c(Nsp,Nlen,Ntime))
total_food <- array(0, dim = c(Nsp,Nlen,Ntime))
N_tot <- array(0, dim = c(Nsp,Nlen,Nyr))
ssb <- array(0, dim = c(Nsp, Nyr))

Zflt <- array(0, dim = c(Nsp, Nflt, Nlen, Ntime))
pred_catchsize_bio <- array(0, dim = c(Nsp, Nflt, Nlen, Ntime))
pred_catchsize_num <- array(0, dim = c(Nsp, Nflt, Nlen, Ntime))
pred_catch <- array(0, dim = c(Nsp, Nflt, Ntime))
pred_catchsize_bio_yr <- array(0, dim = c(Nsp, Nflt, Nlen, Nyr))
pred_catchsize_num_yr <- array(0, dim = c(Nsp, Nflt, Nlen, Nyr))
pred_catch_yr <- array(0, dim = c(Nsp, Nflt, Nyr))
pred_survey_index_yr <- array(0, dim = c(Nsp, Nsurvey, Nyr))
pred_surveysize_num_yr <- array(0, dim = c(Nsp, Nsurvey, Nlen, Nyr))


#calc initial states
N[,,1] <- exp(log_initN)

#loop over time
for (t in 1:Ntime) {
#t <- 1

##year
iyr <- t_yrs[t]

## increment Nvec
if (t > 1) {
  N[,,t] <- N[,,t-1]
}

##recruitment if not t = 1
if (t > 1) {
  #add recruits if first time step of year
  if ((t %% Nstep)-1 == 0) {
    for (isp in 1:Nsp) {

      ##Ricker (equation 10)
      if (rec_type[isp == 4]) {
        recruits[isp, iyr] <- ricker_alpha[isp]*ssb[isp, iyr - rec_lag[isp]] *
          exp(-1*ricker_beta[isp]*ssb[isp, iyr] + beta_recruits[isp]*covar_recruitment[iyr])
      }

      ##Average plus devs (equation 8)
      if (rec_type[isp == 9]) {
        recruits[isp, iyr] <- avg_recruit[isp] * exp(beta_recruits[isp]*covar_recruitment[iyr])
      }

      #add rec devs
      recruits[isp, iyr] <- recruits[isp, iyr]*exp(rec_devs[isp, iyr] - 0.5*sigma_r[isp]*sigma_r[isp])

      #if ((t %% Nstep)-1 == 0)
      N[isp, 1, t] <- N[isp, 1, t] + recruits[isp, iyr]
    }
  }
}

##calc available N? [secondary, for sp not always in area]
##predation mortality
#predation mortality (equation 15)

for (predsp in 1:Nsp) {
  for (predlen in 1:Nlen) {
    total_food[predsp, predlen, t] <- sum(suitability[predsp, predlen, , ]*weight_mid*N[,,t]) +
      other_food[predsp, t]
  }
}

M2[ , , t]  <- 0
for (preysp in 1:Nsp) {
  for (preylen in 1:Nlen) {
    for (predsp in 1:Nsp) {
      for (predlen in 1:Nlen) {
        M2[preysp, preylen, t] <- M2[preysp, preylen, t] +
          daily_intake[predsp, predlen, t] * (365/Nstep) *
          N[predsp, predlen, t] * suitability[predsp, predlen, preysp, preylen] /
          total_food[predsp, predlen, t]
      }
    }
  }
}


##total mortality
#survival equation 11
S[,,t] <- exp(-1*(M1 + M2[,,t] + Frate[,,t]))

##pop_dynamics

#update numbers at age (equation 1)
for (isp in 1:Nsp) {
  Ntemp <- N[isp, , t]
  for (ilen in Nlen:2) {
    Ntemp[ilen] <-  phi[isp, ilen-1, t]*S[isp, ilen-1, t]*N[isp, ilen-1, t] +
      (1 - phi[isp, ilen, t])*S[isp, ilen, t]*N[isp, ilen, t]
  }
  Ntemp[1] <- (1 - phi[isp, 1, t])*S[isp, 1, t]*N[isp, 1, t]
  N[isp,,t] <- Ntemp
  N_tot[isp,,iyr] <- N_tot[isp,,iyr] + N[isp,,t]
}

##SSB
ssb[isp, iyr] <- sum(maturity[isp, , t]*weight_mid[isp,]*N[isp,,t],na.rm = TRUE)

##movement

# model fit
##predicted values for data

#predicted catch (equation 17). [check this - should denominator be all fleets/sizes?]
for (isp in 1:Nsp) {
  for (iflt in 1:Nflt) {
    Zflt[isp,iflt,,t] <- M1[isp] + M2[isp, , t] + fish_sel[iflt, isp,]*FF[isp,iflt,t]
    pred_catchsize_bio[isp, iflt, , t] <- fish_sel[iflt,isp,]*FF[isp,iflt,t]*N[isp,,t]*weight_mid[isp,]/Zflt[isp,iflt,,t]
    pred_catchsize_num[isp, iflt, , t] <- fish_sel[iflt,isp,]*FF[isp,iflt,t]*N[isp,,t]/Zflt[isp,iflt,,t]
    pred_catch[isp, iflt, t] <- sum(pred_catchsize_bio[isp, iflt, , t])
    pred_catchsize_bio_yr[isp, iflt, , iyr] <- pred_catchsize_bio_yr[isp, iflt, , iyr] + pred_catchsize_bio[isp, iflt, , t]
    pred_catchsize_num_yr[isp, iflt, , iyr] <- pred_catchsize_num_yr[isp, iflt, , iyr] + pred_catchsize_num[isp, iflt, , t]
    pred_catch_yr[isp, iflt, iyr] <- pred_catch_yr[isp, iflt, iyr] + pred_catch[isp, iflt, t]
  }
}

}
#create annual quantities

#predicted values for survey quantities
for (isp in 1:Nsp) {
  for (isurv in 1:Nsurvey) {
    for (iyr in 1:Nyr) {
      pred_survey_index_yr[isp, isurv, iyr] <- sum(N_tot[isp, , iyr] * survey_sel[isurv, isp, ] * surv_q[isurv, isp] * weight_mid[isp, ])/ Nstep #in biomass
      pred_surveysize_num_yr[isp, isurv, , iyr] <- N_tot[isp, , iyr] * survey_sel[isurv, isp, ]
      #pred_surveysize_num_yr[isp, isurv, , iyr] <- pred_surveysize_num_yr[isp, isurv, , iyr]/sum(pred_surveysize_num_yr[isp, isurv, , iyr])
      #for (ipredlen in 1:Nsp) {
      #  pred_surveydiet_wt_yr[isp,isurv,ipredlen, ,iyr] <-
      #}
    }
  }
}

##objective function
nll <- 0

#survey index
pred_vec <- pred_survey_index_yr[obs_surveyindex[,1], obs_surveyindex[,2], obs_surveyindex[,3]] + 0.0001 #small constant?
nll <- nll - sum(dnorm(obs_surveyindex[,4], log(pred_vec), obs_surveyindex[,5], log = TRUE))

#survey size comp
for (iobs in 1:nrow(obs_surveysize)) {
  isp <- obs_surveysize[iobs,1]
  isurv <-obs_surveysize[iobs,2]
  iyr <- obs_surveysize[iobs,3]
  temp_obs <- as.numeric(obs_surveysize[iobs,4]*obs_surveysize[iobs,-(1:4)]/sum(obs_surveysize[iobs,-(1:4)]))
  temp_pred <- pred_surveysize_num_yr[isp,isurv,,iyr]/sum(pred_surveysize_num_yr[isp,isurv,,iyr])
  #print(class(temp_obs))
  nll <- nll - dmultinom(temp_obs, prob = temp_pred, log = TRUE)
}

#catch
pred_vec <- pred_catch_yr[obs_catch[,1], obs_catch[,2], obs_catch[,3]] + 0.0001
nll <- nll - sum(dnorm(obs_catch[,4], log(pred_vec), obs_catch[,5], log = TRUE))

#catch size comp
for (iobs in 1:nrow(obs_catchsize)) {
  isp <- obs_catchsize[iobs,1]
  iflt <- obs_catchsize[iobs,2]
  iyr <- obs_catchsize[iobs,3]
  if (sum(obs_catchsize[iobs,-(1:4)])>0) {
  temp_obs <- as.numeric(obs_catchsize[iobs,4]*obs_catchsize[iobs,-(1:4)]/sum(obs_catchsize[iobs,-(1:4)]))
  temp_pred <- pred_catchsize_num_yr[isp,iflt,,iyr]/sum(pred_catchsize_num_yr[isp,iflt,,iyr])
  #print(class(temp_obs))
  nll <- nll - dmultinom(temp_obs, prob = temp_pred, log = TRUE)
  }
}

# #diet composition
# for (iobs in 1:nrow(obs_surveydietcomp)) {
#   ipred <- obs_surveydietcomp[iobs,1]
#   isurv <-obs_surveydietcomp[iobs,2]
#   iyr <- obs_surveydietcomp[iobs,3]
#   ipredlen <- obs_surveydietcomp[iobs,4]
#   temp_obs <- as.numeric(obs_surveydietcomp[iobs,5]*obs_surveydietcomp[iobs,-(1:5)]/sum(obs_surveydietcomp[iobs,-(1:5)]))
#   temp_pred <- pred_surveydiet_wt_yr[isp,isurv,ipredlen,,iyr]/sum(pred_surveydiet_wt_yr[ipred,isurv,ipredlen,,iyr])
#   #print(class(temp_obs))
#   nll <- nll - dmultinom(temp_obs, prob = temp_pred, log = TRUE)
# }

#recruitment
for (isp in 1:Nsp) {
  nll <- nll - sum(dnorm(rec_devs[isp,], 0, sigma_r[isp], log = TRUE))
}


 nll
} #end baby_hydra_fun

baby_hydra_fun(data, parameters)

cmb <- function(f, d) function(p) f(p, d)

# ## Make a function object
# obj <- MakeADFun(cmb(baby_hydra_fun, data), parameters,
#                  #random = c("rdev"),
#                  map = list(logM = factor(NA),
#                             #logfracR1 = factor(NA),
#                             #logith = factor(NA),
#                             logSigmaC = factor(NA), #)) #,
#                             logSigmaI = factor(NA),
#                             logSigmaR = factor(NA)#,
#                             #recdev = factor(rep(NA,nrow(catch)))
#                  ))
# #control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))
#
# ## Call function minimizer
# opt <- nlminb(obj$par, obj$fn, obj$gr, control= list(eval.max = 10000,
#                                                      iter.max = 10000))
#
# ## Get parameter uncertainties and convergence diagnostics
# sdr <- sdreport(obj)
# summary(sdr)
#
