##########################################################################
# natural_history_code.R
#
# Author: Jane Lange
#
# Purpose:
#   Fit a continuous-time multi=state natural-history model by Lange et al. (2024)
#
#     See: 
#     Lange JM, Gogebakan KC, Gulati R, Etzioni R. 
#     Projecting the impact of multi-cancer early detection on late-stage incidence using 
#     multi-state disease modeling. Cancer Epidemiology, Biomarkers & Prevention. 
#     2024 Jun 3;33(6):830-7
#
#   
#   The fitting procedure:
#       (1) Reads SEER-style incidence by age & stage
#       (2) Builds a CTMC rate matrix consistent with OMST and DMST
#       (3) Computes age-specific cause-specific hazards (early vs late)
#       (4) Maximizes the Poisson likelihood of observed counts
#       (5) Returns fitted rate matrix, sojourn times, and diagnostic plots
#
# Inputs:
#       • Incidence data: mid-age, Count, Pop, years, stage
#       • OMST = mean sojourn time
#       • DMST = late-stage mean sojourn time
#       • k = number of states in CTMC 
#
# Outputs:
#       • Fitted rate matrix
#       • Log-likelihood optimal parameters
#       • Predicted incidence vs observed incidence plots
#       • Estimate of EMST
#
###############################################################################


#################################
# libraries
###################################
library(msm)
library(ggplot2)
library(reshape2)
library(nloptr)

# =====================================================================
# read_data()
#   Reads SEER-like age × stage incidence data.
#   Computes:
#       • Person-years   = Pop × years
#       • Rate per 100k  = Count / PY × 1e5
#   Restricts to ages 0 < midage < 80.
# =====================================================================
read_data <- function(filepath, mulitplier = 1) {
  the_data <- read.csv(filepath)
  the_data$PY   <- the_data$Pop * the_data$years
  the_data$rate <- mulitplier * the_data$Count / the_data$PY * 1e5
  the_data$Count <- the_data$Count * mulitplier
  the_data <- the_data[the_data$midage > 0 & the_data$midage < 80, ]
  return(the_data)
}

# =====================================================================
# get_rate_mat()
#
# Construct CTMC rate matrix consistent with:
#       • OMST (mst)
#       • DMST (dmst)
#
# This follows the algebraic constraints described in Lange et al.,
# ensuring that transition intensities produce the desired mean sojourn
# times in early and late preclinical periods.
#
# params = log-transition rates for early preclinical states.
# =====================================================================
get_rate_mat <- function(params, k, mst, dmst) {
  
  rate_mat <- matrix(0, nrow = k, ncol = k)
  
  # early preclinical transitions
  for (i in 1:(k - 3)) {
    rate_mat[i, i + 1] <- exp(params[i])
  }
  
  # late-stage preclinical rate structure (Lange constraints)
  c <- 1 / dmst
  b <- rate_mat[k - 3, k - 2]
  rate_mat[k - 3, k - 1] <- (c + b) / (c * mst) - b
  a <- rate_mat[k - 3, k - 1]
  rate_mat[k - 2, k] <- c
  
  # finalize diagonal
  diag(rate_mat) <- -apply(rate_mat, 1, sum)
  return(rate_mat)
}

# =====================================================================
# cause_spec_hazards()
#
# Given age & rate matrix, compute:
#       • h_L = cause-specific hazard of EARLY-stage diagnosis
#       • h_D = cause-specific hazard of LATE-stage diagnosis
#
# It uses the CTMC forward solution MatrixExp(age).
# =====================================================================
cause_spec_hazards <- function(age, rate_mat, k) {
  
  prob_mat <- MatrixExp(t = age, mat = rate_mat)
  
  prob_l <- prob_mat[1, k - 1]    # early-stage clinical
  prob_d <- prob_mat[1, k]        # late-stage clinical
  denom  <- 1 - (prob_l + prob_d)
  
  # numerator for early-stage hazard
  num_L <- prob_mat[1, k - 3] * rate_mat[k - 3, k - 1]
  h_L   <- num_L / denom
  
  # numerator for late-stage hazard
  num_D <- prob_mat[1, k - 2] * rate_mat[k - 2, k]
  h_D   <- num_D / denom
  
  return(list(h_L = h_L, h_D = h_D,
              dens_L = num_L, dens_D = num_D,
              prob_L = prob_l, prob_D = prob_d))
}

# =====================================================================
# indiv_likelihood()
#
# Poisson likelihood contribution for one age bin:
#       observed_count ~ Poisson(PY × hazard_stage)
#
# hazard_stage = h_L or h_D depending on observed stage.
# =====================================================================
indiv_likelihood <- function(age, observed, PY, stage,
                             rate_mat, k) {
  
  haz_out <- cause_spec_hazards(age = age, rate_mat = rate_mat, k = k)
  thehaz <- if (stage %in% c("local", "Early")) haz_out$h_L else haz_out$h_D
  
  mean_val <- PY * thehaz
  loglike <- observed * log(mean_val) - mean_val
  return(loglike)
}

# =====================================================================
# likelihood_fun()
#
# Wrapper to compute the total log-likelihood across all ages/stages.
# Optimized over params (log transition intensities).
# =====================================================================
likelihood_fun <- function(params, mst, dmst, k,
                           allage, allobserved, allPY, allstage) {
  
  rate_mat <- get_rate_mat(params, k, mst, dmst)
  loglike  <- mapply(allage, allobserved, allPY, allstage,
                     FUN = "indiv_likelihood",
                     MoreArgs = list(rate_mat = rate_mat, k = k))
  
  sum_loglike <- sum(unlist(loglike))
  return(sum_loglike)
}

# =====================================================================
# get_max_LL()
#
# Runs multiple random seeds and selects the best-fitting optimization.
# Uses optim(..., method="L-BFGS") with box constraints for the last rate.
# =====================================================================
get_max_LL <- function(the_data, mst, num_seeds, k, dmst, mean1) {
  
  c <- 1 / dmst
  M <- mst
  if (c < 1 / M) {
    print("Error: dmst must be less than mst")
    return()
  }
  
  out_list <- list()
  for (i in 1:num_seeds) {
    
    out_list[[i]] <- optim(
      par   = -1 * abs(rnorm(n = (k - 3), mean = mean1)),
      fn    = likelihood_fun,
      method = "L-BFGS",
      lower  = rep(-Inf, (k - 3)),
      upper  = c(rep(Inf, (k - 4)), log(c / (c * M - 1))),
      mst = mst, dmst = dmst, k = k,
      allage = the_data$midage,
      allobserved = the_data$Count,
      allPY = the_data$PY,
      allstage = the_data$stage,
      control = list(fnscale = -1, trace = 2, maxit = 20000)
    )
  }
  
  maxLL <- which.max(unlist(lapply(out_list, "[[", "value")))
  return(out_list[[maxLL]])
}

# =====================================================================
# plot_observed_expected()
#
# Plot predicted vs observed incidence curves (early & late).
# Useful diagnostics checking whether the fitted CTMC reproduces
# SEER stage-specific rates by age.
# =====================================================================
plot_observed_expected <- function(the_data, fit, mst, dmst, k) {
  
  b <- exp(fit$par[k - 3])
  c <- 1 / dmst
  a <- (c + b) / (c * mst) - b
  
  sojourn_time  <- b / (c * (a + b)) + 1 / (a + b)
  earlypreclin  <- 1 / (a + b)
  
  rate_mat <- get_rate_mat(fit$par, k = k, mst = mst, dmst = dmst)
  
  pred_haz <- mapply(
    the_data$midage[1:79],
    FUN = "cause_spec_hazards",
    MoreArgs = list(rate_mat = rate_mat, k = k)
  )
  
  the_data$pred <- c(unlist(pred_haz[2, ]) * 1e5,
                     unlist(pred_haz[1, ]) * 1e5)
  
  the_data_m <- melt(the_data, measure.vars = c("pred", "rate"))
  
  p1 <- ggplot(the_data_m, aes(x = midage, y = value)) +
    geom_point(aes(color = variable)) +
    geom_line(aes(color = variable, group = variable)) +
    facet_grid(~stage) +
    xlab("Age") + ylab("Rate per 100,000") +
    ggtitle(paste0(
      "Mean sojourn time = ", round(sojourn_time, 2),
      ", Early preclinical = ", round(earlypreclin, 2)
    ))
  
  return(list(p1 = p1, the_data = the_data_m))
}

# =====================================================================
# summarize_fit()
#
# Compute final sojourn-time outputs for the fitted model.
# =====================================================================
summarize_fit <- function(the_data, fit, mst, dmst, k) {
  
  b <- exp(fit$par[k - 3])
  c <- 1 / dmst
  a <- (c + b) / (c * mst) - b
  
  sojourn_time  <- b / (c * (a + b)) + 1 / (a + b)
  preclinearly  <- 1 / (a + b)
  
  return(list(a = a, b = b,
              sojourn_time = sojourn_time,
              mst = mst,
              latemst = dmst,
              preclinearly = preclinearly))
}

# =====================================================================
# get_fit()
#
# Full wrapper:
#       • Runs repeated optimization
#       • Produces diagnostic plot
#       • Returns fitted rate matrix + sojourn times
#
# This is the function you call elsewhere in your pipeline.
# =====================================================================
get_fit <- function(mst, dmst, the_data, num_seeds, k, mean1) {
  
  the_fit <- tryCatch(
    get_max_LL(the_data = the_data, num_seeds = num_seeds,
               mst = mst, dmst = dmst, k = k, mean1 = mean1),
    error = function(e) e
  )
  
  while (any(class(the_fit) == "error")) {
    the_fit <- tryCatch(
      get_max_LL(the_data = the_data, num_seeds = num_seeds,
                 mst = mst, dmst = dmst, k = k, mean1 = mean1),
      error = function(e) e
    )
  }
  
  outplot <- plot_observed_expected(
    the_data = the_data, fit = the_fit, mst = mst, dmst = dmst, k = k
  )
  
  rate.matrix <- get_rate_mat(params = the_fit$par, k = k, mst = mst, dmst = dmst)
  summary_out <- summarize_fit(the_data, the_fit, mst, dmst, k)
  
  return(list(
    the_fit     = the_fit,
    outplot     = outplot$p1,
    fitted_data = outplot$the_data,
    rate.matrix = rate.matrix,
    summary_out = summary_out
  ))
}



#################################
#libraries
###################################
library(msm)
library(ggplot2)
library(reshape2)
library(nloptr)
##################################
#read in the data
read_data<-function(filepath,mulitplier=1){
  the_data=read.csv(filepath)
  the_data$PY=the_data$Pop*the_data$years
  the_data$rate=mulitplier*the_data$Count/the_data$PY*1e5
  the_data$Count=the_data$Count*mulitplier
  the_data=the_data[the_data$midage > 0 & the_data$midage < 80, ]
  return(the_data)
}

get_rate_mat<-function(params,k,mst,dmst){
  # browser()
  rate_mat=matrix(0,nrow=k, ncol=k)
  for(i in 1:(k-3)){
    rate_mat[i,i+1]=exp(params[i])
  }
  c=1/dmst
  b=rate_mat[k-3,k-2]
  rate_mat[k-3,k-1]=(c+b)/(c*mst)-b
  a=rate_mat[k-3,k-1]
  rate_mat[k-2,k]=c
  print(b/(c*(a+b))+1/(a+b))
  
  diag(rate_mat)=-apply(rate_mat,1,"sum")
  return(rate_mat)
}

#Compute the cause specific hazard##
cause_spec_hazards<-function(age,rate_mat,k){
  # browser()
  prob_mat=MatrixExp(t=age,mat=rate_mat)
  prob_l=prob_mat[1,(k-1)]
  prob_d=prob_mat[1,(k)]
  denom=1-(prob_mat[1,(k-1)]+prob_mat[1,k])
  num_L=prob_mat[1,(k-3)]*rate_mat[(k-3),k-1]
  h_L=num_L/denom
  num_D=prob_mat[1,(k-2)]*rate_mat[(k-2),(k)]
  h_D=num_D/denom
  return(list(h_L=h_L,h_D=h_D,dens_L=num_L, dens_D=num_D,prob_L=prob_l,prob_D=prob_d))
}


likelihood_fun<-function(params, mst, dmst,k,
                         allage,allobserved,allPY,allstage){
  
  rate_mat=get_rate_mat(params,k, mst,dmst)
  loglike=mapply(allage,allobserved,allPY,allstage,FUN="indiv_likelihood",MoreArgs = list(rate_mat=rate_mat,k=k))
  
  # browser()
  out_LL=sum(unlist(loglike))
  return(out_LL)
}

indiv_likelihood<-function(age,observed,PY,stage,
                           rate_mat,k){
  
  haz_out=cause_spec_hazards(age=age,rate_mat=rate_mat,k=k)
  
  if(stage=="local"||stage=="Early"){  
    thehaz=haz_out$h_L
  }else{
    thehaz=haz_out$h_D
  }
  themean=PY*thehaz
  
  loglike=observed*log(themean)-themean
  return(loglike)
}


get_max_LL<-function(the_data,mst, num_seeds,k,dmst,mean1){
  c=1/dmst
  M=mst
  if(c<1/M){print("Error: dmst must be less than mst")
    return()}
  
  out_list=list()
  for(i in 1:num_seeds){
    
    out_list[[i]]=optim(par=-1*abs(rnorm(n=(k-3),mean=mean1)), fn=likelihood_fun,method="L-BFGS",lower=rep(-Inf,(k-3)),upper = c(rep(Inf, (k-4)),
                                                                                                                                 log(c/(c*M-1))),mst=mst,dmst=dmst,k=k, allage=the_data$midage,allobserved=the_data$Count,
                        allPY=the_data$PY,allstage=the_data$stage,control=list(fnscale=-1, trace=2,maxit=20000))
  }
  maxLL=which.max(unlist(lapply(out_list,"[[",c("value"))))
  maxout=out_list[[maxLL]]
  return(maxout)
}

plot_observed_expected<-function(the_data,fit,mst,dmst,k){
  b=exp(fit$par[k-3])
  c=1/dmst
  a=(c+b)/(c*mst)-b
  
  sojourn_time=b/(c*(a+b))+1/(a+b)
  earlypreclin=1/(a+b)
  the_data$stage=factor(the_data$stage)
  rate_mat=get_rate_mat(fit$par,k=k, mst=mst,dmst=dmst)
  pred_cause_spec_hazards<-mapply(the_data$midage[1:79],FUN="cause_spec_hazards",MoreArgs=list(rate_mat=rate_mat,k=k))
  the_data$pred=c(unlist(pred_cause_spec_hazards[2,])*1e5,unlist(pred_cause_spec_hazards[1,])*1e5)
  the_data=melt(the_data,measure.vars = c("pred","rate"))
  #browser()
  p1=ggplot(data=the_data,aes(x=midage,y=value))+geom_point(aes(color=variable))+geom_line(aes(color=variable,group=variable))+
    facet_grid(~stage)+xlab("Age")+ylab("Rate/100,000")+ggtitle(paste("Mean sojourn time=",round(sojourn_time,digits=2),", ",
                                                                      "Mean duration early pre-clinical=", 
                                                                      round(earlypreclin,digits=2),sep=""))
  
  
  
  p2=ggplot(data=subset(the_data,variable=="rate"),aes(x=midage,y=value))+geom_point(aes(color=variable))+geom_line(aes(color=variable,group=variable))+
    facet_grid(~stage)+xlab("Age")+ylab("Rate/100,000")+ggtitle("Observed Cancer Incidence Data (SEER)")
  
  return(list(p1=p1,the_data=the_data))
  
  
}

summarize_fit<-function(the_data,fit,mst,dmst,k){
  
  b=exp(fit$par[k-3])
  c=1/dmst
  a=(c+b)/(c*mst)-b
  
  
  sojourn_time=b/(c*(a+b))+1/(a+b)
  preclinearly=1/(a+b)
  out=list(a=a,b=b,sojourn_time=sojourn_time,mst=mst,latemst=dmst,preclinearly=preclinearly)
  return(out)
}

get_fit<-function(mst=mst,dmst=dmst,the_data=the_data,num_seeds=num_seeds,k=k,mean1=mean1){
  print(k)
  the_fit=tryCatch(get_max_LL(the_data=the_data,num_seeds=num_seeds,mst=mst,dmst=dmst,k=k,mean1=mean1),error = function(e) e)
  while(any(class(the_fit)=="error")){
    the_fit=tryCatch(get_max_LL(the_data=the_data,num_seeds=num_seeds,mst=mst,dmst=dmst,k=k,mean1=mean1),error = function(e) e)
  }
  #browser()
  outplot=plot_observed_expected(the_data=the_data,fit=the_fit,mst=mst,dmst=dmst,k=k)
  rate.matrix=get_rate_mat(params=the_fit$par,k=k,mst=mst,dmst=dmst)
  
  summary_out=summarize_fit(the_data=the_data,fit=the_fit,mst=mst,dmst=dmst,k=k)
  return(list(the_fit=the_fit,outplot=outplot$p1,fitted_data=outplot$the_data,rate.matrix=rate.matrix,summary_out=summary_out))
}