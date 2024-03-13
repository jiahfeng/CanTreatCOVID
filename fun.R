#######################################################
# fun.R
# Haolun Shi & Jiahui Feng
# update 13 Mar 2024
# 
# support functions for the simulation study
# 
# Functions:
#
# sampling.binary                   # simulate binary outcomes based on a given probability
# sampling.ttc                      # simulate survival outcomes
# posterior.binary                  # simulate binary outcomes from a posterior distribution
# posterior.ttc                     # generate event rate from a gamma distribution
# poptest.binary.mediandiff         # test the difference in success probabilities between two groups
#                                     against a pre-specified threshold
# poptest.ttc.mediandiff            # test the difference in median event times between two groups
#                                     against a pre-specified threshold
# poptest                           # compare the success rates of two groups using Bayesian inference
# poptest.ttc                       # compute the posterior probability that the event rate for one group
#                                     is greater than that for another group
# trial.simulate_data               # simulate and add data to a clinical trial structure
# trialdata.simulate_data           # simulate and add trial data to an existing dataset
# trial.add_arm                     # add a new arm to a clinical trial
# trial.compute_superiority_binary  # compute the posterior probabilities of superiority for binary outcomes
#                                     in each arm of a clinical trial
# trial.compute_superiority_ttc     # compute the posterior probabilities of superiority regarding time to recovery (TTC)
#                                     for each arm of a clinical trial
# trial.determine_soc               # dynamically determine the new SoC within a clinical trial
# trial.compute_futility_binary     # compute the posterior probability of futility for binary outcomes
#                                     for each arm of a clinical trial
# trial.compute_futility_ttc        # compute the posterior probabilities of futility regarding TTC
#                                     for each arm of a clinical trial
# trial.determine_futility          # evaluate and update the futility status for non-SoC arms in a clinical trial
#                                     based on binary outcomes
# trial.update_randomization_prob   # update the randomization probabilities for each arm in a clinical trial,
#                                     based on the current futility status of each arm
# trial.print                       # print a summary of a clinical trial
# trial.init                        # initialize a clinical trial with 2 arms
# trial.init3                       # initialize a clinical trial with 3 arms
# trial.init4                       # initialize a clinical trial with 4 arms
########################################################


sampling.binary <- function(p, n) {
  # p: probability
  # n: sample size
  x = rbinom(1, n, p)
  return(list(n = n, x = x))
}


sampling.ttc <- function(time.median, n) {
  # time.median: median time to event
  # n: sample size
  lambda = log(2) / time.median
  time.t = rexp(n, lambda)
  time.c = 14 # censoring
  delta = ifelse(time.t < time.c, 1, 0)
  time.y = ifelse(time.t < time.c, time.t, time.c)
  return(list(delta = delta, time = time.y))
}


posterior.binary <- function(numofsample, x, n) {
  # numofsample: # of samples to generate from the posterior distribution
  # x: # of successes
  # n: # of trials
  
  # Define prior success (ypr1) and failure (ypr2) rates
  ypr1 = 0.5
  ypr2 = 0.5
  
  y2 = x
  n2 = n
  
  return(rbeta(numofsample, ypr1 + y2, ypr2 + n2 - y2))
}


posterior.ttc <- function(numofsample, time, delta) {
  # numofsample: # of samples
  # time: time-to-event data
  # delta: event indicator
  
  a = 0.001
  b = 0.001
  
  alpha1 = a + sum(delta)
  bet1 = b + sum(time)
  return(rgamma(numofsample, shape = alpha1, rate = bet1))
}


poptest.binary.mediandiff <- function(x1, n1, x2, n2, mediandiff) {
  # x1: # of successes in group 1
  # n1: # of trials in group 1
  # x2: # of successes in group 2
  # n2: # of trials in group 2
  # mediandiff: median difference
  
  p.pos1 = posterior.binary(100000, x1, n1)
  p.pos2 = posterior.binary(100000, x2, n2)
  
  return(length(which(p.pos1 - p.pos2 > mediandiff)) / 100000)
}


poptest.ttc.mediandiff <- function(time1, delta1, time2, delta2, mediandiff) {
  # time1: event time for group 1
  # delta1: event indicator for group 1
  # time2: event time for group 2
  # delta2: event indicator for group 2
  # mediandiff: median difference
  

  lambda.pos1 = posterior.ttc(100000, time1, delta1)
  lambda.pos2 = posterior.ttc(100000, time2, delta2)

  med1 = log(2) / lambda.pos1
  med2 = log(2) / lambda.pos2

  return(length(which(med1 - med2 > mediandiff)) / 100000)
}



poptest <- function(x1, n1, y2, n2) {
  # x1: # of successes in group 1
  # n1: # of trials in group 1
  # x2: # of successes in group 2
  # n2: # of trials in group 2
  
  # Initialize prior success and failure rates for both groups to 0.5
  xpr1 = 0.5
  xpr2 = 0.5
  ypr1 = 0.5
  ypr2 = 0.5
  # Calculate the integral of the product of the complement of the Beta cumulative distribution function (CDF) 
  a = try(integrate(function(x) {
    (1 - pbeta(x, xpr1 + x1, xpr2 + n1 - x1)) * dbeta(x, ypr1 + y2, ypr2 + n2 - y2)
  }, lower = 0, upper = 1, rel.tol = 1e-4)$value)
  
  if (class(a) == 'try-error') {
    p.pos1 = posterior.binary(100000, x1, n1)
    p.pos2 = posterior.binary(100000, y2, n2)
    return(length(which(p.pos1 > p.pos2)) / 100000)
  }

  return(a)
}



poptest.ttc <- function(time1, delta1, time2, delta2) {
  # time1: event time for group 1
  # delta1: event indicator for group 1
  # time2: event time for group 2
  # delta2: event indicator for group 2

  a = 0.001
  b = 0.001
  alpha1 = a + sum(delta1)
  bet1 = b + sum(time1)
  
  alpha2 = a + sum(delta2)
  bet2 = b + sum(time2)

  a = try(integrate(function(x) {
    (1 - pgamma(x, shape = alpha1, rate = bet1)) * dgamma(x, shape = alpha2, rate = bet2)
  }, lower = 0, upper = Inf, rel.tol = 1e-4)$value)
  

  if(class(a) == "try-error") {
    lambda.pos1 = posterior.ttc(100000, time1, delta1)
    lambda.pos2 = posterior.ttc(100000, time2, delta2)
    
    return(length(which(lambda.pos1 > lambda.pos2)) / 100000)
  }
  
  return(a)
}


trial.simulate_data <- function(trial, samplesize) {
  # trial: a list representing a clinical trial. It must contain:
  #          - 'number of arm': The total number of arms in the trial.
  #          - 'randomization probability': A vector with the probability of randomization to each arm.
  #          - 'median time to recovery': A vector with the median time to recovery for each arm.
  #          - 'response rate': A vector with the response rate for each arm.
  #          - 'trialdata': A list of datasets for each arm, where each dataset has 'n', 'x', 'delta', and 'time'.
  # samplesize: sample size
  
  num = trial[['number of arm']]
  prob = trial[['randomization probability']]
  
  perarmsize = as.integer(samplesize * prob)
  
  for(i in 1:num) {
    mediantime = simu[['median time to recovery']][i]
    responserate= simu[['response rate']][i]
    trial$trialdata[[i]] = trialdata.simulate_data(trial$trialdata[[i]], perarmsize[i], mediantime, responserate)
  }
  
  return(trial)
}






trialdata.simulate_data <- function(data, samplesize, mediantime, responserate) {
  # data: input data
  # samplesize: # of subjects to simulate 
  # mediantime: the median time to event for the simulated subjects
  # responserate: the response rate (probability of success) for the simulated subjects
  
  if (samplesize == 0) {
    return(data)
  }
  
  res = sampling.binary(responserate, samplesize)
  
  data$n = data$n + res$n
  data$x = data$x + res$x
  
  ttc = sampling.ttc(mediantime, samplesize)
  data$delta = c(data$delta, ttc$delta)
  data$time = c(data$time, ttc$time)

  return(data)
}



trial.add_arm <- function(trial) {
  # trial: a list of representing a clinical trial
  
  print("add_arm")
  
  trial[['number of arm']] = trial[['number of arm']] + 1
  
  trial[['posterior probability of superiority time to recovery']] = 
    c(trial[['posterior probability of superiority time to recovery']], 0)
  trial[['posterior probability of superiority binary']] = 
    c(trial[['posterior probability of superiority binary']], 0)
  

  trial[['posterior probability of futility binary']] = 
    c(trial[['posterior probability of futility binary']], 1)
  trial[['posterior probability of futility time to recovery']] = 
    c(trial[['posterior probability of futility time to recovery']], 1)
  

  trial[['status of non-binding futility for all non-soc arm']] = 
    c(trial[['status of non-binding futility for all non-soc arm']], 0)
  trial[['status of prior-soc-become-a-loser']] = 
    c(trial[['status of prior-soc-become-a-loser']], 0)
  
  trial$trialdata[[(length(trial$trialdata) + 1)]] = list(n = 0, x = 0, time = NULL, delta = NULL)
  
  return(trial)
}


trial.compute_superiority_binary <- function(trial) {
  socindex = trial[['index for the soc']]
  numarm = trial[['number of arm']]

  indexes = seq(1, numarm)
  for (i in indexes[which(indexes != socindex)]) {
    trial[['posterior probability of superiority binary']][i] = 
      1 - poptest(trial$trialdata[[i]]$x, trial$trialdata[[i]]$n, 
                  trial$trialdata[[socindex]]$x, trial$trialdata[[socindex]]$n)
  }
  
  trial[['posterior probability of superiority binary']][socindex] = 0.5
  print('compute_superiority_binary')
  print(trial[['posterior probability of superiority binary']])
  
  return(trial)
}



trial.compute_superiority_ttc <- function(trial) {
  socindex = trial[['index for the soc']]
  numarm = trial[['number of arm']]

  indexes = seq(1, numarm)
  for (i in indexes[which(indexes != socindex)]) {
    trial[['posterior probability of superiority time to recovery']][i] = 
      poptest.ttc(trial$trialdata[[i]]$time, trial$trialdata[[i]]$delta, 
                  trial$trialdata[[socindex]]$time, trial$trialdata[[socindex]]$delta)
  }
  
  trial[['posterior probability of superiority time to recovery']][socindex] = 0
  print('compute_superiority_ttc')
  print(trial[['posterior probability of superiority time to recovery']])

  return(trial)
}


trial.determine_soc <- function(trial) {
  socindex = trial[['index for the soc']]
  print("determine soc (old soc)")
  print(socindex)
  
  prob.ttc = trial[['posterior probability of superiority time to recovery']]
  prob.binary = trial[['posterior probability of superiority binary']]
  numarm = trial[['number of arm']]
  
  exclindex = which(trial[['status of prior-soc-become-a-loser']] == 1)
  if (length(exclindex) == 0) {
    exclindex = c(9999999) # Use an unlikely high number as a placeholder to avoid empty indexing.
  }
  
  indexes = seq(1, numarm)
  
  candidate = NULL
  candidateprob = NULL
  candidatetype = NULL

  for (i in indexes[which(indexes != socindex & !(indexes %in% exclindex))]) {
    if (prob.binary[i] <= simu[['cutoff binary']]) {
      candidate = c(candidate, i)
      candidateprob = c(candidateprob, prob.binary[i])
      candidatetype = c(candidatetype, 0)
    }
    
    if (prob.binary[i] > simu[['cutoff binary']]) {
      candidate = c(candidate, i)
      candidateprob = c(candidateprob, prob.binary[i])
      candidatetype = c(candidatetype, 1)
    }
  }
  
  
  if (length(candidatetype) > 0) {
    inclindex = which(candidatetype == 1)
    if (length(inclindex) > 0) {
      candidate = candidate[inclindex]
      candidatetype = candidatetype[inclindex]
      candidateprob = candidateprob[inclindex]
      
      inclindex = which(candidateprob == max(candidateprob))
      
      candidate = candidate[inclindex]
      candidatetype = candidatetype[inclindex]
      candidateprob = candidateprob[inclindex]
      
      newsocindex = candidate[1]
      trial[['index for the soc']] = newsocindex
      trial[['status of prior-soc-become-a-loser']][socindex] = 1
      trial[['posterior probability of superiority time to recovery']][newsocindex] = 0
      print('change soc')
      print(newsocindex)
    }
  }
  
  print("determine soc (new soc)")
  print(trial[['index for the soc']])
  
  return(trial)
}



trial.compute_futility_binary <- function(trial) {
  socindex = trial[['index for the soc']]
  numarm = trial[['number of arm']]

  indexes = seq(1, numarm)
  for (i in indexes[which(indexes != socindex)]) {
    trial[['posterior probability of futility binary']][i] =
      poptest(trial$trialdata[[i]]$x, trial$trialdata[[i]]$n, 
              trial$trialdata[[socindex]]$x, trial$trialdata[[socindex]]$n)
  }
  
  trial[['posterior probability of futility binary']][socindex] = 1

  print('compute_futility_binary')
  print(trial[['posterior probability of futility binary']])
  
  return(trial)
}


trial.compute_futility_ttc <- function(trial) {

  socindex = trial[['index for the soc']]
  numarm = trial[['number of arm']]
 
  indexes = seq(1, numarm)
  for (i in indexes[which(indexes != socindex)]) {
    trial[['posterior probability of futility time to recovery']][i] = 
      poptest.ttc.mediandiff(trial$trialdata[[socindex]]$time, trial$trialdata[[socindex]]$delta, 
                             trial$trialdata[[i]]$time, trial$trialdata[[i]]$delta, 
                             simu[['futility minimum effect time to recovery']])
  }
  
  trial[['posterior probability of futility time to recovery']][socindex] = 1

  print('compute_futility_ttc')
  print(trial[['posterior probability of futility time to recovery']])
  
  return(trial)
}


trial.determine_futility <- function(trial) {
  socindex = trial[['index for the soc']]
  prob.binary = trial[['posterior probability of futility binary']]
  numarm = trial[['number of arm']]
  
  oldstatus = trial[['status of non-binding futility for all non-soc arm']]

  oldstatus[socindex] = 0
  indexes = seq(1, numarm)
  
  for (i in indexes[which(indexes != socindex)]) {
    if (prob.binary[i] < simu[['futility cutoff binary']]) {
      oldstatus[i] = 1
    }
  }

  trial[['status of non-binding futility for all non-soc arm']] = oldstatus
  
  print('determine_futility')
  print(trial[['status of non-binding futility for all non-soc arm']])

  return(trial)
}


trial.update_randomization_prob <- function(trial) {
  socindex = trial[['index for the soc']]
  numarm = trial[['number of arm']]
  futilitystatus = trial[['status of non-binding futility for all non-soc arm']]

  eligibleindex = which(futilitystatus == 0)
  numarm_e = length(eligibleindex)
  
  trial[['randomization probability']] = rep(1/numarm, numarm)
  trial[['randomization probability']][socindex] = 1/numarm_e

  indexes = seq(1, numarm)
  for (i in indexes) {
    if (futilitystatus[i] == 1) {
      trial[['randomization probability']][i] = 0
    } else if (futilitystatus[i] == 0) {
      trial[['randomization probability']][i] = 1/numarm_e
    }
  }
  
  trial[['randomization probability']] = trial[['randomization probability']] / sum(trial[['randomization probability']])
  
  print('update_randomization_prob')
  print(trial[['randomization probability']])

  return(trial)
}


trial.print <- function(trial) {
  numarm = trial[['number of arm']]

  tbl1 = data.frame(Successes = integer(),
                    Trials = integer(),
                    Success_Rate = numeric(),
                    Median_Time = numeric(),
                    stringsAsFactors = FALSE)

  for (i in 1:numarm) {
    if (trial$trialdata[[i]]$n > 0) {
      successes = trial$trialdata[[i]]$x
      trials = trial$trialdata[[i]]$n
      success_rate = successes / trials
      median_time = ifelse(is.null(trial$trialdata[[i]]$time), NA, median(trial$trialdata[[i]]$time, na.rm = TRUE))
      tbl1 = rbind(tbl1, c(Successes = successes, Trials = trials, Success_Rate = success_rate, Median_Time = median_time))
    }
  }
  
  print('Trial Summary')
  print(tbl1)
  print('------------------------------------------------')
}


trial.init <- function(){
  trial=list()
  trial[['number of arm']] = 2
  trial[['randomization probability']] = c(0.5,0.5)
  trial[['index for the soc']] = 1
  trial[['index for uc']] = 1
  trial[['posterior probability of superiority time to recovery']] = c(0,0)
  trial[['posterior probability of superiority binary']] = c(0,0)
  trial[['posterior probability of futility time to recovery']] = c(1,1)
  trial[['posterior probability of futility binary']] = c(1,1)
  trial[['status of non-binding futility for all non-soc arm']] = c(0,0)
  trial[['status of prior-soc-become-a-loser']] = c(0,0)
  
  
  trial$trialdata=list()
  trial$trialdata[[1]]=list(n=0,x=0,time=NULL,delta=NULL)
  trial$trialdata[[2]]=list(n=0,x=0,time=NULL,delta=NULL)
  
  return(trial)
}



trial.init3 <- function(){
  trial=list()
  trial[['number of arm']] = 3
  trial[['randomization probability']] = rep(1/3,3)
  trial[['index for the soc']] = 1
  trial[['index for uc']] = 1
  trial[['posterior probability of superiority time to recovery']] = c(0,0,0)
  trial[['posterior probability of superiority binary']] = c(0,0,0)
  trial[['posterior probability of futility time to recovery']] = c(1,1,1)
  trial[['posterior probability of futility binary']] = c(1,1,1)
  trial[['status of non-binding futility for all non-soc arm']] = c(0,0,0)
  trial[['status of prior-soc-become-a-loser']] = c(0,0,0)
  
  
  trial$trialdata=list()
  trial$trialdata[[1]]=list(n=0,x=0,time=NULL,delta=NULL)
  trial$trialdata[[2]]=list(n=0,x=0,time=NULL,delta=NULL)
  trial$trialdata[[3]]=list(n=0,x=0,time=NULL,delta=NULL)
  return(trial)
}


trial.init4 <- function(){
  trial=list()
  trial[['number of arm']] = 4
  trial[['randomization probability']] = c(0.25,0.25,0.25,0.25)
  trial[['index for the soc']] = 1
  trial[['index for uc']] = 1
  trial[['posterior probability of superiority time to recovery']] = c(0,0,0,0)
  trial[['posterior probability of superiority binary']] = c(0,0,0,0)
  trial[['posterior probability of futility time to recovery']] = c(1,1,1,1)
  trial[['posterior probability of futility binary']] = c(1,1,1,1)
  trial[['status of non-binding futility for all non-soc arm']] = c(0,0,0,0)
  trial[['status of prior-soc-become-a-loser']] = c(0,0,0,0)
  
  
  trial$trialdata=list()
  trial$trialdata[[1]]=list(n=0,x=0,time=NULL,delta=NULL)
  trial$trialdata[[2]]=list(n=0,x=0,time=NULL,delta=NULL)
  trial$trialdata[[3]]=list(n=0,x=0,time=NULL,delta=NULL)
  trial$trialdata[[4]]=list(n=0,x=0,time=NULL,delta=NULL)
  return(trial)
  
}




