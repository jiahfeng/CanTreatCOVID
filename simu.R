#####################################################
# simu.R
# Haolun Shi & Jiahui Feng
# update Mar. 13, 2024
# 
# clinical trial simulation for CanTreatCOVID program
#
#####################################################

library(dplyr)
library(doMC)
registerDoMC(50)
library(doRNG)
library(data.table)

source(fun.R)

simulation <- function(numarm=3, numsim=1000) {
  numint = 3000 * numarm / 500
  alloutput = foreach(sim=1:numsim, .errorhandling = c('remove'), .combine = rbind) %dorng% {
    if(numarm == 2) {
      trial = trial.init()
    }
    if(numarm == 3) {
      trial = trial.init3()
    }
    if(numarm == 4) {
      trial = trial.init4()
    }

    set.seed(sim)

    reportsoc = c()
    for(interim in 1:numint) {
      trial = trial.simulate_data(trial, 500)  
      trial = trial.compute_superiority_binary(trial)  
      # trial = trial.compute_superiority_ttc(trial)  
      trial = trial.determine_soc(trial)  
      trial = trial.compute_futility_binary(trial)  
      # trial = trial.compute_futility_ttc(trial) 
      trial = trial.determine_futility(trial) 
      trial = trial.update_randomization_prob(trial)  
      
      # if(interim %in% c(1)) {
      #   trial = trial.add_arm(trial)
      #   trial = trial.add_arm_randomization_prob(trial)
      # }

      reportsoc = c(reportsoc, trial[['index for the soc']])
    }
    reportsoc
  }
  
  omat = NULL
  for(i in 1:numint){
    
    output = alloutput[,i]
    oprobs = c(length(which(output==1))/length(output),
               length(which(output==2))/length(output),
               length(which(output==3))/length(output),
               length(which(output==4))/length(output)
    )
    omat = rbind(omat,oprobs)
  }
  
  omat=cbind(seq(1,numint)*500,omat)
  omat 
}


# simulation setting
simu = list()

simu[['response rate']] = c(0.15,0.15,0.15,0.15)
simu[['median time to recovery']] = c(10,10,10,10)
simu[['cutoff time to recovery']] = 0.99
simu[['futility cutoff time to recovery']] = 0.025
simu[['futility minimum effect time to recovery']] = 0.5
simu[['futility minimum effect binary']] = 0.0

print = base::print
print <- function(x){}
cer = 0.15
or = 1

eventrate <- function(cer,or){
  ner = cer * or/(1-cer+cer *or)
  return(ner)
}


simu[['response rate']] = c(0.05,0.05,0.05)  
simu[['cutoff binary']] = 0.9968
simu[['futility cutoff binary']] = 1 - simu[['cutoff binary']]
simulation(3, 1000)


simu[['response rate']] = c(0.05,0.05)  
simulation(2,1000)


outputtable= list()
NUMARM=3
or=1
for(cer in c(0.05,0.1,0.15)){
  numofsup=0
  ner = eventrate(cer,or)
  eventrates = c(rep(cer,NUMARM-numofsup),rep(ner,numofsup))
  simu[['response rate']] = eventrates
  omat = simulation()
  output=list(cer = cer ,or = or, eventrates = eventrates,numofsup=numofsup,omat=omat)
  outputtable[[length(outputtable)+1]] = output
}

for(cer in c(0.05,0.1,0.15)){
  for(or in c(0.5,0.6,0.7,0.8,0.9)){
    for(numofsup in 1:(NUMARM-1)){
      ner = eventrate(cer,or)
      eventrates = c(rep(cer,NUMARM-numofsup),rep(ner,numofsup))
      simu[['response rate']] = eventrates
      omat = simulation()
      output=list(cer = cer ,or = or, eventrates = eventrates,numofsup=numofsup,omat=omat)
      outputtable[[length(outputtable)+1]] = output
      base::print(length(outputtable))
    }
  }
}



outputtable2= list()
NUMARM=2
or=1
for(cer in c(0.05,0.1,0.15)){
  numofsup=0
  ner = eventrate(cer,or)
  eventrates = c(rep(cer,NUMARM-numofsup),rep(ner,numofsup))
  simu[['response rate']] = eventrates
  omat = simulation(2)
  output=list(cer = cer ,or = or, eventrates = eventrates,numofsup=numofsup,omat=omat)
  outputtable2[[length(outputtable2)+1]] = output
}

for(cer in c(0.05,0.1,0.15)){
  for(or in c(0.5,0.6,0.7,0.8,0.9)){
    for(numofsup in 1:(NUMARM-1)){
      ner = eventrate(cer,or)
      eventrates = c(rep(cer,NUMARM-numofsup),rep(ner,numofsup))
      simu[['response rate']] = eventrates
      omat = simulation(2)
      output=list(cer = cer ,or = or, eventrates = eventrates,numofsup=numofsup,omat=omat)
      outputtable2[[length(outputtable2)+1]] = output
    }
  }
}
