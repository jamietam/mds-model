rm(list = ls())
mainDir <- "C:/Users/JT936/Dropbox/GitHub/mds-model"

setwd(file.path(mainDir))
library(openxlsx)
library(reshape)
library(splines)
library(Bhat)
library(data.table)
library(ggplot2)
library(gridExtra)
load("mdseprevs0520.rda")

namethisrun = "_090523_tx" 
folder = paste0("Koop",namethisrun,"/") # name the folder where results will be saved

# Read inputs -------------------------------------------------------------
startyear = 1900 # the burn-in period starting point
endyear = 2100 # 2014 or 2050 , still need census projections through 2065
startage = 0
endage = 99
Ny= endyear - startyear + 1 
Na= endage - startage + 1
emptycompartment <- matrix(0, nrow = Na, ncol = Ny, dimnames=list(c(startage:endage),c(startyear:endyear))) # create matrix of zeroes for compartments
policystart = 2024

allparamsF = read.xlsx("parameters_tx.xlsx",sheet=paste0("model_females"),rowNames=TRUE,colNames=TRUE) # Adjust parameters in this excel file with estimates to be used in the model
paramsF = as.vector(subset(allparamsF,bhat==1)[['estimate']]) # Parameters where bhat = 1 can be estimated by bhat
paramsnamesF = rownames(subset(allparamsF,bhat==1))  
lowervectorF = as.vector(subset(allparamsF,bhat==1)[['lower']]) # lower bounds for bhat parameters
uppervectorF = as.vector(subset(allparamsF,bhat==1)[['upper']]) # upper bounds for bhat parameters

allparamsM = read.xlsx("parameters_tx.xlsx",sheet=paste0("model_males"),rowNames=TRUE,colNames=TRUE) # Adjust parameters in this excel file with estimates to be used in the model
paramsM = as.vector(subset(allparamsM,bhat==1)[['estimate']]) # Parameters where bhat = 1 can be estimated by bhat
paramsnamesM = rownames(subset(allparamsM,bhat==1))  
lowervectorM = as.vector(subset(allparamsM,bhat==1)[['lower']]) # lower bounds for bhat parameters
uppervectorM = as.vector(subset(allparamsM,bhat==1)[['upper']]) # upper bounds for bhat parameters

agerownames<-c("18to25", "26to34", "35to49", "50to64",  "65plus", "total")
agegroupstart <- c(18,26,35,50,65,18)
agegroupend <- c(25,34,49,64,99,99)

# Get model prevs by age group --------------------------------------------
getmodelprevs <- function(numerator,denominator){
  numerator = as.data.frame(numerator)
  denominator = as.data.frame(denominator)
  prevs = NULL
  for (a in c(1:length(agegroupstart))) {
    prevs <- rbind(prevs,colSums(numerator[c((agegroupstart[a]+1):(agegroupend[a]+1)), ],na.rm=TRUE)/
                     colSums(denominator[c((agegroupstart[a]+1):(agegroupend[a]+1)), ],na.rm=TRUE))
  }
  row.names(prevs)<-agerownames 
  return(prevs)
}

# Main model --------------------------------------------------------------
main <- function(getmodelprevs, whichgender, allparamsF, paramsF,paramsnamesF, mpc, txeffcess, util){
  setwd(file.path(mainDir))
  
  pop = read.xlsx("census_data/np2017_d1.xlsx",sheet=whichgender,rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  
  # smk params --------------------------------------------------------------
  
  smk_init_cisnet = read.xlsx("cisnet_smkrates_by_calyear_061923.xlsx",sheet=paste0(whichgender,"_init"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  smk_cess_cisnet = read.xlsx("cisnet_smkrates_by_calyear_061923.xlsx",sheet=paste0(whichgender,"_cess"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  # smk_init_cisnet = read.xlsx("cisnet_smkrates_nhis2018.xlsx",sheet=paste0(whichgender,"_init"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  # smk_cess_cisnet = read.xlsx("cisnet_smkrates_nhis2018.xlsx",sheet=paste0(whichgender,"_cess"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  
  smkinit_SF = matrix(0,100,1)
  smkinit_SF[1:18] = ifelse(allparamsF["smkinit_youthSF","bhat"]==0,allparamsF["smkinit_youthSF","estimate"],paramsF[match("smkinit_youthSF",paramsnames)])
  smkinit_SF[19:35] = ifelse(allparamsF["smkinit_SF_18to34","bhat"]==0,allparamsF["smkinit_SF_18to34","estimate"],paramsF[match("smkinit_SF_18to34",paramsnamesF)])
  smkinit_SF[36:66] = ifelse(allparamsF["smkinit_SF_35to64","bhat"]==0,allparamsF["smkinit_SF_35to64","estimate"],paramsF[match("smkinit_SF_35to64",paramsnamesF)])
  smkinit_SF[67:100] = ifelse(allparamsF["smkinit_SF_65plus","bhat"]==0,allparamsF["smkinit_SF_65plus","estimate"],paramsF[match("smkinit_SF_65plus",paramsnamesF)])
  
  smkcess_SF = matrix(0,100,1)
  smkcess_SF[1:18] = ifelse(allparamsF["smkcess_youthSF","bhat"]==0,allparamsF["smkcess_youthSF","estimate"],paramsF[match("smkcess_youthSF",paramsnamesF)])
  smkcess_SF[19:35] = ifelse(allparamsF["smkcess_SF_18to34","bhat"]==0,allparamsF["smkcess_SF_18to34","estimate"],paramsF[match("smkcess_SF_18to34",paramsnamesF)])
  smkcess_SF[36:66] = ifelse(allparamsF["smkcess_SF_35to64","bhat"]==0,allparamsF["smkcess_SF_35to64","estimate"],paramsF[match("smkcess_SF_35to64",paramsnamesF)])
  smkcess_SF[67:100] = ifelse(allparamsF["smkcess_SF_65plus","bhat"]==0,allparamsF["smkcess_SF_65plus","estimate"],paramsF[match("smkcess_SF_65plus",paramsnamesF)])
  
  smk_init = smk_init_cisnet*smkinit_SF # scale smoking initiation rates ### CHECK DIMENSIONS might be off by one
  smk_cess = smk_cess_cisnet*smkcess_SF # scale smoking cessation rates
  smk_cess_txpop = smk_cess # when mpc ==0

  death_ns = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("ns_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  death_cs = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("cs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  death_fs = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("fs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 

  LE_ns = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_ns_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  LE_cs = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_cs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  LE_fs = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_fs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 

  ucs_minus_uns = death_cs - death_ns
  ufs_minus_uns = death_fs - death_ns
  
  dep1inc = read.xlsx("incidence_eaton.xlsx",sheet=whichgender,rowNames=TRUE, colNames=FALSE, check.names=FALSE)
  
  # model parameters are retrieved from the parameters.xlsx file, if they are bhat=1, then they are read in as a params vector for bhat estimation
  RRdepr_death =  c(rep(1,17),rep(ifelse(allparamsF["RRdepr_death","bhat"]==0,allparamsF["RRdepr_death","estimate"],paramsF[match("RRdepr_death",paramsnamesF)]),Na-18))
  
  # scale all recovery rates by the scaling factor deprecov_SF
  deprecov_rate =  c(rep(ifelse(allparamsF["deprecov_rate","bhat"]==0,allparamsF["deprecov_rate","estimate"],paramsF[match("deprecov_rate",paramsnamesF)]),Na-1))
  deprecov_SF = matrix(1,99,1)
  deprecov_rate = deprecov_rate *  deprecov_SF
  deprecovSF_cs =  c(rep(ifelse(allparamsF["deprecovSF_cs","bhat"]==0,allparamsF["deprecovSF_cs","estimate"],paramsF[match("deprecovSF_cs",paramsnamesF)]),Na-1))
  deprecovSF_fs =  c(rep(ifelse(allparamsF["deprecovSF_fs","bhat"]==0,allparamsF["deprecovSF_fs","estimate"],paramsF[match("deprecovSF_fs",paramsnamesF)]),Na-1))
  
  deprinc = rep(ifelse(allparamsF["depr_inc","bhat"]==0,allparamsF["depr_inc","estimate"],paramsF[match("depr_inc",paramsnamesF)]),Na-1)
  deprinc_SF = matrix(1,99,1)
  depr_inc = deprinc * deprinc_SF
  
  depinc1 = ifelse(allparamsF["depinc1","bhat"]==0,allparamsF["depinc1","estimate"],paramsF[match("depinc1",paramsnamesF)])
  depinc2 = ifelse(allparamsF["depinc2","bhat"]==0,allparamsF["depinc2","estimate"],paramsF[match("depinc2",paramsnamesF)])
  depinc3 = ifelse(allparamsF["depinc3","bhat"]==0,allparamsF["depinc3","estimate"],paramsF[match("depinc3",paramsnamesF)])
  
  inc_SF = ifelse(allparamsF["inc_SF","bhat"]==0,allparamsF["inc_SF","estimate"],paramsF[match("inc_SF",paramsnamesF)])
  
  if (whichgender=="females"){
    MP=ns(0:21,knots=c(13,18)) ### Matrix of X's
    Rps=predict(MP,21)[1,] ## Predicts y value given a set of X's for age 22  # 0.0051285304 anchor at age 22
    y=c()
    for (j in 0:21){
      y=c(y,0.0051285304*exp(sum((MP[j,]-Rps)*c(depinc1,depinc2,depinc3)))) ## multiply by coefficients and sum , exp makes it positive for incidence
    }
    dep1_inc=dep1inc
    dep1_inc[0:22,]<-y
    dep1_inc[0:12,]<-rep(0,12) # assumes no 1st MDE before age 12
    scaleddep1_inc = dep1_inc
    scaleddep1_inc[0:26,]<-scaleddep1_inc[0:26,]*inc_SF 
  }
    
  if (whichgender=="males"){
    MP=ns(0:28,knots=c(13,18))
    Rps=predict(MP,28)[1,] ## Predicts y value given a set of X's for age 22  # 0.0019072600 anchor at age 29
    y=c()
    for (j in 0:28){
      y=c(y,0.0019072600*exp(sum((MP[j,]-Rps)*c(depinc1,depinc2,depinc3)))) ## multiply by coefficients and sum , exp makes it positive for incidence
    }
    dep1_inc=dep1inc
    dep1_inc[0:29,]<-y
    dep1_inc[0:12,]<-rep(0,12) # assumes no 1st MDE before age 12
    scaleddep1_inc = dep1_inc
    scaleddep1_inc[0:26,]<-scaleddep1_inc[0:26,]*inc_SF # multiply all incidence probabilities by scaling factor for ages <=25
  }
  
  # Age-group categorical forgetting probabilities  -------------------------
  forget = matrix(NA,99,1)
  forget[1:17] = rep(0,17) # assumes no forgetting before age 18
  forget[18:25] = ifelse(allparamsF["forget1","bhat"]==0,allparamsF["forget1","estimate"],paramsF[match("forget1",paramsnamesF)])
  forget[26:34] = ifelse(allparamsF["forget2","bhat"]==0,allparamsF["forget2","estimate"],paramsF[match("forget2",paramsnamesF)])
  forget[35:49] = ifelse(allparamsF["forget3","bhat"]==0,allparamsF["forget3","estimate"],paramsF[match("forget3",paramsnamesF)])
  forget[50:64] = ifelse(allparamsF["forget4","bhat"]==0,allparamsF["forget4","estimate"],paramsF[match("forget4",paramsnamesF)])
  forget[65:99] = ifelse(allparamsF["forget5","bhat"]==0,allparamsF["forget5","estimate"],paramsF[match("forget5",paramsnamesF)])

  
  # smkdep effects parameters  ------------------------------------------
  RRcs_dep1 = c(rep(ifelse(allparamsF["RRcs_dep1","bhat"]==0,allparamsF["RRcs_dep1","estimate"],paramsF[match("RRcs_dep1",paramsnamesF)]),Na-1)) # Note: RR estimate comes from adult survey but applies to youth in model
  RRfs_dep1 = c(rep(ifelse(allparamsF["RRfs_dep1","bhat"]==0,allparamsF["RRfs_dep1","estimate"],paramsF[match("RRfs_dep1",paramsnamesF)]),Na-1))
  Efs_depr =  c(rep(ifelse(allparamsF["Efs_depr","bhat"]==0,allparamsF["Efs_depr","estimate"],paramsF[match("Efs_depr",paramsnamesF)]),Na-1))  
  Ecs_depr =  c(rep(ifelse(allparamsF["Ecs_depr","bhat"]==0,allparamsF["Ecs_depr","estimate"],paramsF[match("Ecs_depr",paramsnamesF)]),Na-1))  
  ORhdep_quit =  c(rep(ifelse(allparamsF["ORhdep_quit","bhat"]==0,allparamsF["ORhdep_quit","estimate"],paramsF[match("ORhdep_quit",paramsnamesF)]),Na-1))
  Edepr_smkinit =  c(rep(ifelse(allparamsF["Edepr_smkinit","bhat"]==0,allparamsF["Edepr_smkinit","estimate"],paramsF[match("Edepr_smkinit",paramsnamesF)]),Na-1))
  
  ORhdep_quit_txpop = ORhdep_quit # when mpc = 0
  
  # mental health service utilization parameters ----------------------------
  if ('ahltmde_18to25' %in% paramsnamesF){ # if doing bhat estimation for utilization rates
    proptx = matrix(NA,99,1)
    # proptx[1:17] = rep(0,17) # assumes no one is in tx before age 18
    proptx[1:25] = ifelse(allparamsF["ahltmde_18to25","bhat"]==0,allparamsF["ahltmde_18to25","estimate"],paramsF[match("ahltmde_18to25",paramsnamesF)])
    proptx[26:34] = ifelse(allparamsF["ahltmde_26to34","bhat"]==0,allparamsF["ahltmde_26to34","estimate"],paramsF[match("ahltmde_26to34",paramsnamesF)])
    proptx[35:49] = ifelse(allparamsF["ahltmde_35to49","bhat"]==0,allparamsF["ahltmde_35to49","estimate"],paramsF[match("ahltmde_35to49",paramsnamesF)])
    proptx[50:64] = ifelse(allparamsF["ahltmde_50to64","bhat"]==0,allparamsF["ahltmde_50to64","estimate"],paramsF[match("ahltmde_50to64",paramsnamesF)])
    proptx[65:99] = ifelse(allparamsF["ahltmde_65plus","bhat"]==0,allparamsF["ahltmde_65plus","estimate"],paramsF[match("ahltmde_65plus",paramsnamesF)])
  }
  if (!('ahltmde_18to25'%in% paramsnamesF)){
    utilrates = read.xlsx("parameters_tx.xlsx",sheet=paste0("util_", whichgender),rowNames=TRUE,colNames=TRUE)
    colnames(utilrates) <- c("1.0","1.1","1.2","100")
    proptx = matrix(NA,99,1)
    # proptx[1:17] = utilrates["ahltmde_18to25",util] #rep(0,17) # assumes same proportion is in tx before age 18
    proptx[1:25] = utilrates["ahltmde_18to25",util]
    proptx[26:34] = utilrates["ahltmde_26to34",util]
    proptx[35:49] = utilrates["ahltmde_35to49",util]
    proptx[50:64] = utilrates["ahltmde_50to64",util]
    proptx[65:99] = utilrates["ahltmde_65plus",util]
  }

  # Compartments / state variables ------------------------------------------
  
  # Initialize population - Model compartments are organized with age 0-99 as rows, and year as columns
  matrix.names<-c('ns_nevdep', 'ns_dep1','ns_fdep', 'ns_dep2','ns_recall',
                  'cs_nevdep', 'cs_dep1_tx','cs_dep1_notx','cs_fdep', 'cs_dep2_tx','cs_dep2_notx','cs_recall',
                  'fs_nevdep', 'fs_dep1','fs_fdep', 'fs_dep2','fs_recall','ns_deaths', 'cs_deaths', 'fs_deaths')
  for (name in matrix.names) assign(name,emptycompartment)

  ns_nevdep[paste(startage),1:Ny] <- as.matrix(pop[paste(startage),1:Ny])   # Takes empty compartment and populates the top row of the matrix with the number of 0-yrolds

  # Run it ------------------------------------------------------------------

  for (y in c((startyear+1):(endyear))){
    py = paste(y - 1)
    if (y<policystart){ 
      txeff = c(rep(1, Na-1))
    } else {
      txeff = c(rep(1,17),rep(txeffcess,Na-1-17)) # Treatment effect applies to those ages 18+,
    }
	  if (mpc==1 & y>=policystart){ # when mpc == 1. Under mpc scenario, 100% quitting among treated smkdep pop starting in policy year
      ORhdep_quit_txpop = 1 # so that all smokers in treatment quit immediately (1-1=0)
      smk_cess_txpop[19:100,120:201]=1 # Treatment effect applies to those ages 18+, 2019-2100 - need to include 2019 because of 'py'
    }
	  if (y>=2016){
      usethis = scaleddep1_inc # from 2016-2100, use scaled incidence probabilities
    } else {
      usethis = dep1_inc # all other years, use the original incidence probabilities						
    }
    
    # Never MDE
    ns_nevdep[2:Na,paste(y)] <- ns_nevdep[1:Na-1,py]*(1-smk_init[(startage+1):endage,py])*(1-usethis[(startage+1):(endage),])*(1-death_ns[(startage+1):endage,py])
    cs_nevdep[2:Na,paste(y)] <- cs_nevdep[1:Na-1,py]*(1-smk_cess[(startage+1):endage,py])*(1-RRcs_dep1*usethis[(startage+1):(endage),])*(1-death_cs[(startage+1):endage,py]) + ns_nevdep[1:Na-1,py]*(smk_init[(startage+1):endage,py])
    fs_nevdep[2:Na,paste(y)] <- fs_nevdep[1:Na-1,py]*(1-RRfs_dep1*usethis[(startage+1):(endage),])*(1-death_fs[(startage+1):endage,py]) + cs_nevdep[1:Na-1,py]*(smk_cess[(startage+1):endage,py])
    
    # Current MDE (1st episode)
    
    ns_dep1[2:Na,paste(y)] <- ns_dep1[1:Na-1,py]*(1-Edepr_smkinit*smk_init[(startage+1):endage,py])*(1-deprecov_rate)*(1-RRdepr_death*death_ns[(startage+1):endage,py]) + ns_nevdep[1:Na-1,py]*(usethis[(startage+1):(endage),])
    
    cs_dep1_tx[2:Na,paste(y)] <-  cs_dep1_tx[1:Na-1,py]*(1-ORhdep_quit_txpop*txeff*smk_cess_txpop[(startage+1):endage,py])*(1-deprecovSF_cs*deprecov_rate)*(1-RRdepr_death*death_cs[(startage+1):endage,py]) +  proptx*cs_nevdep[1:Na-1,py]*(RRcs_dep1*usethis[(startage+1):(endage),]) + proptx*ns_dep1[1:Na-1,py]*(Edepr_smkinit*smk_init[(startage+1):endage,py])
    cs_dep1_notx[2:Na,paste(y)] <- cs_dep1_notx[1:Na-1,py]*(1-ORhdep_quit*smk_cess[(startage+1):endage,py])*(1-deprecovSF_cs*deprecov_rate)*(1-RRdepr_death*death_cs[(startage+1):endage,py]) +  (1-proptx)*cs_nevdep[1:Na-1,py]*(RRcs_dep1*usethis[(startage+1):(endage),]) + (1-proptx)*ns_dep1[1:Na-1,py]*(Edepr_smkinit*smk_init[(startage+1):endage,py])
    
    cs_dep1 = cs_dep1_tx +  cs_dep1_notx
    
    fs_dep1[2:Na,paste(y)] <- fs_dep1[1:Na-1,py]*(1-deprecovSF_fs*deprecov_rate)*(1-RRdepr_death*death_fs[(startage+1):endage,py]) + fs_nevdep[1:Na-1,py]*(RRfs_dep1*usethis[(startage+1):(endage),]) + cs_dep1_notx[1:Na-1,py]*(ORhdep_quit*smk_cess[(startage+1):endage,py]) + cs_dep1_tx[1:Na-1,py]*(ORhdep_quit_txpop*txeff*smk_cess_txpop[(startage+1):endage,py])
    
    
    # Current MDE (Recurrent episode)
    ns_dep2[2:Na,paste(y)] <- ns_dep2[1:Na-1,py]*(1-Edepr_smkinit*smk_init[(startage+1):endage,py])*(1-deprecov_rate)*(1-RRdepr_death*death_ns[(startage+1):endage,py]) + ns_fdep[1:Na-1,py]*(depr_inc[(startage+1):(endage),]) + ns_recall[1:Na-1,py]*(depr_inc[(startage+1):(endage),])
    
    cs_dep2_tx[2:Na,paste(y)] <- cs_dep2_tx[1:Na-1,py] *(1-ORhdep_quit_txpop*txeff*smk_cess_txpop[(startage+1):endage,py])*(1-deprecovSF_cs*deprecov_rate)*(1-RRdepr_death*death_cs[(startage+1):endage,py]) + proptx*cs_fdep[1:Na-1,py]*(Ecs_depr*depr_inc[(startage+1):(endage),]) + proptx*ns_dep2[1:Na-1,py]*(Edepr_smkinit*smk_init[(startage+1):endage,py]) + proptx*cs_recall[1:Na-1,py]*(depr_inc[(startage+1):(endage),])
    
    cs_dep2_notx[2:Na,paste(y)] <- cs_dep2_notx[1:Na-1,py] *(1-ORhdep_quit*smk_cess[(startage+1):endage,py])*(1-deprecovSF_cs*deprecov_rate)*(1-RRdepr_death*death_cs[(startage+1):endage,py]) +  (1-proptx)*cs_fdep[1:Na-1,py]*(Ecs_depr*depr_inc[(startage+1):(endage),]) + (1-proptx)*ns_dep2[1:Na-1,py]*(Edepr_smkinit*smk_init[(startage+1):endage,py]) + (1-proptx)*cs_recall[1:Na-1,py]*(depr_inc[(startage+1):(endage),])
    
    cs_dep2 = cs_dep2_tx +  cs_dep2_notx
    
    fs_dep2[2:Na,paste(y)] <- fs_dep2[1:Na-1,py]*(1-deprecovSF_fs*deprecov_rate)*(1-RRdepr_death*death_fs[(startage+1):endage,py]) + fs_fdep[1:Na-1,py]*(Efs_depr*depr_inc[(startage+1):(endage),]) + cs_dep2_tx[1:Na-1,py]*(ORhdep_quit_txpop*txeff*smk_cess_txpop[(startage+1):endage,py]) +cs_dep2_notx[1:Na-1,py]*(ORhdep_quit*smk_cess[(startage+1):endage,py])+fs_recall[1:Na-1,py]*(depr_inc[(startage+1):(endage),])
    
    # Former MDE - Recovery with risk of recurrence  
    ns_fdep[2:Na,paste(y)] <- ns_fdep[1:Na-1,py]*(1-smk_init[(startage+1):endage,py])*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_ns[(startage+1):endage,py])*(1-forget[(startage+1):(endage)]) + ns_dep1[1:Na-1,py]*(deprecov_rate) + ns_dep2[1:Na-1,py]*(deprecov_rate)
    
    cs_fdep[2:Na,paste(y)] <- cs_fdep[1:Na-1,py]*(1-ORhdep_quit*smk_cess[(startage+1):endage,py])*(1-Ecs_depr*depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_cs[(startage+1):endage,py])*(1-forget[(startage+1):(endage)])  + cs_dep1_tx[1:Na-1,py]*(deprecovSF_cs*deprecov_rate) + cs_dep2_tx[1:Na-1,py]*(deprecovSF_cs*deprecov_rate) + cs_dep1_notx[1:Na-1,py]*(deprecovSF_cs*deprecov_rate) + cs_dep2_notx[1:Na-1,py]*(deprecovSF_cs*deprecov_rate) + ns_fdep[1:Na-1,py]*(smk_init[(startage+1):endage,py])
    
    fs_fdep[2:Na,paste(y)] <- fs_fdep[1:Na-1,py]*(1-Efs_depr*depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_fs[(startage+1):endage,py])*(1-forget[(startage+1):(endage)]) + fs_dep1[1:Na-1,py]*(deprecovSF_fs*deprecov_rate) + fs_dep2[1:Na-1,py]*(deprecovSF_fs*deprecov_rate) + cs_fdep[1:Na-1,py]*(ORhdep_quit*smk_cess[(startage+1):endage,py])
    
    # Recall Error 
    ns_recall[2:Na,paste(y)] <- ns_recall[1:Na-1,py]*(1-smk_init[(startage+1):endage,py])*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_ns[(startage+1):endage,py]) + ns_fdep[1:Na-1,py]*(forget[(startage+1):(endage)])
  
  	cs_recall[2:Na,paste(y)] <- cs_recall[1:Na-1,py]*(1-ORhdep_quit*smk_cess[(startage+1):endage,py])*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_cs[(startage+1):endage,py]) + ns_recall[1:Na-1,py]*(smk_init[(startage+1):endage,py]) + cs_fdep[1:Na-1,py]*(forget[(startage+1):(endage)])
    
    fs_recall[2:Na,paste(y)] <- fs_recall[1:Na-1,py]*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_fs[(startage+1):endage,py])+ cs_recall[1:Na-1,py]*(ORhdep_quit*smk_cess[(startage+1):endage,py]) + fs_fdep[1:Na-1,py]*(forget[(startage+1):(endage)])
    
  }

  # Get inputs --------------------------------------------------------------
  
  initiation <- data.frame(cbind(smk_init_cisnet[,"2018"],smk_init_cisnet[,"2018"],smkinit_SF,smk_init["2018"],c(0:99),c(Edepr_smkinit,Edepr_smkinit[99])))
  colnames(initiation) <- c("2018", "cisnet2018","SF","scaledrates","age","Edepr_smkinit")
  
  cessation <- data.frame(cbind(smk_cess_cisnet[,"2018"],smk_cess_cisnet[,"2018"],smkcess_SF,smk_cess["2018"],c(0:99),c(ORhdep_quit,Edepr_smkinit[99])))
  colnames(cessation) <- c("2018", "cisnet2018","SF","scaledrates","age","ORhdep_quit")
  
  incidence <- data.frame(cbind(dep1inc[-1,],dep1_inc[-1,],scaleddep1_inc[-1,],deprinc,deprinc_SF,data.frame(depr_inc),c(1:99)))
  colnames(incidence) <- c("dep1inc","splines","scaleddep1inc_yr", "deprinc","SF_dep2", "scaledrates_dep2", "age")
  
  recovery <- data.frame(cbind(deprecov_rate, deprecov_SF,deprecov_rate,deprecovSF_fs,deprecovSF_cs,c(1:99)))
  colnames(recovery) <- c("dep1recov","SF","scaledrates","formersmokers","currentsmokers", "age")
  
  forget <- data.frame(cbind(forget,c(1:99)))
  colnames(forget) <- c("forget_prob","age")
  
  # Get population counts for each subpopulation ----------------------------
  
  nspop = (ns_nevdep + ns_recall + ns_dep2 + ns_fdep + ns_dep1) 
  cspop = (cs_nevdep + cs_recall + cs_dep2 + cs_fdep + cs_dep1) 
  fspop = (fs_nevdep + fs_recall + fs_dep2 + fs_fdep + fs_dep1) 
  
  totalpop = nspop+cspop+fspop
  
  ns_dep <- (ns_dep2+ns_dep1) 
  cs_dep <- (cs_dep2+cs_dep1) 
  fs_dep <- (fs_dep2+fs_dep1) 
  
  nevdeppop <- (ns_nevdep + cs_nevdep + fs_nevdep)
  fdeppop <- (ns_fdep + cs_fdep + fs_fdep + ns_recall + cs_recall + fs_recall)
  deppop <- (ns_dep + cs_dep + fs_dep)
  recallpop <- ns_recall + cs_recall + fs_recall
  
  cs_everdep = cs_recall + cs_dep2 + cs_fdep + cs_dep1
  fs_everdep = fs_recall + fs_dep2 + fs_fdep + fs_dep1
  ns_everdep = ns_recall + ns_dep2 + ns_fdep + ns_dep1
  
  cs_dep_tx = cs_dep1_tx + cs_dep2_tx
  cs_dep_notx = cs_dep1_notx + cs_dep2_notx
  
  ns_ndep <- (ns_nevdep+ns_fdep+ns_recall) 
  cs_ndep <- (cs_nevdep+cs_fdep+cs_recall) 
  fs_ndep <- (fs_nevdep+fs_fdep+fs_recall) 
  
  ndeppop = ns_ndep+cs_ndep+fs_ndep
    
  # Get prevalences for model fitting ---------------------------------------
  
  # s1 <- getmodelprevs(ns_nevdep,nevdeppop) # smoker prevalence among the entire never depressed population
  # s2 <- getmodelprevs(cs_nevdep,nevdeppop)
  # s3 <- getmodelprevs(fs_nevdep,nevdeppop)
  # s4 <- getmodelprevs(ns_fdep,fdeppop) # smoker prevalence among the population w/ history of depression
  # s5 <- getmodelprevs(cs_fdep,fdeppop)
  # s6 <- getmodelprevs(fs_fdep,fdeppop)
  # s7 <- getmodelprevs(ns_dep,deppop)  # smoker prevalence among the depressed population 
  # s8 <- getmodelprevs(cs_dep,deppop)
  # s9 <- getmodelprevs(fs_dep,deppop)
  # modelsmkprevdata <- list(s1,s2,s3,s4,s5,s6,s7,s8,s9)
  # 
  # d1 <- getmodelprevs(ns_nevdep,nspop) # dep prevalence among never smokers
  # d2 <- getmodelprevs(ns_fdep,nspop)
  # d3 <- getmodelprevs(ns_dep,nspop)
  # d4 <- getmodelprevs(cs_nevdep,cspop) # dep prevalence among current smokers
  # d5 <- getmodelprevs(cs_fdep,cspop) 
  # d6 <- getmodelprevs(cs_dep,cspop)
  # d7 <- getmodelprevs(fs_nevdep,fspop) # dep prevalence among former smokers
  # d8 <- getmodelprevs(fs_fdep,fspop)
  # d9 <- getmodelprevs(fs_dep,fspop) 
  # modeldepprevdata <- list(d1,d2,d3,d4,d5,d6,d7,d8,d9)
  
  s1 <- getmodelprevs(ns_dep,deppop)  # smoker prevalence among the depressed population 
  s2 <- getmodelprevs(cs_dep,deppop)
  s3 <- getmodelprevs(fs_dep,deppop)
  s4 <- getmodelprevs(nspop,totalpop)
  s5 <- getmodelprevs(cspop,totalpop)
  s6 <- getmodelprevs(fspop,totalpop)
  s7 <- getmodelprevs(cs_ndep,ndeppop)
  modelsmkprevdata <- list(s1,s2,s3,s4,s5,s6,s7)
  
  d1 <- getmodelprevs(ns_dep,nspop) # dep prevalence among never smokers
  d2 <- getmodelprevs(cs_dep,cspop) # dep prevalence among current smokers
  d3 <- getmodelprevs(fs_dep,fspop) # dep prevalence among former smokers
  d4<- getmodelprevs(deppop,totalpop)
  modeldepprevdata <- list(d1,d2,d3,d4)
  
  return(list(modelsmkprevdata,modeldepprevdata,
              initiation,cessation,incidence,recovery,forget,
              totalpop, 
              nspop, cspop, fspop, 
              nevdeppop, deppop, fdeppop,recallpop, 
              ns_nevdep, cs_nevdep, fs_nevdep, 
              ns_everdep, cs_everdep, fs_everdep,
              ns_dep, cs_dep, fs_dep,
              ns_fdep, cs_fdep, fs_fdep, 
              ns_recall, cs_recall, fs_recall,
              cs_dep_tx, cs_dep_notx
              ))
}

# Get NSDUH prevs by age group --------------------------------------------
# Generates the corresponding NSDUH prevalences so that we can compare them to output from the getmodelprevs function
getnsduhprevs <- function(mdseprevs,whichgender,numpop,denompop){
  nsduhdata = NULL
  nsduhdata <- melt(mdseprevs,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev") )
  nsduhdata <- cast(nsduhdata, group~survey_year,mean,subset=gender==whichgender&subpopulation==denompop&status==numpop)
  nsduhdata <- nsduhdata[c(-1)] # remove group name column
  row.names(nsduhdata)<-agerownames
  return(nsduhdata)
}

getnsduhprevsCI <- function(mdseprevs,whichgender,numpop,denompop){
  nsduhdatalow = NULL
  nsduhdatalow <- melt(mdseprevs,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev_lowCI") )
  nsduhdatalow <- cast(nsduhdatalow, group~survey_year,mean,subset=gender==whichgender&subpopulation==denompop&status==numpop)
  nsduhdatalow <- nsduhdatalow[c(-1)] # remove group name column
  row.names(nsduhdatalow)<-agerownames
  
  nsduhdatahigh = NULL
  nsduhdatahigh <- melt(mdseprevs,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev_highCI") )
  nsduhdatahigh <- cast(nsduhdatahigh, group~survey_year,mean,subset=gender==whichgender&subpopulation==denompop&status==numpop)
  nsduhdatahigh <- nsduhdatahigh[c(-1)] # remove group name column
  row.names(nsduhdatahigh)<-agerownames
  return(list(nsduhdatalow, nsduhdatahigh))
}

# Get model sum of squares ------------------------------------------------
getsumdiffs = function(out, whichgender){
  modelsmkdata=out[[1]]
  modeldepdata=out[[2]]
  years<- c("X2005","X2006","X2007","X2008","X2009","X2010","X2011","X2012","X2013","X2014","X2015","X2016","X2017","X2018","X2019","X2020") # only look at output for years where NSDUH data are available
  years2 <- c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020")

  smkdiffs <- rowSums(subset(as.data.frame(modelsmkdata[1]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender, "neversmoker","nevdeppop"),select=years2))^2 + # sum of squares
    rowSums(subset(as.data.frame(modelsmkdata[2]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"currentsmoker","nevdeppop"),select=years2))^2 +
    rowSums(subset(as.data.frame(modelsmkdata[3]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"formersmoker","nevdeppop"),select=years2))^2 +

    rowSums(subset(as.data.frame(modelsmkdata[4]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"neversmoker","notdeppop"),select=years2))^2 +
    rowSums(subset(as.data.frame(modelsmkdata[5]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"currentsmoker","notdeppop"),select=years2))^2 +
    rowSums(subset(as.data.frame(modelsmkdata[6]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"formersmoker","notdeppop"),select=years2))^2 +
    rowSums(subset(as.data.frame(modelsmkdata[7]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"neversmoker","deppop"),select=years2))^2 +
    rowSums(subset(as.data.frame(modelsmkdata[8]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"currentsmoker","deppop"),select=years2))^2 +
    rowSums(subset(as.data.frame(modelsmkdata[9]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"formersmoker","deppop"),select=years2))^2
  depdiffs <- rowSums(subset(as.data.frame(modeldepdata[1]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"nevdep","neversmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[2]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"notdep","neversmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[3]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"dep","neversmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[4]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"nevdep","currentsmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[5]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"notdep","currentsmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[6]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"dep","currentsmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[7]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"nevdep","formersmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[8]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"notdep","formersmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[9]),select=years) - subset(getnsduhprevs(mdseprevs,whichgender,"dep","formersmokers"),select=years2))^2

  return(c(sum(smkdiffs),sum(depdiffs)))
}


getsumdiffs_util = function(out, whichgender,util){
  if (whichgender=="females"){
    ahltmde <-	as.numeric(util)*c(0.4955417, 0.6229463, 0.7067596, 0.7326142, 0.6762736, 0.6576061) 
  }
  if (whichgender=="males"){
    ahltmde <- as.numeric(util)*c(0.3750369,0.4482582, 0.5409039,0.6605244,0.5895799,0.5128584 ) 
  }
  modelsmkdata=out[[1]]
  txdiffs <- abs(subset(as.data.frame(modelsmkdata[13]),select="X2018")[1:6,] - ahltmde) #absolute value of diffs
  return(sum(txdiffs))
}

# Bhat estimation ---------------------------------------------------------

# xF <- list(label=paramsnamesF, est=paramsF,low=lowervectorF,upp=uppervectorF) # est = parameter starting values
# xM <- list(label=paramsnamesM, est=paramsM,low=lowervectorM,upp=uppervectorM) # est = parameter starting values
#  
# ML_bhatF=function(paramsF){
#   out = main(getmodelprevs,"females",allparamsF,paramsF,paramsnamesF, mpc = 0, txeffcess = 1.0,util="1.0")
#   # LL = sum(getsumdiffs(out,"females"))   #Least squares
#   LL = sum(getsumdiffs_util(out,"females",util="1.0"))   #Least squares
#   cat(LL,paramsF,'\n')
#   return(LL)
# }
# 
# # resbhatF=dfp(xF,ML_bhatF)
# # paramsF=resbhatF$est
#  
# ML_bhatM=function(paramsM){
#   out = main(getmodelprevs,"males",allparamsM,paramsM,paramsnamesM, mpc = 0, txeffcess = 1.0,util="1.0")
#   # LL = sum(getsumdiffs(out,"males"))   #Least squares
#   LL = sum(getsumdiffs_util(out,"males",util="1.0"))
#   cat(LL,paramsM,'\n')
#   return(LL)
# }

# resbhatM=dfp(xM,ML_bhatM)
# paramsM=resbhatM$est

dir.create(file.path(mainDir, folder), showWarnings = FALSE)
setwd(file.path(mainDir,folder)) # save all output to this subdirectory

## Get model data following calibration 
outF = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF,mpc=0,txeffcess=1.0,util="1.0" )
outM = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, mpc=0,txeffcess=1.0,util="1.0" )
years2 <- c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020")

ns_deppopF = cbind(as.data.frame(outF[[1]][[1]][,years2]),"neversmoker","deppop","Women")
cs_deppopF = cbind(as.data.frame(outF[[1]][[2]][,years2]),"currentsmoker","deppop","Women")
fs_deppopF = cbind(as.data.frame(outF[[1]][[3]][,years2]),"formersmoker","deppop","Women")
ns_totalpopF = cbind(as.data.frame(outF[[1]][[4]][,years2]),"neversmoker","totalpop","Women")
cs_totalpopF = cbind(as.data.frame(outF[[1]][[5]][,years2]),"currentsmoker","totalpop","Women")
fs_totalpopF = cbind(as.data.frame(outF[[1]][[6]][,years2]),"formersmoker","totalpop","Women")

dep_nspopF = cbind(as.data.frame(outF[[2]][[1]][,years2]),"dep","neversmokers","Women") # dep prevalence among never smokers
dep_cspopF = cbind(as.data.frame(outF[[2]][[2]][,years2]),"dep","currentsmokers","Women")  # dep prevalence among current smokers
dep_fspopF = cbind(as.data.frame(outF[[2]][[3]][,years2]),"dep","formersmokers","Women")  # dep prevalence among former smokers
dep_totalpopF = cbind(as.data.frame(outF[[2]][[4]][,years2]),"dep","totalpop","Women")

ns_deppopM = cbind(as.data.frame(outM[[1]][[1]][,years2]),"neversmoker","deppop","Men")
cs_deppopM = cbind(as.data.frame(outM[[1]][[2]][,years2]),"currentsmoker","deppop","Men")
fs_deppopM = cbind(as.data.frame(outM[[1]][[3]][,years2]),"formersmoker","deppop","Men")
ns_totalpopM = cbind(as.data.frame(outM[[1]][[4]][,years2]),"neversmoker","totalpop","Men")
cs_totalpopM = cbind(as.data.frame(outM[[1]][[5]][,years2]),"currentsmoker","totalpop","Men")
fs_totalpopM = cbind(as.data.frame(outM[[1]][[6]][,years2]),"formersmoker","totalpop","Men")

dep_nspopM = cbind(as.data.frame(outM[[2]][[1]][,years2]),"dep","neversmokers","Men") # dep prevalence among never smokers
dep_cspopM = cbind(as.data.frame(outM[[2]][[2]][,years2]),"dep","currentsmokers","Men")  # dep prevalence among current smokers
dep_fspopM = cbind(as.data.frame(outM[[2]][[3]][,years2]),"dep","formersmokers","Men")  # dep prevalence among former smokers
dep_totalpopM = cbind(as.data.frame(outM[[2]][[4]][,years2]),"dep","totalpop","Men")

dfs_list = list(ns_deppopF,cs_deppopF,fs_deppopF,ns_totalpopF,cs_totalpopF,fs_totalpopF,dep_nspopF,dep_cspopF,dep_fspopF,dep_totalpopF,
                ns_deppopM,cs_deppopM,fs_deppopM,ns_totalpopM,cs_totalpopM,fs_totalpopM,dep_nspopM,dep_cspopM,dep_fspopM,dep_totalpopM)
dfs_list = lapply(dfs_list,function(x) cbind(x,agerownames))

outdata = data.frame(rbindlist(dfs_list,use.names=FALSE))
colnames(outdata) <- c(2005:2020,"status","subpopulation","gender", "age")
model = reshape2::melt(outdata,id.vars=c("age","status","subpopulation","gender"))
model$variable <- as.numeric(levels(model$variable))[model$variable] # converts factor to numeric

Figcalibdata <- function(whichgender,groupstatus,popgroup){
  
  nsduhdata <- subset(mdseprevs,gender==whichgender & status==groupstatus & subpopulation==popgroup)
  modeldata <- subset(model,gender==whichgender & status==groupstatus & subpopulation==popgroup)
  
  fig <-  ggplot() + geom_line(data = modeldata, aes(x=variable, y= value, colour=age))+
    geom_pointrange(data=nsduhdata, aes(x = survey_year, y = prev, colour=age, ymin=prev_lowCI, ymax=prev_highCI))+
    scale_y_continuous(name="Prevalence (%)",limits=c(0,1),breaks=seq(0,1,0.1)) +
    scale_x_continuous(name="Year",limits=c(2005,2020),breaks=seq(2005,2020, 1))  +
    labs(title= paste0(groupstatus, "prev, ",popgroup,", ",whichgender))+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.position="bottom",legend.title = element_blank())
  
  figtotal <-  ggplot() + geom_line(data = subset(modeldata,age=="total"), aes(x=variable, y= value))+
    geom_pointrange(data=subset(nsduhdata,age=="total"), aes(x = survey_year, y = prev, ymin=prev_lowCI, ymax=prev_highCI))+
    scale_y_continuous(name="Prevalence (%)",limits=c(0,1),breaks=seq(0,1,0.1)) +
    scale_x_continuous(name="Year",limits=c(2005,2020),breaks=seq(2005,2020, 1))  +
    labs(title= paste0(groupstatus, "prev, ",popgroup,", ",whichgender))+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.position="bottom",legend.title = element_blank())
  return(list(fig,figtotal))
}

a <- Figcalibdata("Women","neversmoker","deppop")
b <- Figcalibdata("Women","currentsmoker","deppop")
c <- Figcalibdata("Women","formersmoker","deppop")

d <- Figcalibdata("Women","neversmoker","totalpop")
e <- Figcalibdata("Women","currentsmoker","totalpop")
f <- Figcalibdata("Women","formersmoker","totalpop")

g <- Figcalibdata("Women","dep","neversmokers")
h <- Figcalibdata("Women","dep","currentsmokers")
i <- Figcalibdata("Women","dep","formersmokers")
j <- Figcalibdata("Women","dep","totalpop")

k <- Figcalibdata("Men","neversmoker","deppop")
l <- Figcalibdata("Men","currentsmoker","deppop")
m <- Figcalibdata("Men","formersmoker","deppop")

n <- Figcalibdata("Men","neversmoker","totalpop")
o <- Figcalibdata("Men","currentsmoker","totalpop")
p <- Figcalibdata("Men","formersmoker","totalpop")

q <- Figcalibdata("Men","dep","neversmokers")
r <- Figcalibdata("Men","dep","currentsmokers")
s <- Figcalibdata("Men","dep","formersmokers")
t <- Figcalibdata("Men","dep","totalpop")

pdf(file = paste0(folder,"modelcalib_", format(Sys.Date(),"%m-%d-%y"),".pdf"),width=10, height=6, onefile=TRUE)
plot.new()
text(.5, 0.9, "Calibrated model - Women", font=2, cex=1.5)
# grid.table(c(paramsF,getsumdiffs(outF,"females")),rows=c(paramsnamesF,"Smk Fit", "Dep Fit"))
grid.arrange(a[[1]], b[[1]],c[[1]],ncol=3)
grid.arrange(a[[2]], b[[2]],c[[2]],ncol=3)
grid.arrange(d[[1]], e[[1]],f[[1]],ncol=3)
grid.arrange(d[[2]], e[[2]],f[[2]],ncol=3)
grid.arrange(g[[1]], h[[1]],i[[1]],ncol=3)
grid.arrange(j[[1]], j[[2]],ncol=2)
plot.new()
text(.5, 0.9, "Calibrated model - Men", font=2, cex=1.5)
# grid.table(c(paramsM,getsumdiffs(outM,"males")),rows=c(paramsnamesM,"Smk Fit", "Dep Fit"))
grid.arrange(k[[1]], l[[1]],m[[1]],ncol=3)
grid.arrange(k[[2]], l[[2]],m[[2]],ncol=3)
grid.arrange(n[[1]], o[[1]],p[[1]],ncol=3)
grid.arrange(n[[2]], o[[2]],p[[2]],ncol=3)
grid.arrange(q[[1]], r[[1]],s[[1]],ncol=3)
grid.arrange(t[[1]], t[[2]],ncol=2)
dev.off()
