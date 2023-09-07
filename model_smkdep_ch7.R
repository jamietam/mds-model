rm(list = ls())
# mainDir <- "C:/Users/mauro/Dropbox/SGR Chapter 7 Mental Health Substance Use/Data analysis/mds-model"
mainDir <- "C:/Users/JT936/Dropbox/SGR Chapter 7 Mental Health Substance Use/Data Analysis/mds-model"

setwd(file.path(mainDir))
library(openxlsx)
library(reshape)
library(splines)
library(Bhat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
load("mdseprevs0520.rda")
depsmkprevs_by_year = subset(mdseprevs, survey_year<2020)


date = "083023" # name the folder where results will be saved
namethisrun = "_checkSADprops" # name this run

folder = paste0(date,namethisrun,"/") # name the folder where results will be saved

# Read inputs -------------------------------------------------------------
startyear = 1900 # the burn-in period starting point
endyear = 2100 # 2014 or 2050 , still need census projections through 2065
startage = 0
endage = 99
Ny= endyear - startyear + 1 
Na= endage - startage + 1
emptycompartment <- matrix(0, nrow = Na, ncol = Ny, dimnames=list(c(startage:endage),c(startyear:endyear))) # create matrix of zeroes for compartments
policystart = 2024

allparamsF = read.xlsx("parameters_tx_v2.xlsx",sheet=paste0("model_females"),rowNames=TRUE,colNames=TRUE) # Adjust parameters in this excel file with estimates to be used in the model
paramsF = as.vector(subset(allparamsF,bhat==1)[['estimate']]) # Parameters where bhat = 1 can be estimated by bhat
paramsnamesF = rownames(subset(allparamsF,bhat==1))  
lowervectorF = as.vector(subset(allparamsF,bhat==1)[['lower']]) # lower bounds for bhat parameters
uppervectorF = as.vector(subset(allparamsF,bhat==1)[['upper']]) # upper bounds for bhat parameters

allparamsM = read.xlsx("parameters_tx_v2.xlsx",sheet=paste0("model_males"),rowNames=TRUE,colNames=TRUE) # Adjust parameters in this excel file with estimates to be used in the model
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

getnsduhprevs <- function(depsmkprevs_by_year,whichgender,numpop,denompop){
  nsduhdata = NULL
  nsduhdata <- reshape2::melt(depsmkprevs_by_year,id.vars=c("group","survey_year","age","sex","subpopulation","status"),measure.vars = c("prev") )
  nsduhdata <- cast(nsduhdata, group~survey_year,mean,subset=sex==whichgender&subpopulation==denompop&status==numpop)
  nsduhdata <- nsduhdata[c(-1)] # remove group name column
  row.names(nsduhdata)<-agerownames
  return(nsduhdata)
}
# Main model --------------------------------------------------------------
main <- function(getmodelprevs, whichgender, allparamsF, paramsF,paramsnamesF, 
                 initeff_dep,initeff_ndep,cesseff_dep, cesseff_ndep, equity_init, equity_cess){
  setwd(file.path(mainDir))
  
  pop = read.xlsx("census_data/np2017_d1.xlsx",sheet=whichgender,rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  
  # smk params --------------------------------------------------------------
  
  # smk_init_cisnet = read.xlsx("cisnet_smkrates_nhis2018.xlsx",sheet=paste0(whichgender,"_init"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  # smk_cess_cisnet = read.xlsx("cisnet_smkrates_nhis2018.xlsx",sheet=paste0(whichgender,"_cess"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  smk_init_cisnet = read.xlsx("cisnet_smkrates_by_calyear_061923.xlsx",sheet=paste0(whichgender,"_init"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  smk_cess_cisnet = read.xlsx("cisnet_smkrates_by_calyear_061923.xlsx",sheet=paste0(whichgender,"_cess"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  
  smkinit_SF = matrix(0,100,1)
  smkinit_SF[1:18] = ifelse(allparamsF["smkinit_youthSF","bhat"]==0,allparamsF["smkinit_youthSF","estimate"],paramsF[match("smkinit_youthSF",paramsnamesF)])
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
  
  death_ns = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("ns_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  death_cs = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("cs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  death_fs = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("fs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  
  LE_ns_dep = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_ns_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  LE_cs_dep = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_cs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  LE_fs_dep = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_fs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
 
  LE_ns_nodep = LE_ns_dep
  LE_cs_nodep = LE_cs_dep
  LE_fs_nodep = LE_fs_dep
  
  dep1inc = read.xlsx("incidence_eaton.xlsx",sheet=whichgender,rowNames=TRUE, colNames=FALSE, check.names=FALSE)
  
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
  
  # Compartments / state variables ------------------------------------------
  
  # Initialize population - Model compartments are organized with age 0-99 as rows, and year as columns
  matrix.names<-c('ns_nevdep', 'ns_dep1','ns_fdep', 'ns_dep2','ns_recall',
                  'cs_nevdep', 'cs_dep1','cs_fdep', 'cs_dep2','cs_recall',
                  'fs_nevdep', 'fs_dep1','fs_fdep', 'fs_dep2','fs_recall','ns_deaths', 'cs_deaths', 'fs_deaths')
  for (name in matrix.names) assign(name,emptycompartment)
  
  ns_nevdep[paste(startage),1:Ny] <- as.matrix(pop[paste(startage),1:Ny])   # Takes empty compartment and populates the top row of the matrix with the number of 0-yrolds
  
  # Run it ------------------------------------------------------------------
  
  for (y in c((startyear+1):(endyear))){
    py = paste(y - 1)
    if (y<policystart){
      txeffdep_init = 1
      txeffndep_init = 1
      txeffdep_cess = 1
      txeffndep_cess = 1
    } else {
      txeffdep_init = initeff_dep*equity_init
      txeffndep_init = initeff_ndep
      txeffdep_cess = cesseff_dep*equity_cess
      txeffndep_cess = cesseff_ndep
    }
    if (y>=2016 & y <= 2100){
      usethis = scaleddep1_inc # from 2016-2100, use scaled incidence probabilities
    } else {
      usethis = dep1_inc # all other years, use the original incidence probabilities						
    }
    
    # Never MDE
    ns_nevdep[2:Na,paste(y)] <- ns_nevdep[1:Na-1,py]*(1-txeffndep_init*smk_init[(startage+1):endage,py])*(1-usethis[(startage+1):(endage),])*(1-death_ns[(startage+1):endage,py])
    cs_nevdep[2:Na,paste(y)] <- cs_nevdep[1:Na-1,py]*(1-txeffndep_cess*smk_cess[(startage+1):endage,py])*(1-RRcs_dep1*usethis[(startage+1):(endage),])*(1-death_cs[(startage+1):endage,py]) + ns_nevdep[1:Na-1,py]*(txeffndep_init*smk_init[(startage+1):endage,py])
    fs_nevdep[2:Na,paste(y)] <- fs_nevdep[1:Na-1,py]*(1-RRfs_dep1*usethis[(startage+1):(endage),])*(1-death_fs[(startage+1):endage,py]) + cs_nevdep[1:Na-1,py]*txeffndep_cess*(smk_cess[(startage+1):endage,py])
    
    # Current MDE (1st episode)
    ns_dep1[2:Na,paste(y)] <- ns_dep1[1:Na-1,py]*(1-Edepr_smkinit*txeffdep_init*smk_init[(startage+1):endage,py])*(1-deprecov_rate)*(1-RRdepr_death*death_ns[(startage+1):endage,py]) + ns_nevdep[1:Na-1,py]*(usethis[(startage+1):(endage),])
    cs_dep1[2:Na,paste(y)] <- cs_dep1[1:Na-1,py]*(1-ORhdep_quit*txeffdep_cess*smk_cess[(startage+1):endage,py])*(1-deprecovSF_cs*deprecov_rate)*(1-RRdepr_death*death_cs[(startage+1):endage,py]) +  cs_nevdep[1:Na-1,py]*(RRcs_dep1*usethis[(startage+1):(endage),]) + ns_dep1[1:Na-1,py]*(Edepr_smkinit*txeffdep_init*smk_init[(startage+1):endage,py])
    fs_dep1[2:Na,paste(y)] <- fs_dep1[1:Na-1,py]*(1-deprecovSF_fs*deprecov_rate)*(1-RRdepr_death*death_fs[(startage+1):endage,py]) + fs_nevdep[1:Na-1,py]*(RRfs_dep1*usethis[(startage+1):(endage),]) + cs_dep1[1:Na-1,py]*(ORhdep_quit*txeffdep_cess*smk_cess[(startage+1):endage,py])
    
    # Current MDE (Recurrent episode)
    ns_dep2[2:Na,paste(y)] <- ns_dep2[1:Na-1,py]*(1-Edepr_smkinit*txeffdep_init*smk_init[(startage+1):endage,py])*(1-deprecov_rate)*(1-RRdepr_death*death_ns[(startage+1):endage,py]) + ns_fdep[1:Na-1,py]*(depr_inc[(startage+1):(endage),]) + ns_recall[1:Na-1,py]*(depr_inc[(startage+1):(endage),])
    cs_dep2[2:Na,paste(y)] <- cs_dep2[1:Na-1,py] *(1-ORhdep_quit*txeffdep_cess*smk_cess[(startage+1):endage,py])*(1-deprecovSF_cs*deprecov_rate)*(1-RRdepr_death*death_cs[(startage+1):endage,py]) + cs_fdep[1:Na-1,py]*(Ecs_depr*depr_inc[(startage+1):(endage),]) + ns_dep2[1:Na-1,py]*(Edepr_smkinit*txeffdep_init*smk_init[(startage+1):endage,py]) + cs_recall[1:Na-1,py]*(depr_inc[(startage+1):(endage),])
    fs_dep2[2:Na,paste(y)] <- fs_dep2[1:Na-1,py]*(1-deprecovSF_fs*deprecov_rate)*(1-RRdepr_death*death_fs[(startage+1):endage,py]) + fs_fdep[1:Na-1,py]*(Efs_depr*depr_inc[(startage+1):(endage),]) + cs_dep2[1:Na-1,py]*(ORhdep_quit*txeffdep_cess*smk_cess[(startage+1):endage,py]) + fs_recall[1:Na-1,py]*(depr_inc[(startage+1):(endage),])
    
    # Former MDE - Recovery with risk of recurrence  
    ns_fdep[2:Na,paste(y)] <- ns_fdep[1:Na-1,py]*(1-txeffndep_init*smk_init[(startage+1):endage,py])*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_ns[(startage+1):endage,py])*(1-forget[(startage+1):(endage)]) + ns_dep1[1:Na-1,py]*(deprecov_rate) + ns_dep2[1:Na-1,py]*(deprecov_rate)
    cs_fdep[2:Na,paste(y)] <- cs_fdep[1:Na-1,py]*(1-ORhdep_quit*txeffndep_cess*smk_cess[(startage+1):endage,py])*(1-Ecs_depr*depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_cs[(startage+1):endage,py])*(1-forget[(startage+1):(endage)])  + cs_dep1[1:Na-1,py]*(deprecovSF_cs*deprecov_rate) + cs_dep2[1:Na-1,py]*(deprecovSF_cs*deprecov_rate) + ns_fdep[1:Na-1,py]*(txeffndep_init*smk_init[(startage+1):endage,py])
    fs_fdep[2:Na,paste(y)] <- fs_fdep[1:Na-1,py]*(1-Efs_depr*depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_fs[(startage+1):endage,py])*(1-forget[(startage+1):(endage)]) + fs_dep1[1:Na-1,py]*(deprecovSF_fs*deprecov_rate) + fs_dep2[1:Na-1,py]*(deprecovSF_fs*deprecov_rate) + cs_fdep[1:Na-1,py]*(ORhdep_quit*txeffndep_cess*smk_cess[(startage+1):endage,py])
    
    # Recall Error 
    ns_recall[2:Na,paste(y)] <- ns_recall[1:Na-1,py]*(1-txeffndep_init*smk_init[(startage+1):endage,py])*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_ns[(startage+1):endage,py]) + ns_fdep[1:Na-1,py]*(forget[(startage+1):(endage)])
    cs_recall[2:Na,paste(y)] <- cs_recall[1:Na-1,py]*(1-ORhdep_quit*txeffndep_cess*smk_cess[(startage+1):endage,py])*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_cs[(startage+1):endage,py]) + ns_recall[1:Na-1,py]*(txeffndep_init*smk_init[(startage+1):endage,py]) + cs_fdep[1:Na-1,py]*(forget[(startage+1):(endage)])
    fs_recall[2:Na,paste(y)] <- fs_recall[1:Na-1,py]*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_fs[(startage+1):endage,py])+ cs_recall[1:Na-1,py]*(ORhdep_quit*txeffndep_cess*smk_cess[(startage+1):endage,py]) + fs_fdep[1:Na-1,py]*(forget[(startage+1):(endage)])
    
  }
  
  # Get inputs --------------------------------------------------------------
  
  initiation <- data.frame(cbind(smk_init_cisnet[,"2019"],smk_init_cisnet[,"2019"],smkinit_SF,smk_init["2019"],c(0:99),c(Edepr_smkinit,Edepr_smkinit[99])))
  colnames(initiation) <- c("2019", "cisnet2019","SF","scaledrates","age","Edepr_smkinit")
  
  cessation <- data.frame(cbind(smk_cess_cisnet[,"2019"],smk_cess_cisnet[,"2019"],smkcess_SF,smk_cess["2019"],c(0:99),c(ORhdep_quit,Edepr_smkinit[99])))
  colnames(cessation) <- c("2019", "cisnet2019","SF","scaledrates","age","ORhdep_quit")
  
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
  
  deppop = ns_dep+cs_dep+fs_dep
  
  ns_ndep <- (ns_nevdep+ns_fdep+ns_recall) 
  cs_ndep <- (cs_nevdep+cs_fdep+cs_recall) 
  fs_ndep <- (fs_nevdep+fs_fdep+fs_recall) 
  
  ndeppop = ns_ndep+cs_ndep+fs_ndep
  

  # Smoking attributable mortality
  SADdep = cs_dep *(death_cs[,paste(c(startyear:endyear))]-death_ns[,paste(c(startyear:endyear))])+fs_dep * (death_fs[,paste(c(startyear:endyear))]-death_ns[,paste(c(startyear:endyear))])
  SADndep = cs_ndep *(death_cs[,paste(c(startyear:endyear))]-death_ns[,paste(c(startyear:endyear))])+fs_ndep * (death_fs[,paste(c(startyear:endyear))]-death_ns[,paste(c(startyear:endyear))])
  ALLdeathsdep = ns_dep * death_ns[,paste(c(startyear:endyear))] +cs_dep *death_cs[,paste(c(startyear:endyear))] + fs_dep * death_fs[,paste(c(startyear:endyear))]
  ALLdeathsndep = ns_ndep * death_ns[,paste(c(startyear:endyear))] +cs_ndep *death_cs[,paste(c(startyear:endyear))] + fs_ndep * death_fs[,paste(c(startyear:endyear))]
  
  browser()
  data = cbind(cs_dep[,"2023"],cs_dep[,"2040"],cs_dep[,"2060"],cs_dep[,"2080"],cs_dep[,"2100"],  fs_dep[,"2023"],fs_dep[,"2040"],fs_dep[,"2060"],fs_dep[,"2080"],fs_dep[,"2100"],
  cs_ndep[,"2023"],cs_ndep[,"2040"],cs_ndep[,"2060"],cs_ndep[,"2080"],cs_ndep[,"2100"],fs_ndep[,"2023"],fs_ndep[,"2040"],fs_ndep[,"2060"],fs_ndep[,"2080"],fs_ndep[,"2100"],
  ns_dep[,"2023"],ns_dep[,"2040"],ns_dep[,"2060"],ns_dep[,"2080"],ns_dep[,"2100"],ns_ndep[,"2023"],ns_ndep[,"2040"],ns_ndep[,"2060"],ns_ndep[,"2080"],ns_ndep[,"2100"]  )
  colnames(data) <- c("cs_dep2023","2040","2060","2080","2100","fs_dep2023","2040","2060","2080","2100", "cs_ndep2023","2040","2060","2080","2100","fs_ndep2023","2040","2060","2080","2100",
                      "ns_dep2023","2040","2060","2080","2100","ns_ndep2023","2040","2060","2080","2100")
  write.csv(data,"popsizesF.csv")
  # 
  # depdata <- cbind(cs_dep[,c("2023","2100")],fs_dep[,"2023"],ns_dep[,c("2023","2100")death_cs[,"2023"],death_ns[,"2023"],
  # Years of Life Lost
  YLLdep =  LE_ns_dep[,paste(c(2023:endyear))] * (cs_dep[,paste(c(2023:endyear))]*ucs_minus_uns[,paste(c(2023:endyear))] + fs_dep[,paste(c(2023:endyear))]*ufs_minus_uns[,paste(c(2023:endyear))])
  YLLndep = LE_ns_nodep[,paste(c(2023:endyear))] * (cs_nevdep[,paste(c(2023:endyear))]*ucs_minus_uns[,paste(c(2023:endyear))] + fs_nevdep[,paste(c(2023:endyear))]*ufs_minus_uns[,paste(c(2023:endyear))])

  # Get prevalences for model fitting ---------------------------------------

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
  modeldepdata <- list(d1,d2,d3,d4)
  
  PRequity = (s2["total","2100"]/s7["total" ,"2100"]-1)^2
  
  return(list(modelsmkprevdata, PRequity, SADdep, SADndep, ALLdeathsdep, ALLdeathsndep, YLLdep, YLLndep))
  #return(list(modelsmkprevdata, modeldepdata)) #PRequity, SADdep, SADndep, ALLdeathsdep, ALLdeathsndep, YLLdep, YLLndep))
}

# Bhat estimation ---------------------------------------------------------

getsumdiffs = function(out, assignedsex){
  modelsmkdata=out[[1]]
  modeldepdata=out[[2]]
  years<- c("X2005","X2006","X2007","X2008","X2009","X2010","X2011","X2012","X2013","X2014","X2015","X2016","X2017","X2018","X2019") # only look at output for years where NSDUH data are available
  years2 <- c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019")
  
  smkdiffs <- 
    rowSums(subset(as.data.frame(modelsmkdata[1]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"neversmoker","deppop"),select=years2))^2 +
    rowSums(subset(as.data.frame(modelsmkdata[2]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"currentsmoker","deppop"),select=years2))^2 +    
    rowSums(subset(as.data.frame(modelsmkdata[3]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"formersmoker","deppop"),select=years2))^2+    
    rowSums(subset(as.data.frame(modelsmkdata[4]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"neversmoker","totalpop"),select=years2))^2+    
    rowSums(subset(as.data.frame(modelsmkdata[5]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"currentsmoker","totalpop"),select=years2))^2+    
    rowSums(subset(as.data.frame(modelsmkdata[6]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"formersmoker","totalpop"),select=years2))^2
  
  depdiffs <- 
    rowSums(subset(as.data.frame(modeldepdata[1]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"dep","neversmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[2]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"dep","currentsmokers"),select=years2))^2 +
    rowSums(subset(as.data.frame(modeldepdata[3]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"dep","formersmokers"),select=years2))^2+
    rowSums(subset(as.data.frame(modeldepdata[4]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,assignedsex,"dep","totalpop"),select=years2))^2
  
  return(c(sum(smkdiffs),sum(depdiffs)))
}

xF <- list(label=paramsnamesF, est=paramsF,low=lowervectorF,upp=uppervectorF) # est = parameter starting values
xM <- list(label=paramsnamesM, est=paramsM,low=lowervectorM,upp=uppervectorM) # est = parameter starting values

ML_bhatF=function(paramsF){
  out = main(getmodelprevs,"females",allparamsF,paramsF,paramsnamesF,1,1,1,1,paramsF[[1]],paramsF[[2]])
  #LL = sum(getsumdiffs(out,"females"))
  LL = out[[2]]
  cat(LL,paramsF,'\n')
  return(LL)
}
resbhatF=dfp(xF,ML_bhatF)
paramsF=resbhatF$est


ML_bhatM=function(paramsM){
  out = main(getmodelprevs,"males",allparamsM,paramsM,paramsnamesM,1,1,1,1,paramsM[[1]],paramsM[[2]])
  #LL = sum(getsumdiffs(out,"males"))   #Least squares
  LL = out[[2]]
  cat(LL,paramsM,'\n')
  return(LL)
}
resbhatM=dfp(xM,ML_bhatM)
paramsM=resbhatM$est

#
# Calibration visualization -------------------------------------------------------

dir.create(file.path(mainDir, folder), showWarnings = FALSE)
setwd(file.path(mainDir,folder)) # save all output to this subdirectory

## Get model data following calibration 
outF = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, 1,1,1,1,1,1)
outM = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, 1,1,1,1,1,1)
years2 <- c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019")
  
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
colnames(outdata) <- c(2005:2019,"status","subpopulation","gender", "age")
model = reshape2::melt(outdata,id.vars=c("age","status","subpopulation","gender"))
model$variable <- as.numeric(levels(model$variable))[model$variable] # converts factor to numeric

Figcalibdata <- function(whichgender,groupstatus,popgroup){

  nsduhdata <- subset(depsmkprevs_by_year,gender==whichgender & status==groupstatus & subpopulation==popgroup)
  modeldata <- subset(model,gender==whichgender & status==groupstatus & subpopulation==popgroup)
  
  fig <-  ggplot() + geom_line(data = modeldata, aes(x=variable, y= value, colour=age))+
    geom_pointrange(data=nsduhdata, aes(x = survey_year, y = prev, colour=age, ymin=prev_lowCI, ymax=prev_highCI))+
    scale_y_continuous(name="Prevalence (%)",limits=c(0,1),breaks=seq(0,1,0.1)) +
    scale_x_continuous(name="Year",limits=c(2005,2019),breaks=seq(2005,2019, 1))  +
    labs(title= paste0(groupstatus, "prev, ",popgroup,", ",whichgender))+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.position="bottom",legend.title = element_blank())
  
  figtotal <-  ggplot() + geom_line(data = subset(modeldata,age=="total"), aes(x=variable, y= value))+
    geom_pointrange(data=subset(nsduhdata,age=="total"), aes(x = survey_year, y = prev, ymin=prev_lowCI, ymax=prev_highCI))+
    scale_y_continuous(name="Prevalence (%)",limits=c(0,1),breaks=seq(0,1,0.1)) +
    scale_x_continuous(name="Year",limits=c(2005,2019),breaks=seq(2005,2019, 1))  +
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

pdf(file = paste0(folder,"SGR_smkdep_", format(Sys.Date(),"%m-%d-%y"),".pdf"),width=10, height=6, onefile=TRUE)
plot.new()
text(.5, 0.9, "Calibration parameters - Women", font=2, cex=1.5)
grid.table(c(paramsF,getsumdiffs(outF,"females")),rows=c(paramsnamesF,"Smk Fit", "Dep Fit"))
grid.arrange(a[[1]], b[[1]],c[[1]],ncol=3)
grid.arrange(a[[2]], b[[2]],c[[2]],ncol=3)
grid.arrange(d[[1]], e[[1]],f[[1]],ncol=3)
grid.arrange(d[[2]], e[[2]],f[[2]],ncol=3)
grid.arrange(g[[1]], h[[1]],i[[1]],ncol=3)
grid.arrange(j[[1]], j[[2]],ncol=2)
plot.new()
text(.5, 0.9, "Calibration parameters - Men", font=2, cex=1.5)
grid.table(c(paramsM,getsumdiffs(outM,"males")),rows=c(paramsnamesM,"Smk Fit", "Dep Fit"))
grid.arrange(k[[1]], l[[1]],m[[1]],ncol=3)
grid.arrange(k[[2]], l[[2]],m[[2]],ncol=3)
grid.arrange(n[[1]], o[[1]],p[[1]],ncol=3)
grid.arrange(n[[2]], o[[2]],p[[2]],ncol=3)
grid.arrange(q[[1]], r[[1]],s[[1]],ncol=3)
grid.arrange(t[[1]], t[[2]],ncol=2)
dev.off()


#### added 0318
dir.create(file.path(mainDir, folder), showWarnings = FALSE)
setwd(file.path(mainDir,folder)) # save all output to this subdirectory


txeffset = list(c(1,1,1,1,1,1,1,1), 
                c(0.8,0.8,1.2,1.2,1,1,1,1),
                c(1,1,1,1,0.23456919,6.4208368,0.53074971, 3.2352955))
scenarios = c("Status Quo", "Equality","Equity")


getprevsFM <- function(initeff_dep,initeff_ndep,cesseff_dep, cesseff_ndep, 
                       equity_init_F, equity_cess_F, equity_init_M, equity_cess_M){
  
  outF = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, 
              initeff_dep,initeff_ndep,cesseff_dep, cesseff_ndep, equity_init_F, equity_cess_F)
  
  cs_ndeppopF = outF[[1]][[7]]
  cs_deppopF = outF[[1]][[2]]
  
  outM = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, 
              initeff_dep,initeff_ndep,cesseff_dep, cesseff_ndep, equity_init_M, equity_cess_M)
  
  cs_ndeppopM = outM[[1]][[7]]
  cs_deppopM = outM[[1]][[2]]
  
  return(list(cs_deppopF, cs_deppopM, cs_ndeppopF,cs_ndeppopM))
}

for(t in 1:length(txeffset)){ # Baseline and policy scenarios
  allprevs = getprevsFM(txeffset[[t]][1],txeffset[[t]][2],txeffset[[t]][3],txeffset[[t]][4],
                        txeffset[[t]][5],txeffset[[t]][6],txeffset[[t]][7],txeffset[[t]][8])
  assign(paste0("cs_deppop","F",(t-1)), allprevs[[1]])
  assign(paste0("cs_deppop","M",(t-1)), allprevs[[2]])
  assign(paste0("cs_ndeppop","F",(t-1)), allprevs[[3]])
  assign(paste0("cs_ndeppop","M",(t-1)), allprevs[[4]])
}


Figmodelnsduh <- function(whichgender,dfs,smkgroup,popgroup,popgroup2, title){
  dfs = lapply(dfs, function(x) {as.data.frame(x)})
  dfs = Map(cbind,dfs,scenario=as.list(scenarios))
  
  fullset = bind_rows(dfs, .id = NULL)
  fullset$age <-rep(agerownames,length(scenarios))
  fullset$status <- c(rep("Current MD",18),rep("No Current MD",18))
  model <- rbind(reshape2::melt(fullset,id.vars=c("age","scenario","status")))
  
  colnames(model)[5] <- "modelvalue"
  model$variable <- as.numeric(levels(model$variable))[model$variable] # converts factor to numeric
  model$scenario  <- factor(model$scenario , levels = scenarios[1:length(scenarios)] )
  
  nsduh1 <- getnsduhprevs(depsmkprevs_by_year,whichgender,smkgroup,popgroup)
  nsduh1$age<-agerownames
  nsduh1$status<-"Current MD"
  nsduh1 <- reshape2::melt(as.data.frame(nsduh1),id.vars=c("age","status"))
  
  nsduh2 <- getnsduhprevs(depsmkprevs_by_year,whichgender,smkgroup,popgroup2)
  nsduh2$age<-agerownames
  nsduh2$status<-"No Current MD"
  nsduh2 <- reshape2::melt(as.data.frame(nsduh2),id.vars=c("age","status"))
  
  nsduh <- rbind(nsduh1,nsduh2)
  colnames(nsduh)[4] <- c("nsduhvalue")
  nsduh$variable <- as.numeric(levels(nsduh$variable))[nsduh$variable] # converts factor to numeric
  
  fig <-  ggplot() + geom_line(data = subset(model,age=="total"), aes(x=variable, y= modelvalue*100, colour=status,linetype=scenario))+
    geom_point(data=subset(nsduh, age=="total"), aes(x = variable, y = nsduhvalue*100, colour=status, shape="National Survey on Drug Use and Health"))+
    # scale_color_viridis(discrete=TRUE, begin=0, end = 0.75, option="magma") +
    scale_color_grey(start=0, end = 0.75) +
    scale_y_continuous(name="Smoking prevalence (%)",limits=c(0,50),breaks=seq(0,50,5)) +
    scale_x_continuous(name="Year",limits=c(2005,2100),breaks=c(2005,seq(2020,2100, 10)))  +
    labs(title = title)+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
  # theme(axis.text.x=element_text(angle=60, hjust=1), legend.title = element_blank())
  return(fig)
}

dfsF = list(cs_deppopF0,cs_deppopF1,cs_deppopF2,cs_ndeppopF0, cs_ndeppopF1, cs_ndeppopF2)
dfsM = list(cs_deppopM0,cs_deppopM1,cs_deppopM2,cs_ndeppopM0, cs_ndeppopM1, cs_ndeppopM2)

Figure3a <- Figmodelnsduh("females",dfsF,"currentsmoker","deppop","nevdeppop", "A) Women with Current MD vs. No Current MD")
Figure3b <- Figmodelnsduh("males",dfsM,"currentsmoker","deppop","nevdeppop", "B) Men with Current MD vs. No Current MD")

library(grid)
library(gridExtra)
grid_arrange_shared_legend <- function(plots,columns,titletext) {
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(arrangeGrob(grobs= lapply(plots, function(x)
    x + theme(legend.position="none", plot.title = element_text(size = rel(0.8)))),ncol=columns),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight),
    top=textGrob(titletext,just="top", vjust=1,check.overlap=TRUE,gp=gpar(fontsize=9, fontface="bold"))
  )
}
#pdf(file = paste0("Modelsmkprev_current_non", namethisrun,".pdf"),width=10, height=6, onefile=FALSE)
jpeg(filename = paste0("Modelsmkprev_current_non",namethisrun, ".jpg"),width=13, height=6, units ="in", res=1000)

grid_arrange_shared_legend(list(Figure3a, Figure3b),2,"")
dev.off()

# RESULTS -----------------------------------------------------------------
outF0 = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, 1,1,1,1,1,1) # runs model using parameters specified in excel sheet OR using bhat estimates
outM0 = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, 1,1,1,1,1,1)

outF1 = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, 0.8,0.8,1.2,1.2,1,1) # runs model using parameters specified in excel sheet OR using bhat estimates
outM1 = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, 0.8,0.8,1.2,1.2,1,1)

outF2 = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, 1,1,1,1,0.23456919,6.4208368) # runs model using parameters specified in excel sheet OR using bhat estimates
outM2 = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, 1,1,1,1, 0.53074971,3.2352955)

#outF9 = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, mphr = 1) # runs model using parameters specified in excel sheet OR using bhat estimates
#outM9 = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, mphr = 1)



library(ggplot2)
library(reshape)
library(grid)
library(gridBase)
library(gridExtra)
library(plyr)
library(Hmisc)
library(scales)
library(xlsx)

xaxisbreaks = c(2023,seq(2023,2100,10)) # specify the ticks on the x-axis of your results plots
minyear = 2023
maxyear = 2100

dir.create(file.path(mainDir, date), showWarnings = FALSE)
setwd(file.path(mainDir, date)) # save all output to this subdirectory

# Read in each model's outputs ---------------------------------------------------
scenarios = c("Status quo", "Equality","Equity")

outlistF = list(outF0,outF1,outF2)
outlistM = list(outM0,outM1,outM2)

### Read in model data 
for (i in seq(0,length(outlistF)-1)){
  outF = outlistF[[i+1]]
  
  assign(paste0("modelsmkdata","F",i), outF[[1]])
  
  assign(paste0("SADdep","F",i), data.frame(outF[[3]],check.names=FALSE))
  assign(paste0("ALLdeathsdep","F",i), data.frame(outF[[5]],check.names=FALSE))
  assign(paste0("SADndep","F",i), data.frame(outF[[4]],check.names=FALSE))
  assign(paste0("ALLdeathsndep","F",i), data.frame(outF[[6]],check.names=FALSE))
  
  assign(paste0("YLLdep","F",i), data.frame(outF[[7]],check.names=FALSE))
  assign(paste0("YLLndep","F",i), data.frame(outF[[8]],check.names=FALSE))
  
  
  outM = outlistM[[i+1]]
  
  assign(paste0("modelsmkdata","M",i), outM[[1]])
  
  assign(paste0("SADdep","M",i), data.frame(outM[[3]],check.names=FALSE))
  assign(paste0("ALLdeathsdep","M",i), data.frame(outM[[5]],check.names=FALSE))
  assign(paste0("SADndep","M",i), data.frame(outM[[4]],check.names=FALSE))
  assign(paste0("ALLdeathsndep","M",i), data.frame(outM[[6]],check.names=FALSE))
  
  assign(paste0("YLLdep","M",i), data.frame(outM[[7]],check.names=FALSE))
  assign(paste0("YLLndep","M",i), data.frame(outM[[8]],check.names=FALSE))
  
  
  for (g in c("F","M")){
    # assign(paste0("ns_ndeppop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[1],check.names=FALSE)) ### WRONG ASSIGNMENT
    # assign(paste0("cs_ndeppop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[2],check.names=FALSE))
    # assign(paste0("fs_ndeppop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[3],check.names=FALSE))
    # assign(paste0("ns_deppop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[4],check.names=FALSE))
    # assign(paste0("cs_deppop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[5],check.names=FALSE))
    # assign(paste0("fs_deppop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[6],check.names=FALSE))
    # assign(paste0("cs_totalpop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[7],check.names=FALSE))
    assign(paste0("cs_deppop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[2],check.names=FALSE)) # s2 <- getmodelprevs(cs_dep,deppop)
    assign(paste0("cs_ndeppop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[7],check.names=FALSE)) # s7 <- getmodelprevs(cs_ndep,ndeppop)
    assign(paste0("cs_totalpop",g,i), data.frame(get(paste0("modelsmkdata",g,i))[5],check.names=FALSE)) # s5 <- getmodelprevs(cspop,totalpop)
  }
}



## Reproduce Table 1 in Tam (2020) Modeling smoking-attributable mortality among adults with major depression in the United States - Preventive Medicine 
# Table 1. ------------------------------------------------------------------------
scenariotable <- function(SADdep, SADndep,YLLdep, YLLndep, ALLdeathsdep, ALLdeathsndep,
                          cs_deppop,cs_ndeppop, cs_totalpop, g){
  
  totalSAD <- SADdep + SADndep
  totaldeaths <- ALLdeathsdep + ALLdeathsndep
  totalYLL <- YLLdep + YLLndep
  
  prop =round(c(colSums(SADdep)["2023"]/colSums(ALLdeathsdep)["2023"]*100, colSums(SADdep)["2100"]/colSums(ALLdeathsdep)["2100"]*100,
                colSums(SADndep)["2023"]/colSums(ALLdeathsndep)["2023"]*100, colSums(SADndep)["2100"]/colSums(ALLdeathsndep)["2100"]*100,
                colSums(totalSAD)["2023"]/colSums(totaldeaths)["2023"]*100,colSums(totalSAD)["2100"]/colSums(totaldeaths)["2100"]*100),1)
  annSAD = round(c(colSums(SADdep)["2023"],colSums(SADdep)["2100"],colSums(SADndep)["2023"],colSums(SADndep)["2100"],
                   colSums(totalSAD)["2023"], colSums(totalSAD)["2100"]),-3)  
  
  cSAD = round(c(colSums(SADdep)["2023"], sum(colSums(SADdep)[124:201]),colSums(SADndep)["2023"], sum(colSums(SADndep)[124:201]),
                 colSums(totalSAD)["2023"], sum(colSums(totalSAD)[124:201])),-3)
  
  annYLL = round(c(colSums(YLLdep)["2023"],colSums(YLLdep)["2100"],colSums(YLLndep)["2023"],colSums(YLLndep)["2100"],
                   colSums(totalYLL)["2023"], colSums(totalYLL)["2100"]),-3)  
  cYLL = round(c(colSums(YLLdep)["2023"], sum(colSums(YLLdep)[1:78]),colSums(YLLndep)["2023"], sum(colSums(YLLndep)[1:78]),
                 colSums(totalYLL)["2023"], sum(colSums(totalYLL)[1:78])),-3)
  
  smkprev = c( cs_deppop["2023"][6,1]*100, cs_deppop["2100"][6,1]*100, cs_ndeppop["2023"][6,1]*100 , cs_ndeppop["2100"][6,1]*100,  cs_totalpop["2023"][6,1]*100, cs_totalpop["2100"][6,1]*100)
  prevchange = c(NA, (cs_deppop["2023"][6,1] - cs_deppop["2100"][6,1])/cs_deppop["2023"][6,1]*100, NA,  
                 (cs_ndeppop["2023"][6,1] - cs_ndeppop["2100"][6,1])/cs_ndeppop["2023"][6,1]*100, NA, 
                 (cs_totalpop["2023"][6,1] - cs_totalpop["2100"][6,1])/cs_totalpop["2023"][6,1]*100)
  PR = c( cs_deppop["2023"][6,1]/cs_ndeppop["2023"][6,1], cs_deppop["2100"][6,1]/cs_ndeppop["2100"][6,1], NA, NA, NA, NA )
  
  Table1 <- rbind(sprintf("%.1f",smkprev), sprintf("%.1f",prop),format(annSAD,big.mark=","),format(cSAD,big.mark=","),format(annYLL,big.mark=","),format(cYLL,big.mark=","), sprintf("%.2f",PR),sprintf("%.1f",prevchange) )
  
  colnames(Table1) = c(paste0("dep2023"),paste0("dep2100"), paste0("ndep2023"), paste0("ndep2100"), paste0("total2023"), paste0("total2100"))
  
  return(Table1)
}

baselineF = cbind(scenarios[1],"Females", scenariotable(SADdepF0, SADndepF0, YLLdepF0, YLLndepF0, ALLdeathsdepF0, ALLdeathsndepF0,cs_deppopF0,cs_ndeppopF0,cs_totalpopF0,"F"))
baselineM = cbind(scenarios[1],"Males",scenariotable(SADdepM0, SADndepM0, YLLdepM0, YLLndepM0, ALLdeathsdepM0, ALLdeathsndepM0,cs_deppopM0,cs_ndeppopM0,cs_totalpopM0,"M"))
baseline = rbind(baselineM, baselineF)

equalityF = cbind(scenarios[2],"Females", scenariotable(SADdepF1, SADndepF1, YLLdepF1, YLLndepF1, ALLdeathsdepF1, ALLdeathsndepF1,cs_deppopF1,cs_ndeppopF1,cs_totalpopF1,"F"))
equalityM = cbind(scenarios[2],"Males", scenariotable(SADdepM1, SADndepM1, YLLdepM1, YLLndepM1, ALLdeathsdepM1, ALLdeathsndepM1,cs_deppopM1,cs_ndeppopM1,cs_totalpopM1,"M"))
equality = rbind(equalityM, equalityF)

equityF = cbind(scenarios[3],"Females", scenariotable(SADdepF2, SADndepF2, YLLdepF2, YLLndepF2, ALLdeathsdepF2, ALLdeathsndepF2,cs_deppopF2,cs_ndeppopF2,cs_totalpopF2,"F"))
equityM = cbind(scenarios[3],"Males", scenariotable(SADdepM2, SADndepM2, YLLdepM2, YLLndepM2, ALLdeathsdepM2, ALLdeathsndepM2,cs_deppopM2,cs_ndeppopM2,cs_totalpopM2,"M"))
equity = rbind(equityM, equityF)

data <- rbind(baseline,equality,equity)
rownames(data)= rep(c("smkprevM","propM","annSADM","cSADM", "annYLLM","cYLLM","PRM","prevchangeM","smkprevF","propF","annSADF","cSADF", "annYLLF","cYLLF","PRF","prevchangeF"),3)

write.csv(data,file=paste0("Tables",namethisrun,".csv"))
library(rio)

write.xlsx(baseline, file = paste0("Tables",namethisrun,".xlsx"),
           sheetName = "statusquo", append = TRUE)
export(baseline, "clipboard",  col.names = TRUE, row.names = TRUE)

write.xlsx(equality, file = paste0("Tables",namethisrun,".xlsx"),
           sheetName="equality scenario", append=TRUE)
export(equality, "clipboard",  col.names = TRUE, row.names = TRUE)

write.xlsx(equity, file = paste0("Tables",namethisrun,".xlsx"),
           sheetName="equity scenario", append=TRUE)
export(equity, "clipboard",  col.names = TRUE, row.names = TRUE)


## tables for figures
male_spd_model = as.data.frame(t(cs_deppopM0[106:201][6,]))
male_spd_model_equality = as.data.frame(t(cs_deppopM1[106:201][6,]))
male_spd_model_equity = as.data.frame(t(cs_deppopM2[106:201][6,]))

male_nospd_model = as.data.frame(t(cs_ndeppopM0[106:201][6,]))
male_nospd_model_equality = as.data.frame(t(cs_ndeppopM1[106:201][6,]))
male_nospd_model_equity = as.data.frame(t(cs_ndeppopM2[106:201][6,]))

female_spd_model = as.data.frame(t(cs_deppopF0[106:201][6,]))
female_spd_model_equality = as.data.frame(t(cs_deppopF1[106:201][6,]))
female_spd_model_equity = as.data.frame(t(cs_deppopF2[106:201][6,]))

female_nospd_model = as.data.frame(t(cs_ndeppopF0[106:201][6,]))
female_nospd_model_equality = as.data.frame(t(cs_ndeppopF1[106:201][6,]))
female_nospd_model_equity = as.data.frame(t(cs_ndeppopF2[106:201][6,]))

write.csv(male_spd_model,"male_spd_model.csv")
write.csv(male_spd_model_equality,"male_spd_model_equality.csv")
write.csv(male_spd_model_equity,"male_spd_model_equity.csv")

write.csv(male_nospd_model,"male_nospd_model.csv")
write.csv(male_nospd_model_equality,"male_nospd_model_equality.csv")
write.csv(male_nospd_model_equity,"male_nospd_model_equity.csv")

write.csv(female_spd_model,"female_spd_model.csv")
write.csv(female_spd_model_equality,"female_spd_model_equality.csv")
write.csv(female_spd_model_equity,"female_spd_model_equity.csv")

write.csv(female_nospd_model,"female_nospd_model.csv")
write.csv(female_nospd_model_equality,"female_nospd_model_equality.csv")
write.csv(female_nospd_model_equity,"female_nospd_model_equity.csv")

library(rio)
female_nospd <- matrix(data = NA,nrow = 104,ncol = 3)
female_nospd[1:104,1] <- female_nospd_model[1:104,1]
female_nospd[1:104,2] <- female_nospd_model_equality[1:104,1]
female_nospd[1:104,3] <- female_nospd_model_equity[1:104,1]
export(female_nospd, "clipboard",col.names = FALSE)

female_spd <- matrix(data = NA,nrow = 104,ncol = 3)
female_spd[1:104,1] <- female_spd_model[1:104,1]
female_spd[1:104,2] <- female_spd_model_equality[1:104,1]
female_spd[1:104,3] <- female_spd_model_equity[1:104,1]
export(female_spd, "clipboard",col.names = FALSE)

male_nospd <- matrix(data = NA,nrow = 104,ncol = 3)
male_nospd[1:104,1] <- male_nospd_model[1:104,1]
male_nospd[1:104,2] <- male_nospd_model_equality[1:104,1]
male_nospd[1:104,3] <- male_nospd_model_equity[1:104,1]
export(male_nospd, "clipboard",col.names = FALSE)

male_spd <- matrix(data = NA,nrow = 104,ncol = 3)
male_spd[1:104,1] <- male_spd_model[1:104,1]
male_spd[1:104,2] <- male_spd_model_equality[1:104,1]
male_spd[1:104,3] <- male_spd_model_equity[1:104,1]
export(male_spd, "clipboard",col.names = FALSE)

write.csv(depsmkprevs_by_year,"depsmkprevs_by_year.csv")
