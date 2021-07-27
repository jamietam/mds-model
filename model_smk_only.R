rm(list = ls())
mainDir <- "C:/Users/jamietam/Dropbox/GitHub/smk-dep-model"
setwd(file.path(mainDir))
library(openxlsx)
library(reshape)
library(Hmisc)
load("C:/Users/jamietam/Dropbox/GitHub/smk-dep-model/depsmkprevs_2005-2018_v2.rda")

date = "smkonly_newcisnetrates" # name the folder where results will be saved
namethisrun = "nhis2018" # name this run

startyear = 1900 # the burn-in period starting point
endyear = 2018 # census projections through 2060
startage = 0
endage = 99
Ny= endyear - startyear + 1 
Na= endage - startage + 1
emptycompartment <- matrix(0, nrow = Na, ncol = Ny, dimnames=list(c(startage:endage),c(startyear:endyear))) # create matrix of zeroes for compartments

allparamsF = read.xlsx("parameters_tx.xlsx",sheet=paste0("smk_females"),rowNames=TRUE,colNames=TRUE) # Adjust parameters in this excel file with estimates to be used in the model
paramsF = as.vector(subset(allparamsF,bhat==1)[['estimate']]) # Parameters where bhat = 1 can be estimated by bhat

paramsnamesF = rownames(subset(allparamsF,bhat==1))  
lowervectorF = as.vector(subset(allparamsF,bhat==1)[['lower']]) # lower bounds for bhat parameters
uppervectorF = as.vector(subset(allparamsF,bhat==1)[['upper']]) # upper bounds for bhat parameters

allparamsM = read.xlsx("parameters_tx.xlsx",sheet=paste0("smk_males"),rowNames=TRUE,colNames=TRUE) # Adjust parameters in this excel file with estimates to be used in the model
paramsM = as.vector(subset(allparamsM,bhat==1)[['estimate']]) # Parameters where bhat = 1 can be estimated by bhat

paramsnamesM = rownames(subset(allparamsM,bhat==1))  
lowervectorM = as.vector(subset(allparamsM,bhat==1)[['lower']]) # lower bounds for bhat parameters
uppervectorM = as.vector(subset(allparamsM,bhat==1)[['upper']]) # upper bounds for bhat parameters


# Get model prevs by age group --------------------------------------------
getmodelprevs <- function(numerator,denominator){
  numerator = as.data.frame(numerator)
  denominator = as.data.frame(denominator)
  agegroupstart <- c(18,26,35,50,65,18)
  agegroupend <- c(25,34,49,64,99,99)
  prevs = NULL
  for (a in c(1:length(agegroupstart))) {
    prevs <- rbind(prevs,colSums(numerator[c((agegroupstart[a]+1):(agegroupend[a]+1)), ],na.rm=TRUE)/
                     colSums(denominator[c((agegroupstart[a]+1):(agegroupend[a]+1)), ],na.rm=TRUE))
  }
  row.names(prevs)<-c("18to25", "26to34", "35to49", "50to64",  "65plus", "total")
  return(prevs)
}

# Main model --------------------------------------------------------------
main <- function(getmodelprevs,params,whichgender,allparams,paramsnames){
  
  pop = read.xlsx("census_data/np2017_d1.xlsx",sheet=whichgender,rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  smk_init_cisnet = read.xlsx("cisnet_smkrates_nhis2018.xlsx",sheet=paste0(whichgender,"_init"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  smk_cess_cisnet = read.xlsx("cisnet_smkrates_nhis2018.xlsx",sheet=paste0(whichgender,"_cess"),rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  death_ns = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("ns_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  death_cs = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("cs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  death_fs = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("fs_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 

  # model parameters are retrieved from the parameters.xlsx file, if they are bhat=1, then they are read in as a params vector for bhat estimation
  smkinit_SF = matrix(0,100,1)
  smkinit_SF[1:18] = ifelse(allparams["smkinit_youthSF","bhat"]==0,allparams["smkinit_youthSF","estimate"],params[match("smkinit_youthSF",paramsnames)])
  smkinit_SF[19:35] = ifelse(allparams["smkinit_SF_18to34","bhat"]==0,allparams["smkinit_SF_18to34","estimate"],params[match("smkinit_SF_18to34",paramsnames)])
  smkinit_SF[36:66] = ifelse(allparams["smkinit_SF_35to64","bhat"]==0,allparams["smkinit_SF_35to64","estimate"],params[match("smkinit_SF_35to64",paramsnames)])
  smkinit_SF[67:100] = ifelse(allparams["smkinit_SF_65plus","bhat"]==0,allparams["smkinit_SF_65plus","estimate"],params[match("smkinit_SF_65plus",paramsnames)])
  
  smkcess_SF = matrix(0,100,1)
  smkcess_SF[1:18] = ifelse(allparams["smkcess_youthSF","bhat"]==0,allparams["smkcess_youthSF","estimate"],params[match("smkcess_youthSF",paramsnames)])
  smkcess_SF[19:35] = ifelse(allparams["smkcess_SF_18to34","bhat"]==0,allparams["smkcess_SF_18to34","estimate"],params[match("smkcess_SF_18to34",paramsnames)])
  smkcess_SF[36:66] = ifelse(allparams["smkcess_SF_35to64","bhat"]==0,allparams["smkcess_SF_35to64","estimate"],params[match("smkcess_SF_35to64",paramsnames)])
  smkcess_SF[67:100] = ifelse(allparams["smkcess_SF_65plus","bhat"]==0,allparams["smkcess_SF_65plus","estimate"],params[match("smkcess_SF_65plus",paramsnames)])
  
  smk_init = smk_init_cisnet*smkinit_SF # scale smoking initiation rates ### CHECK DIMENSIONS might be off by one
  smk_cess = smk_cess_cisnet*smkcess_SF # scale smoking cessation rates
  
  # Compartments / state variables ------------------------------------------
  
  # Initialize population - Model compartments are organized with age 0-99 as rows, and year as columns
  matrix.names<-c('ns', 'cs', 'fs')
  for (name in matrix.names) assign(name,emptycompartment)

  ns[paste(startage),1:Ny] <- as.matrix(pop[paste(startage),1:Ny])   # Takes empty compartment and populates the top row of the matrix with the number of 0-yrolds
  
  for (y in c((startyear+1):endyear)){
    py = paste(y - 1)
    ns[2:Na,paste(y)] <- ns[1:Na-1,py]*(1-smk_init[(startage+1):endage,py])*(1-death_ns[(startage+1):endage,py])
    cs[2:Na,paste(y)] <- cs[1:Na-1,py]*(1-smk_cess[(startage+1):endage,py])*(1-death_cs[(startage+1):endage,py]) + (smk_init[(startage+1):endage,py])*ns[1:Na-1,py]
    fs[2:Na,paste(y)] <- fs[1:Na-1,py]*(1-death_fs[(startage+1):endage,py]) + (smk_cess[(startage+1):endage,py])*cs[1:Na-1,py]
      
  }
  
  # Get population counts for each subpopulation ----------------------------
  totalpop = ns+cs+fs
  
  s1 <- getmodelprevs(ns,totalpop) # never smoker prevalence among the entire never depressed population
  s2 <- getmodelprevs(cs,totalpop)
  s3 <- getmodelprevs(fs,totalpop)
  
  # modeldepprevdata <- list(d1,d2,d3,e1,h1,h2,d4,d5,d6)
  initiation <- data.frame(cbind(smk_init_cisnet["2018"],smkinit_SF,smk_init["2018"],c(0:99)))
  cessation <- data.frame(cbind(smk_cess_cisnet["2018"],smkcess_SF,smk_cess["2018"],c(0:99)))
  
  colnames(initiation) <- c("cisnet","SF","scaledrates","age")
  colnames(cessation) <- c("cisnet","SF","scaledrates","age")
  
  
  modelsmkprevdata <- list(s1,s2,s3,initiation,cessation)
  
  return(modelsmkprevdata)
}

# Get NSDUH prevs by age group --------------------------------------------
getnsduhprevs <- function(depsmkprevs_by_year,assignedsex,numpop,denompop){
  nsduhdata = NULL
  nsduhdata <- melt(depsmkprevs_by_year,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev") )
  nsduhdata <- cast(nsduhdata, group~survey_year,mean,subset=gender==assignedsex&subpopulation==denompop&status==numpop)
  nsduhdata <- nsduhdata[c(-1)] # remove group name column
  row.names(nsduhdata)<-c("18to25", "26to34", "35to49", "50to64",  "65plus", "total")
  return(nsduhdata)
}

# Get model sum of squares ------------------------------------------------
getsumdiffs = function(out,whichgender){
  years<- c("X2005","X2006","X2007","X2008","X2009","X2010","X2011","X2012","X2013","X2014","X2015","X2016","X2017") # only look at output for years where NSDUH data are available
  years2 <- c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017") 
  smkdiffs <- rowSums(subset(as.data.frame(out[1]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,whichgender, "neversmoker","totalpop"),select=years2))^2 + # sum of squares
    rowSums(subset(as.data.frame(out[2]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,whichgender ,"currentsmoker","totalpop"),select=years2))^2 +
    rowSums(subset(as.data.frame(out[3]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,whichgender, "formersmoker","totalpop"),select=years2))^2
  return(sum(smkdiffs))
}


# Bhat parameter estimation -----------------------------------------------

library(Bhat)
xF <- list(label=paramsnamesF, est=paramsF,low=lowervectorF,upp=uppervectorF) # est = parameter starting values
xM <- list(label=paramsnamesM, est=paramsM,low=lowervectorM,upp=uppervectorM)

ML_bhatF=function(paramsF){
  outF = main(getmodelprevs, paramsF, "females",allparamsF,paramsnamesF)
  LL = sum(getsumdiffs(outF,"females"))   #Least squares
  cat(LL,paramsF,'\n')
  return(LL)
}

resbhatF=dfp(xF,ML_bhatF)
paramsF=resbhatF$est  # replace params with the bhat parameter estimates

ML_bhatM=function(paramsM){
  outM = main(getmodelprevs, paramsM, "males",allparamsM,paramsnamesM)
  LL = sum(getsumdiffs(outM,"males"))   #Least squares
  cat(LL,paramsM,'\n')
  return(LL)
}

resbhatM=dfp(xM,ML_bhatM)
paramsM=resbhatM$est  # replace params with the bhat parameter estimates

# Run the model with parameter estimates
outF = main(getmodelprevs, paramsF, "females",allparamsF,paramsnamesF)
outM = main(getmodelprevs, paramsM, "males",allparamsM,paramsnamesM) # runs model using parameters specified in excel sheet OR using bhat estimates

# Data visualization and results figures -------------------------------------------------------------------

library(ggplot2)
library(reshape)
library(grid)
library(gridBase)
library(gridExtra)
library(plyr)

xaxisbreaks = seq(2005,2018,1) # specify the ticks on the x-axis of your results plots

dir.create(file.path(mainDir, date), showWarnings = FALSE)

theme_set( theme_light(base_size = 20))

namethisfile = paste0(namethisrun,".pdf")

s1F <- data.frame(outF[1],check.names = FALSE) # ns,totalpop
s2F <- data.frame(outF[2],check.names = FALSE) # cs,totalpop
s3F <- data.frame(outF[3],check.names = FALSE) # fs,totalpop
initiationF = data.frame(outF[4])
cessationF = data.frame(outF[5])
s1M <- data.frame(outM[1],check.names = FALSE) # ns,totalpop
s2M <- data.frame(outM[2],check.names = FALSE) # cs,totalpop
s3M <- data.frame(outM[3],check.names = FALSE) # fs,totalpop
initiationM = data.frame(outM[4])
cessationM = data.frame(outM[5])

# Combine plots with shared legend
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

getnsduhprevsCI <- function(depsmkprevs_by_year,assignedsex,numpop,denompop){
  nsduhdatalow = NULL
  nsduhdatalow <- melt(depsmkprevs_by_year,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev_lowCI") )
  nsduhdatalow <- cast(nsduhdatalow, group~survey_year,mean,subset=gender==assignedsex&subpopulation==denompop&status==numpop)
  nsduhdatalow <- nsduhdatalow[c(-1)] # remove group name column
  row.names(nsduhdatalow)<-c("18to25", "26to34", "35to49", "50to64",  "65plus", "total")
  nsduhdatahigh = NULL
  nsduhdatahigh <- melt(depsmkprevs_by_year,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev_highCI") )
  nsduhdatahigh <- cast(nsduhdatahigh, group~survey_year,mean,subset=gender==assignedsex&subpopulation==denompop&status==numpop)
  nsduhdatahigh <- nsduhdatahigh[c(-1)] # remove group name column
  row.names(nsduhdatahigh)<-c("18to25", "26to34", "35to49", "50to64",  "65plus", "total")
  return(list(nsduhdatalow, nsduhdatahigh))
}

# Compare model data with NSDUH data 2005-2018 ----------------------------

comparison <- function(s1F,s2F,s3F,prev1, prev2, prev3, population,agegroup,getnsduhprevs,assignedsex){
  s1F$status <- prev1
  s2F$status <- prev2
  s3F$status <- prev3
  s1F$age<-rownames(s1F)
  s2F$age<-rownames(s2F)
  s3F$age<-rownames(s3F)
  model <- rbind(s1F,s2F)
  model <- rbind(model,s3F)
  model$grp <- paste0(model$status, "_model")
  model <-melt(as.data.frame(model),id.vars=c("grp","age","status"))
  colnames(model)[5] <- "modelvalue"
  s1Fnsduh <- getnsduhprevs(depsmkprevs_by_year,assignedsex,prev1,population)
  s1Fnsduhlow <- getnsduhprevsCI(depsmkprevs_by_year,assignedsex,prev1,population)[[1]]
  s1Fnsduhhigh <- getnsduhprevsCI(depsmkprevs_by_year,assignedsex,prev1,population)[[2]]
  s1Fnsduh$status <- prev1
  s1Fnsduhlow$status <- prev1
  s1Fnsduhhigh$status <- prev1
  s1Fnsduh$age<-rownames(s1Fnsduh)
  s1Fnsduhlow$age<-rownames(s1Fnsduhlow)
  s1Fnsduhhigh$age<-rownames(s1Fnsduhhigh)

  s2Fnsduh <- getnsduhprevs(depsmkprevs_by_year,assignedsex,prev2,population)
  s2Fnsduhlow <- getnsduhprevsCI(depsmkprevs_by_year,assignedsex,prev2,population)[[1]]
  s2Fnsduhhigh <- getnsduhprevsCI(depsmkprevs_by_year,assignedsex,prev2,population)[[2]]
  s2Fnsduh$status <- prev2
  s2Fnsduhlow$status <- prev2
  s2Fnsduhhigh$status <- prev2
  s2Fnsduh$age<-rownames(s2Fnsduh)
  s2Fnsduhlow$age<-rownames(s2Fnsduhlow)
  s2Fnsduhhigh$age<-rownames(s2Fnsduhhigh)

  s3Fnsduh <- getnsduhprevs(depsmkprevs_by_year,assignedsex,prev3,population)
  s3Fnsduhlow <- getnsduhprevsCI(depsmkprevs_by_year,assignedsex,prev3,population)[[1]]
  s3Fnsduhhigh <- getnsduhprevsCI(depsmkprevs_by_year,assignedsex,prev3,population)[[2]]
  s3Fnsduh$status <- prev3
  s3Fnsduhlow$status <- prev3
  s3Fnsduhhigh$status <- prev3
  s3Fnsduh$age<-rownames(s3Fnsduh)
  s3Fnsduhlow$age<-rownames(s3Fnsduhlow)
  s3Fnsduhhigh$age<-rownames(s3Fnsduhhigh)

  nsduh <- rbind(s1Fnsduh,s2Fnsduh)
  nsduh <- rbind(nsduh,s3Fnsduh)
  nsduh$grp <- paste0(nsduh$status, "_nsduh")

  nsduhlow <- rbind(s1Fnsduhlow,s2Fnsduhlow)
  nsduhlow <- rbind(nsduhlow,s3Fnsduhlow)
  nsduhlow$grp <- paste0(nsduhlow$status, "_nsduhlow")

  nsduhhigh <- rbind(s1Fnsduhhigh,s2Fnsduhhigh)
  nsduhhigh <- rbind(nsduhhigh,s3Fnsduhhigh)
  nsduhhigh$grp <- paste0(nsduhhigh$status, "_nsduhhigh")

  nsduh<-melt(as.data.frame(nsduh),id.vars=c("grp","age","status")) ## Add CIs to this dataframe
  nsduhlow<-melt(as.data.frame(nsduhlow),id.vars=c("grp","age","status")) ## Add CIs to this dataframe
  nsduhhigh<-melt(as.data.frame(nsduhhigh),id.vars=c("grp","age","status")) ## Add CIs to this dataframe

  colnames(nsduh)[5] <- "nsduhvalue"
  colnames(nsduhlow)[5] <- "nsduhvalue_low"
  colnames(nsduhhigh)[5] <- "nsduhvalue_high"

  nsduh$nsduh_low <- nsduhlow$nsduhvalue_low
  nsduh$nsduh_high <- nsduhhigh$nsduhvalue_high

  comparedata <- merge(nsduh, model, by=c("age","status","variable"),all.y=TRUE)
  comparedata$variable <- as.numeric(levels(comparedata$variable))[comparedata$variable] # converts factor to numeric
  g <- ggplot(subset(comparedata,age==agegroup&variable>=2005)) +
    geom_line(aes(x=variable, y= modelvalue,colour=grp.y ))+
    geom_pointrange(aes(x = variable, y = nsduhvalue, ymin=nsduh_low, ymax = nsduh_high, colour=grp.x))  +
    labs(title=paste0(assignedsex)) +
    scale_y_continuous(name="prevalence",limits=c(0,1),breaks=seq(0,1,0.1)) +
    # scale_x_continuous(name="year",limits=c(2000,endyear),breaks=seq(2000,endyear,5))  +
    scale_x_continuous(name="year",limits=c(2005,endyear),breaks=xaxisbreaks)  +
    theme(axis.text.x=element_text(angle=60, hjust=1))+

    xlab("year")+
    scale_colour_manual(name="",values = c("red", "red", "black", "black", "blue", "blue"),
                        guide = guide_legend(override.aes = list(
                          linetype = c("solid", "blank","solid","blank","solid","blank"),
                          shape = c(NA,16,NA,16,NA,16))))#+
  return(g)
}

tF <- comparison(s1F,s2F,s3F,"neversmoker", "currentsmoker", "formersmoker", "totalpop","total",getnsduhprevs,"females")
tM <- comparison(s1M,s2M,s3M,"neversmoker", "currentsmoker", "formersmoker", "totalpop","total",getnsduhprevs,"males")

# Combine males and females into one figure
s2M$Sex <- "Males"
s2F$Sex <- "Females"
s2M$age<-rownames(s2M)
s2F$age<-rownames(s2F)
model <- rbind(s2M,s2F)
model <-melt(as.data.frame(model),id.vars=c("Sex","age"))
colnames(model)[5] <- "Model"
s2Mnsduh <- getnsduhprevs(depsmkprevs_by_year,"males","currentsmoker","totalpop")
s2Mnsduhlow <- getnsduhprevsCI(depsmkprevs_by_year,"males","currentsmoker","totalpop")[[1]]
s2Mnsduhhigh <- getnsduhprevsCI(depsmkprevs_by_year,"males","currentsmoker","totalpop")[[2]]
s2Mnsduh$Sex <- "Males"
s2Mnsduhlow$Sex <- "Males"
s2Mnsduhhigh$Sex <- "Males"
s2Mnsduh$age<-rownames(s2Mnsduh)
s2Mnsduhlow$age<-rownames(s2Mnsduhlow)
s2Mnsduhhigh$age<-rownames(s2Mnsduhhigh)

s2Fnsduh <- getnsduhprevs(depsmkprevs_by_year,"females","currentsmoker","totalpop")
s2Fnsduhlow <- getnsduhprevsCI(depsmkprevs_by_year,"females","currentsmoker","totalpop")[[1]]
s2Fnsduhhigh <- getnsduhprevsCI(depsmkprevs_by_year,"females","currentsmoker","totalpop")[[2]]
s2Fnsduh$Sex <- "Females"
s2Fnsduhlow$Sex <- "Females"
s2Fnsduhhigh$Sex <- "Females"
s2Fnsduh$age<-rownames(s2Fnsduh)
s2Fnsduhlow$age<-rownames(s2Fnsduhlow)
s2Fnsduhhigh$age<-rownames(s2Fnsduhhigh)

nsduh <- rbind(s2Fnsduh,s2Mnsduh)
nsduh$grp <- "NSDUH"

nsduhlow <- rbind(s2Fnsduhlow,s2Mnsduhlow)
nsduhlow$grp <- "NSDUH_low"

nsduhhigh <- rbind(s2Fnsduhhigh,s2Mnsduhhigh)
nsduhhigh$grp <- "NSDUH_high"

nsduh<-melt(as.data.frame(nsduh),id.vars=c("grp","age","Sex")) ## Add CIs to this dataframe
nsduhlow<-melt(as.data.frame(nsduhlow),id.vars=c("grp","age","Sex")) ## Add CIs to this dataframe
nsduhhigh<-melt(as.data.frame(nsduhhigh),id.vars=c("grp","age","Sex")) ## Add CIs to this dataframe

colnames(nsduh)[5] <- "nsduhvalue"
colnames(nsduhlow)[5] <- "nsduhvalue_low"
colnames(nsduhhigh)[5] <- "nsduhvalue_high"

nsduh$nsduh_low <- nsduhlow$nsduhvalue_low
nsduh$nsduh_high <- nsduhhigh$nsduhvalue_high

comparedata <- merge(nsduh, model, by=c("age","Sex","variable"),all.y=TRUE)
comparedata$variable <- as.numeric(levels(comparedata$variable))[comparedata$variable] # converts factor to numeric
g <- ggplot(subset(comparedata,age=="total"&variable>=2005)) +
  geom_line(aes(x=variable, y= value*100,colour=Sex,linetype=Sex ))+
  geom_pointrange(aes(x = variable, y = nsduhvalue*100, ymin=nsduh_low*100, ymax = nsduh_high*100, colour=Sex))  +
  scale_y_continuous(name="Prevalence (%)",limits=c(15,35),breaks=seq(15,35,5)) +
  scale_x_continuous(name="Year",limits=c(2005,endyear),breaks=xaxisbreaks)  +
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.title = element_blank())


jpeg(filename = paste0("FigS4_smkonly_", namethisrun,".jpg"),width=6, height=6, units ="in", res=1000)
g
dev.off()

# Generate plots of prevalence by age over time ----------------------------------

plotprevbyage <- function(df,status,nsduhstatus,subpopulation,whichgender){
  thissubset <- depsmkprevs_by_year[depsmkprevs_by_year$status==nsduhstatus & depsmkprevs_by_year$gender==whichgender & depsmkprevs_by_year$subpopulation==subpopulation & depsmkprevs_by_year$age!="total",]
  nsduhdata <- melt(thissubset, id.vars = c("age","survey_year"),measure.vars = c("prev"))
  nsduhdatalow <- melt(thissubset, id.vars = c("age","survey_year"),measure.vars = c("prev_lowCI"))
  nsduhdatahigh <- melt(thissubset, id.vars = c("age","survey_year"),measure.vars = c("prev_highCI"))
  colnames(nsduhdata)[4] <- "nsduh_prev"
  nsduhdata$nsduh_low <- nsduhdatalow$value
  nsduhdata$nsduh_high <- nsduhdatahigh$value

  df$age<-rownames(df)
  df <- melt(subset(df, df$age!="total"), id.vars=c("age"))
  colnames(df)[2] <- "survey_year"
  colnames(df)[3] <- "model_prev"
  comparedata <- merge(nsduhdata, df, by=c("age","survey_year"))
  comparedata$survey_year <- as.numeric(comparedata$survey_year)
  comparedata$model_prev <- as.numeric(comparedata$model_prev)
  g <- ggplot(comparedata, aes(x=survey_year, y=model_prev, group=age)) +
    geom_line(aes(x=survey_year, y= model_prev,colour=age ))+
    geom_pointrange(aes(x = survey_year, y = nsduh_prev, ymin=nsduh_low, ymax = nsduh_high, colour=age))  +
    scale_y_continuous(name="prevalence",limits=c(0,1),breaks=seq(0,1,0.05)) +
    scale_x_continuous(name="Year",limits=c(2005,endyear),breaks=xaxisbreaks)  +
    theme(axis.text.x=element_text(angle=45, hjust=1)) #+
    labs(title=paste0(status),color="age")# +
  return(g)
}

s1plotF <- plotprevbyage(s1F, "Never smoker", "neversmoker", "totalpop","females")
s2plotF <- plotprevbyage(s2F, "Current smoker", "currentsmoker", "totalpop","females")
s3plotF <- plotprevbyage(s3F, "Former smoker", "formersmoker", "totalpop","females")
s1plotM <- plotprevbyage(s1M, "Never smoker", "neversmoker", "totalpop","males")
s2plotM <- plotprevbyage(s2M, "Current smoker", "currentsmoker", "totalpop","males")
s3plotM <- plotprevbyage(s3M, "Former smoker", "formersmoker", "totalpop","males")


# Plot initiation and cessation probabilities -----------------------------

initplotM <- ggplot(initiationM, aes(age,cisnet)) + 
  geom_point(aes(color="CISNET")) +
  geom_line(aes(x=age, y= scaledrates,color="model" ))+
  scale_x_continuous(name="Age",breaks=c(0,18,35,65,100)) +
  scale_y_continuous(name="Annual probability",limits = c(0,0.10),breaks=seq(0,0.1,0.01))+
  labs(title=paste0("Males"))+
  theme(legend.position = c(1, 1),legend.justification = c(1, 1))+
  scale_colour_manual(name="",values = c("black", "black"), guide = guide_legend(override.aes = list(
    linetype = c("blank", "solid"), shape = c(16,NA))))


initplotF <- ggplot(initiationF, aes(age,cisnet)) + 
  geom_point(aes(color="CISNET")) +
  geom_line(aes(x=age, y= scaledrates,color="model" ))+
  scale_x_continuous(name="Age",breaks=c(0,18,35,65,100)) +
  scale_y_continuous(name="Annual probability",limits = c(0,0.10),breaks=seq(0,0.1,0.01))+
  labs(title=paste0("Females"))+
  theme(legend.position = c(1, 1),legend.justification = c(1, 1))+
  scale_colour_manual(name="",values = c("black", "black"), guide = guide_legend(override.aes = list(
    linetype = c("blank", "solid"), shape = c(16,NA))))


cessplotF<- ggplot(cessationF, aes(age,cisnet)) + 
  geom_point(aes(color="CISNET")) +
  geom_line(aes(x=age, y= scaledrates,color="model" ))+
  scale_x_continuous(name="Age",breaks=c(0,18,35,65,100)) +
  scale_y_continuous(name="Annual probability",limits = c(0,0.10),breaks=seq(0,0.1,0.01))+
  labs(title=paste0("Females"))+
  theme(legend.position = c(1, 1),legend.justification = c(1, 1))+
  scale_colour_manual(name="",values = c("black", "black"), guide = guide_legend(override.aes = list(
    linetype = c("blank", "solid"), shape = c(16,NA))))

cessplotM<- ggplot(cessationM, aes(age,cisnet)) + 
  geom_point(aes(color="CISNET")) +
  geom_line(aes(x=age, y= scaledrates,color="model" ))+
  scale_x_continuous(name="Age",breaks=c(0,18,35,65,100)) +
  scale_y_continuous(name="Annual probability",limits = c(0,0.10),breaks=seq(0,0.1,0.01))+
  labs(title=paste0("Males"))+
  theme(legend.position = c(1, 1),legend.justification = c(1, 1))+
  scale_colour_manual(name="",values = c("black", "black"), guide = guide_legend(override.aes = list(
    linetype = c("blank", "solid"), shape = c(16,NA))))


jpeg(filename = paste0("scaledinit_", namethisrun,".jpg"),width=10, height=6, units ="in", res=500)
grid_arrange_shared_legend(list(initplotF, initplotM),2,"")
dev.off()

jpeg(filename = paste0("scaledcess_", namethisrun,".jpg"),width=10, height=6, units ="in", res=500)
grid_arrange_shared_legend(list(cessplotF, cessplotM),2,"")
dev.off()

# Create pdf with output figures ---------------------------------
pdf(file = namethisfile, width = 11, height = 8.5)
grid_arrange_shared_legend(list(tF,tM),2,"Smoking Prevalence - total population ages 18+")
grid_arrange_shared_legend(list(s1plotF, s2plotF, s3plotF),3,"Smoking prevalence in FEMALE adult population")
grid_arrange_shared_legend(list(s1plotM, s2plotM, s3plotM),3,"Smoking prevalence in MALE adult population")
grid_arrange_shared_legend(list(initplotF, initplotM),2,"Initiation probabilities")
grid_arrange_shared_legend(list(cessplotF, cessplotM),2,"Cessation probabilities")
dev.off()