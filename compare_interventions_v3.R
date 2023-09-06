library(ggplot2)
library(reshape)
library(grid)
library(gridBase)
library(gridExtra)
library(plyr)
library(dplyr)
library(Hmisc)
library(scales)
library(openxlsx)
library(viridis)

load("C:/Users/JT936/Dropbox/Analysis/NSDUH/mdetx_depsmk_2010-2018.rda")
load("C:/Users/JT936/Dropbox/GitHub/mds-model/depsmkprevs_2005-2018_v2.rda")

death_nsF = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("ns_females"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
death_csF = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("cs_females"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
death_fsF = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("fs_females"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 

LE_csF = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_cs_females"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
LE_fsF = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_fs_females"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 

death_nsM = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("ns_males"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
death_csM = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("cs_males"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
death_fsM = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("fs_males"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 

LE_csM = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_cs_males"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
LE_fsM = read.xlsx("cisnet_deathrates.xlsx",sheet=paste0("LE_fs_males"),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 

xaxisbreaks = c(2018,seq(2020,2100,10)) # specify the ticks on the x-axis of your results plots
minyear = 2018
maxyear = 2100
txeffset = c(1.0, 1.137, 1.588, 2.0, 2.50,3.0)
utilset = c("1.0","1.1","1.2","100")

theme_set( theme_light(base_size = 14))
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
  row.names(nsduhdatalow)<-agerownames
  
  nsduhdatahigh = NULL
  nsduhdatahigh <- melt(depsmkprevs_by_year,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev_highCI") )
  nsduhdatahigh <- cast(nsduhdatahigh, group~survey_year,mean,subset=gender==assignedsex&subpopulation==denompop&status==numpop)
  nsduhdatahigh <- nsduhdatahigh[c(-1)] # remove group name column
  row.names(nsduhdatahigh)<-agerownames
  return(list(nsduhdatalow, nsduhdatahigh))
}

# Figure 2 MH Treatment utilization ------------------------------------------------

for (u in utilset){
  outF = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, mpc=0, txeffcess = 1.0, util=u)
  outM = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, mpc=0, txeffcess = 1.10, util=u)
  
  assign(paste0("tx","M",u), as.data.frame(getmodelprevs(outM[[31]],outM[[23]]))) # cs_dep_tx / cs_dep
  assign(paste0("tx","F",u), as.data.frame(getmodelprevs(outF[[31]],outF[[23]])))
}

mhutil <- function(txM1.0,txM1.1,txM1.2,txM100,whichgender,title){    
  df <- subset(mdetx_depsmk, (mdetx_depsmk$subpopulation=="depsmkpop" & mdetx_depsmk$status=="ahltmde" & mdetx_depsmk$survey_year=="2010-2018" & mdetx_depsmk$gender==tolower(whichgender)))
  df$variable =as.factor("2018")
  df$scenario ="Baseline"
  
  txM1.0$age=row.names(txM1.0)
  txM1.0$scenario="Baseline"
  
  txM1.1$age=row.names(txM1.1)
  txM1.1$scenario="Increase by 10%"
  
  txM1.2$age=row.names(txM1.2)
  txM1.2$scenario="Increase by 20%"
  
  txM100$age=row.names(txM100)
  txM100$scenario="100%"
  
  tx <- rbind(melt(txM1.0[c("2018","age","scenario")]), melt(txM1.1[c("2018","age","scenario")]), melt(txM1.2[c("2018","age","scenario")]), melt(txM100[c("2018","age","scenario")]) )
  tx$scenario  <- factor(tx$scenario , levels = c("Baseline", "Increase by 10%", "Increase by 20%","100%"))
  txdf = merge(tx,df[,c("variable","age","scenario","prev","prev_lowCI","prev_highCI")], by=c("variable","age","scenario"),all=TRUE)
  
  p<- ggplot(subset(txdf,age!="total")) + 
    geom_bar(aes(x=age,y=value*100, fill=scenario, group=scenario),position=position_dodge(), stat="identity")+
    geom_pointrange(aes(x = age, y = prev*100, ymin=prev_lowCI*100, ymax = prev_highCI*100, group=scenario,shape="NSDUH"),position=position_dodge(width=0.9))  +
    scale_fill_viridis(discrete=TRUE, begin=0.75, end = 0, option="viridis" ) +
    scale_y_continuous(name="Percent (%)",limits=c(0,100),breaks=seq(0,100,10))+
    xlab("Age group")+
    labs(title=title) +
    theme(axis.text.x=element_text(angle=30, hjust=1),legend.position="right",legend.title=element_blank()) 
  return(p)
}

uF <- mhutil(txF1.0,txF1.1,txF1.2,txF100,"Females","A) Women with current MD")
uM <- mhutil(txM1.0,txM1.1,txM1.2,txM100,"Males", "B) Men with current MD")

pdf(file = paste0("Fig2_mhutil", namethisrun,".pdf"),width=10, height=6,onefile = FALSE)
grid_arrange_shared_legend(list(uF, uM),2,"")
dev.off()
# plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "",
#      main = paste0(t),
#      cex.lab=0.75, cex.axis=0.75, cex.main=0.75, cex.sub=0.75)
# plot.new()
# plot.window(c(-1,1), c(-1,1))
# 
# s <- "NSDUH = National Survey of Drug Use and Health; MDE = Major depressive episode; Black dots represent the percent of smokers with current MDE who report that they saw a health professional for their MDE within the past year in the NSDUH 2010-2018; vertical black lines are their corresponding 95% confidence intervals. Bars represented calibrated model estimates of mental health service utilization among smokers with MD by scenario: baseline, relative 10% increase in utilization, relative 20% increase in utilization, and 100% complete utilization scenarios. See Supplement for details."
# t <- "Figure 2. Mental health service utilization scenarios among smokers with current MDE"
# rectangleWidth <- 2
# 
# n <- nchar(s)
# for(i in n:1) {
#   wrappeds <- paste0(strwrap(s, i), collapse = "\n")
#   if(strwidth(wrappeds) < rectangleWidth) break
# }
# 
# textHeight <- strheight(wrappeds)       
# text(-1,0, wrappeds, pos=4)


# Figure 3 prevalence projections with NSDUH comparison -------------------------------------------------

getprevsFM <- function(m, t,u){
  outF = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, mpc=m, txeffcess = t, util= u)
  
  cs_depF = outF[[23]]
  deppopF = outF[[13]]
  
  cs_fdepF = outF[[26]]
  fdeppopF = outF[[14]]
  
  cs_deppopF = getmodelprevs(cs_depF,deppopF)
  cs_fdeppopF = getmodelprevs(cs_depF,fdeppopF)
  
  outM = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, mpc=m, txeffcess = t, util= u)
  
  cs_depM = outM[[23]]
  deppopM = outM[[13]]
  
  cs_fdepM = outM[[26]]
  fdeppopM = outM[[14]]
  
  cs_deppopM = getmodelprevs(cs_depM,deppopM)
  cs_fdeppopM = getmodelprevs(cs_depM,fdeppopM)
  
  return(list(cs_deppopF, cs_deppopM, cs_fdeppopF, cs_fdeppopM))
}

scenarios = c("Baseline", "Any Tx","Pharm Tx","100% increase","150% increase","200% increase","Maximum Potential Cessation")

for(t in 1:length(txeffset)){ # Baseline and treatment scenarios
    allprevs = getprevsFM(0,txeffset[t],"1.0")
    assign(paste0("cs_deppop","F",(t-1)), allprevs[[1]])
    assign(paste0("cs_deppop","M",(t-1)), allprevs[[2]])
    assign(paste0("cs_fdeppop","F",(t-1)), allprevs[[3]])
    assign(paste0("cs_fdeppop","M",(t-1)), allprevs[[4]])
}
allprevs = getprevsFM(1,1.0,"1.0") # MPC scenario
assign(paste0("cs_deppop","F",length(scenarios)-1), allprevs[[1]])
assign(paste0("cs_deppop","M",length(scenarios)-1), allprevs[[2]])
assign(paste0("cs_fdeppop","F",length(scenarios)-1), allprevs[[3]])
assign(paste0("cs_fdeppop","M",length(scenarios)-1), allprevs[[4]])

Figmodelnsduh <- function(whichgender,dfs,smkgroup,popgroup, title){
  dfs = lapply(dfs, function(x) {as.data.frame(x)})
  dfs = Map(cbind,dfs,scenario=as.list(scenarios))
  
  fullset = bind_rows(dfs, .id = NULL)
  fullset$age <-rep(agerownames,length(scenarios))

  model <- rbind(melt(fullset,id.vars=c("age","scenario")))
  
  colnames(model)[4] <- "modelvalue"
  model$variable <- as.numeric(levels(model$variable))[model$variable] # converts factor to numeric
  model$scenario  <- factor(model$scenario , levels = scenarios[1:length(scenarios)] )
  
  
  nsduh <- getnsduhprevs(depsmkprevs_by_year,whichgender,smkgroup,popgroup)
  nsduh$age<-agerownames
  nsduh <- melt(as.data.frame(nsduh),id.vars=c("age"))
  colnames(nsduh)[3] <- c("nsduhvalue")
  nsduh$variable <- as.numeric(levels(nsduh$variable))[nsduh$variable] # converts factor to numeric
  
  fig <- ggplot() +
    geom_point(data=subset(nsduh, age=="total"), aes(x = variable, y = nsduhvalue*100, shape="National Survey on Drug Use and Health"))+
    geom_line(data = subset(model,age=="total"), aes(x=variable, y= modelvalue*100, colour=scenario,linetype=scenario, group=rev(scenario)))+
    # scale_color_viridis(discrete=TRUE, begin=0, end = 0.75 , option="magma") +
    scale_color_grey(start=0, end = 0.75) +
    scale_y_continuous(name="Smoking prevalence (%)",limits=c(0,50),breaks=seq(0,50,5)) +
    scale_x_continuous(name="Year",limits=c(2005,2100),breaks=c(2005,seq(2020,2100, 10)))  +
    labs(title=title)+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank())
    # theme(axis.text.x=element_text(angle=60, hjust=1), legend.title = element_blank())
  return(fig)
}

dfsF = list(cs_deppopF0,cs_deppopF1,cs_deppopF2,cs_deppopF3,cs_deppopF4,cs_deppopF5, cs_deppopF6)
Figure3a <- Figmodelnsduh("females",dfsF,"currentsmoker","deppop","A) Women with Major Depressive Episodes")
  
dfsM = list(cs_deppopM0,cs_deppopM1,cs_deppopM2,cs_deppopM3,cs_deppopM4,cs_deppopM5, cs_deppopM6)
Figure3b <- Figmodelnsduh("males",dfsM,"currentsmoker","deppop", "B) Men with Major Depressive Episodes")

pdf(file = paste0("Fig3_smkprev", namethisrun,"_bw_spelledout.pdf"),width=10, height=6, onefile=FALSE)
grid_arrange_shared_legend(list(Figure3a, Figure3b),2,"")
dev.off()

# dfsFfdep = list(cs_fdeppopF0,cs_fdeppopF1,cs_fdeppopF2,cs_fdeppopF3,cs_fdeppopF4,cs_fdeppopF5, cs_fdeppopF6)
# Figmodelnsduh("females",dfsFfdep,"currentsmoker","notdeppop","Former MDE")
# dfsMfdep = list(cs_fdeppopM0,cs_fdeppopM1,cs_fdeppopM2,cs_fdeppopM3,cs_fdeppopM4,cs_fdeppopM5, cs_fdeppopM6)
# Figmodelnsduh("males",dfsMfdep,"currentsmoker","notdeppop","Former MDE")

# Tables - scenario outcomes -------------------------------------------

eachscen <- function(outF0, y1,y2, death_nsF, death_csF, death_fsF){
  cs_everdepF0 = outF0[[20]]
  fs_everdepF0 = outF0[[21]]
  SADF0 = cs_everdepF0[,(y1-1899):(y2-1899)]*(death_csF[,(y1-1899):(y2-1899)]-death_nsF[,(y1-1899):(y2-1899)])+fs_everdepF0[,(y1-1899):(y2-1899)]*(death_fsF[,(y1-1899):(y2-1899)]-death_nsF[,(y1-1899):(y2-1899)])
  cSAD = sum(SADF0,na.rm=TRUE)
  
  cs_depF0 = outF0[[23]]
  fs_depF0 = outF0[[24]]
  depSADF0 = cs_depF0[,(y1-1899):(y2-1899)]*(death_csF[,(y1-1899):(y2-1899)]-death_nsF[,(y1-1899):(y2-1899)])+fs_depF0[,(y1-1899):(y2-1899)]*(death_fsF[,(y1-1899):(y2-1899)]-death_nsF[,(y1-1899):(y2-1899)])
  cSADdep = sum(depSADF0,na.rm=TRUE)
  
  cs_fdepF0 = outF0[[26]]
  fs_fdepF0 = outF0[[27]]
  fdepSADF0 = cs_fdepF0[,(y1-1899):(y2-1899)]*(death_csF[,(y1-1899):(y2-1899)]-death_nsF[,(y1-1899):(y2-1899)])+fs_fdepF0[,(y1-1899):(y2-1899)]*(death_fsF[,(y1-1899):(y2-1899)]-death_nsF[,(y1-1899):(y2-1899)])
  cSADfdep = sum(fdepSADF0,na.rm=TRUE)
  
  cs_recalldepF0 = outF0[[29]]
  fs_recalldepF0 = outF0[[30]]
  recallSADF0 = cs_recalldepF0[,(y1-1899):(y2-1899)]*(death_csF[,(y1-1899):(y2-1899)]-death_nsF[,(y1-1899):(y2-1899)])+fs_recalldepF0[,(y1-1899):(y2-1899)]*(death_fsF[,(y1-1899):(y2-1899)]-death_nsF[,(y1-1899):(y2-1899)])
  cSADrecall = sum(recallSADF0,na.rm=TRUE)
  
  cs_depF0 = outF0[[23]]
  deppopF0 = outF0[[13]]
  prev_dep = round(getmodelprevs(cs_depF0,deppopF0)[6, paste(seq(y1,y2,20))]*100,1)
    
  cs_fdepF0 = outF0[[26]]
  fdeppopF0 = outF0[[14]]
  prev_fdep = round(getmodelprevs(cs_fdepF0,fdeppopF0)[6, paste(seq(y1,y2,20))]*100,1)
  
  cs_recall = outF0[[29]]
  recallpopF0 = outF0[[15]]
  prev_recall = round(getmodelprevs(cs_recall,recallpopF0)[6, paste(seq(y1,y2,20))]*100,1)
  
  return(list(cSAD, SADF0, cSADdep, depSADF0,cSADfdep, fdepSADF0,cSADrecall,recallSADF0,prev_dep, prev_fdep, prev_recall))
}

outF0 = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, mpc = 0, txeffcess=1.0, util="1.0") 
outM0 = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, mpc = 0, txeffcess=1.0, util="1.0")

SADF0 <- eachscen(outF0, 2020,2100, death_nsF, death_csF, death_fsF)[[2]]
SADM0 <- eachscen(outM0, 2020,2100, death_nsM, death_csM, death_fsM)[[2]]

depSADF0 <- eachscen(outF0, 2020,2100, death_nsF, death_csF, death_fsF)[[4]] # current MDE
depSADM0 <- eachscen(outM0, 2020,2100, death_nsM, death_csM, death_fsM)[[4]]

fdepSADF0 <- eachscen(outF0, 2020,2100, death_nsF, death_csF, death_fsF)[[6]] # former MDE
fdepSADM0 <- eachscen(outM0, 2020,2100, death_nsM, death_csM, death_fsM)[[6]]

recallSADF0 <- eachscen(outF0, 2020,2100, death_nsF, death_csF, death_fsF)[[8]] # recall MDE
recallSADM0 <- eachscen(outM0, 2020,2100, death_nsM, death_csM, death_fsM)[[8]]

SADtable = NULL
LYGtable = NULL
depSADtable = NULL
depLYGtable = NULL
fdepSADtable = NULL
fdepLYGtable = NULL
recallSADtable = NULL
recallLYGtable = NULL

depprevtable =  NULL
fdepprevtable =  NULL
recallprevtable =  NULL

for (t in txeffset){
  for (u in utilset){
    for (m in c(0:1)){
      if(m==1&t>1.0) next
      outF = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, mpc=m, txeffcess = t, util=u)
      outM = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, mpc=m, txeffcess = t, util=u)
      
      scenresultsF = eachscen(outF, 2020,2100, death_nsF, death_csF, death_fsF)
      scenresultsM = eachscen(outM, 2020,2100, death_nsM, death_csM, death_fsM)
      
      SADtable = rbind(SADtable,c("females",m,t,u,scenresultsF[[1]]),c("males",m,t,u,scenresultsM[[1]]))
      LYGtable = rbind(LYGtable,
                        c("females",m,t,u,round(sum((LE_fsF[,(2020-1899):(2100-1899)]-LE_csF[,(2020-1899):(2100-1899)])*(SADF0-scenresultsF[[2]]),0))),
                        c("males",m,t,u,round(sum((LE_fsM[,(2020-1899):(2100-1899)]-LE_csM[,(2020-1899):(2100-1899)])*(SADM0-scenresultsM[[2]]),0))))
      
      depSADtable = rbind(depSADtable,c("females",m,t,u,scenresultsF[[3]]),c("males",m,t,u,scenresultsM[[3]]))
      depLYGtable = rbind(depLYGtable,
                          c("females",m,t,u,round(sum((LE_fsF[,(2020-1899):(2100-1899)]-LE_csF[,(2020-1899):(2100-1899)])*(depSADF0-scenresultsF[[4]]),0))),
                          c("males",m,t,u,round(sum((LE_fsM[,(2020-1899):(2100-1899)]-LE_csM[,(2020-1899):(2100-1899)])*(depSADM0-scenresultsM[[4]]),0))))

      fdepSADtable = rbind(fdepSADtable,c("females",m,t,u,scenresultsF[[5]]),c("males",m,t,u,scenresultsM[[5]]))
      fdepLYGtable = rbind(fdepLYGtable, ## SOMETHING IS WRONG WITH THESE LYG NUMBERS
                           c("females",m,t,u,round(sum((LE_fsF[,(2020-1899):(2100-1899)]-LE_csF[,(2020-1899):(2100-1899)])*(fdepSADF0-scenresultsF[[6]]),0))),
                           c("males",m,t,u,round(sum((LE_fsM[,(2020-1899):(2100-1899)]-LE_csM[,(2020-1899):(2100-1899)])*(fdepSADM0-scenresultsM[[6]]),0))))

      recallSADtable = rbind(recallSADtable,c("females",m,t,u,scenresultsF[[7]]),c("males",m,t,u,scenresultsM[[7]]))
      recallLYGtable = rbind(recallLYGtable,
                             c("females",m,t,u,round(sum((LE_fsF[,(2020-1899):(2100-1899)]-LE_csF[,(2020-1899):(2100-1899)])*(recallSADF0-scenresultsF[[8]]),0))),
                             c("males",m,t,u,round(sum((LE_fsM[,(2020-1899):(2100-1899)]-LE_csM[,(2020-1899):(2100-1899)])*(recallSADM0-scenresultsM[[8]]),0))))

      depprevtable = rbind(depprevtable,c("females",m,t,u,scenresultsF[[9]]),c("males",m,t,u,scenresultsM[[9]]))
      fdepprevtable = rbind(fdepprevtable,c("females",m,t,u,scenresultsF[[10]]),c("males",m,t,u,scenresultsM[[10]]))
      recallprevtable = rbind(recallprevtable,c("females",m,t,u,scenresultsF[[11]]),c("males",m,t,u,scenresultsM[[11]]))
    }
  }
}

dfs = list(SADtable,LYGtable,depSADtable,depLYGtable,fdepSADtable,fdepLYGtable,recallSADtable,recallLYGtable,depprevtable,fdepprevtable, recallprevtable)
dfs = lapply(dfs, function(x) {as.data.frame(x)})
colnames(dfs[[1]]) <- c("sex","mpc","txeff","util","SAD averted") #everdep
colnames(dfs[[2]]) <- c("sex","mpc","txeff","util","LYG")
colnames(dfs[[3]]) <- c("sex","mpc","txeff","util","SAD averted") #dep
colnames(dfs[[4]]) <- c("sex","mpc","txeff","util","LYG")
colnames(dfs[[5]]) <- c("sex","mpc","txeff","util","SAD averted") #fdep
colnames(dfs[[6]]) <- c("sex","mpc","txeff","util","LYG")
colnames(dfs[[7]]) <- c("sex","mpc","txeff","util","SAD averted") #recall
colnames(dfs[[8]]) <- c("sex","mpc","txeff","util","LYG")
colnames(dfs[[9]]) <- c("sex","mpc","txeff","util",paste0("depprev",seq(2020,2100,20)))
colnames(dfs[[10]]) <- c("sex","mpc","txeff","util",paste0("fdepprev",seq(2020,2100,20)))
colnames(dfs[[11]]) <- c("sex","mpc","txeff","util",paste0("recallprev",seq(2020,2100,20)))
rewrite <- function(df){ # rewrite the txeff column to include MPC scenarios
  df$txeff[df$mpc==1]="mpc"
  return(df)
}
dfs = lapply(dfs, function(x) {rewrite(x)})

# Table 2
createTable2 <- function(SADdf,LYGdf){
  melted = melt(SADdf, id.vars = c("txeff","util","sex"), measure.vars = c("SAD averted"))
  melted$value = as.numeric(melted$value)
  Table2a = cast(melted,variable+sex+util~txeff)
  Table2a[,4:10] = round(Table2a[,4]-Table2a[,4:10])
  
  melted = melt(LYGdf, id.vars = c("txeff","util","sex"), measure.vars = c("LYG"))
  melted$value = as.numeric(melted$value)
  Table2b = cast(melted,variable+sex+util~txeff)
  
  Table2 = rbind(Table2a, Table2b)
  Table2 = Table2[,-c(4)] #get rid of baseline
  colnames(Table2) = c("outcome","sex","util", "Any Tx","Pharm Tx","100% increase","150% increase","200% increase","MPC")
  Table2[,10:14]=round(Table2[,4:8]/Table2[,9]*100,1)
  Table2[,4:9] = format(Table2[,4:9],big.mark=",", trim=TRUE)
  
  for (r in 1:nrow(Table2))  {
    Table2[r,10:14] = paste0(format(Table2[r,4:9],big.mark=",", trim=TRUE)," (",Table2[r,10:14],")")
  }
  return(Table2)  
}


createPrevTable<-function(prevdf,whichprev){
  melted = melt(prevdf, id.vars = c("txeff","util","sex"), measure.vars = c(paste0(whichprev,seq(2020,2100,20))))
  melted$value = as.numeric(melted$value)
  eTable1F = cast(subset(melted,sex=="females"),txeff+variable~util)
  eTable1M = cast(subset(melted,sex=="males"),txeff+variable~util)
  eTable1 = cbind(eTable1F,eTable1M[3:6])
  colnames(eTable1) = c("TxF","yearF",paste0(utilset,"F"), paste0(utilset,"M"))
  return(eTable1)
}

Table2 <- createTable2(dfs[[1]],dfs[[2]])
eTable4 <- createTable2(dfs[[3]],dfs[[4]])
eTable5 <- createTable2(dfs[[5]],dfs[[6]])
eTable6 <- createTable2(dfs[[7]],dfs[[8]])
eTable1 <- createPrevTable(dfs[[9]],"depprev") # Smoking prevalence among Current MDE
eTable2 <- createPrevTable(dfs[[10]],"fdepprev")
eTable3 <- createPrevTable(dfs[[11]],"recallprev")

wb <- createWorkbook()
addWorksheet(wb, "Table2")
writeDataTable(wb, "Table2", Table2, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb, "eTable4")
writeDataTable(wb, "eTable4", eTable4, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb, "eTable5")
writeDataTable(wb, "eTable5", eTable5, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb, "eTable6")
writeDataTable(wb, "eTable6", eTable6, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb, "eTable 1")
writeDataTable(wb, "eTable 1", eTable1, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb, "eTable2")
writeDataTable(wb, "eTable2", eTable2, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb, "eTable3")
writeDataTable(wb, "eTable3", eTable3, colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, paste0("Tables",namethisrun,".xlsx"), overwrite = TRUE)