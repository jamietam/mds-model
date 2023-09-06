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

load("mdetx_depsmk_2010-2018.rda")
load("mdseprevs0520.rda")

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
  
  p<- ggplot(data=txdf) + 
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

pdf(file = paste0(folder,"Fig2_mhutil", namethisrun,".pdf"),width=10, height=6,onefile = FALSE)
grid_arrange_shared_legend(list(uF, uM),2,"")
dev.off()


# Figure 3 prevalence projections with NSDUH comparison -------------------------------------------------

getprevsFM <- function(m, t,u){
  outF = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, mpc=m, txeffcess = t, util= u)
  
  cs_depF = outF[[23]]
  deppopF = outF[[13]]
  
  cs_deppopF = getmodelprevs(cs_depF,deppopF)

  outM = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, mpc=m, txeffcess = t, util= u)
  
  cs_depM = outM[[23]]
  deppopM = outM[[13]]
  
  cs_deppopM = getmodelprevs(cs_depM,deppopM)

  return(list(cs_deppopF, cs_deppopM))
}

scenarios = c("Baseline", "Any Tx","Pharm Tx","100% increase","150% increase","200% increase","Maximum Potential Cessation")

for(t in 1:length(txeffset)){ # Baseline and treatment scenarios
    allprevs = getprevsFM(0,txeffset[t],"1.0")
    assign(paste0("cs_deppop","F",(t-1)), allprevs[[1]])
    assign(paste0("cs_deppop","M",(t-1)), allprevs[[2]])
}
allprevs = getprevsFM(1,1.0,"1.0") # MPC scenario
assign(paste0("cs_deppop","F",length(scenarios)-1), allprevs[[1]])
assign(paste0("cs_deppop","M",length(scenarios)-1), allprevs[[2]])

Figmodelnsduh <- function(whichgender,dfs,smkgroup,popgroup, title){
  dfs = lapply(dfs, function(x) {as.data.frame(x)})
  dfs = Map(cbind,dfs,scenario=as.list(scenarios))
  
  fullset = bind_rows(dfs, .id = NULL)
  fullset$age <-rep(agerownames,length(scenarios))

  model <- melt(fullset,id.vars=c("age","scenario"))
  
  colnames(model)[4] <- "modelvalue"
  model$variable <- as.numeric(levels(model$variable))[model$variable] # converts factor to numeric
  model$scenario  <- factor(model$scenario , levels = scenarios[1:length(scenarios)] )
  
  
  nsduh <- getnsduhprevs(mdseprevs,whichgender,smkgroup,popgroup)
  nsduh$age<-agerownames
  nsduh <- melt(as.data.frame(nsduh),id.vars=c("age"))
  colnames(nsduh)[3] <- c("nsduhvalue")
  nsduh$variable <- as.numeric(levels(nsduh$variable))[nsduh$variable] # converts factor to numeric
  
  fig <- ggplot() +
    geom_point(data=subset(nsduh, age=="total"), aes(x = variable, y = nsduhvalue*100, shape="National Survey on Drug Use and Health"))+
    geom_line(data = subset(model,age=="total"), aes(x=variable, y= modelvalue*100, colour=scenario,group=rev(scenario)))+
    scale_color_viridis(discrete=TRUE, begin=0, end = 0.75 , option="magma") +
    scale_y_continuous(name="Smoking prevalence (%)",limits=c(0,50),breaks=seq(0,50,5)) +
    scale_x_continuous(name="Year",limits=c(2005,endyear),breaks=c(2005,seq(2020,endyear, 10)))  +
    labs(title=title)+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.title = element_blank())
  return(fig)
}

dfsF = list(cs_deppopF0,cs_deppopF1,cs_deppopF2,cs_deppopF3,cs_deppopF4,cs_deppopF5, cs_deppopF6)
Figure3a <- Figmodelnsduh("Women",dfsF,"currentsmoker","deppop","A) Women with Current MDE")
  
dfsM = list(cs_deppopM0,cs_deppopM1,cs_deppopM2,cs_deppopM3,cs_deppopM4,cs_deppopM5, cs_deppopM6)
Figure3b <- Figmodelnsduh("Men",dfsM,"currentsmoker","deppop", "B) Men with Current MDE")

pdf(file = paste0(folder,"Fig3_smkprev", namethisrun,".pdf"),width=10, height=6, onefile=FALSE)
grid_arrange_shared_legend(list(Figure3a, Figure3b),2,"")
dev.off()

jpeg(file = paste0(folder,"Fig3_smkprev", namethisrun,".jpg"),width=10, height=6, units="in", res=500)
grid_arrange_shared_legend(list(Figure3a, Figure3b),2,"")
dev.off()

# Tables - scenario outcomes -------------------------------------------

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

eachscen <- function(outF0, y1,y2, death_nsF, death_csF, death_fsF){

  cs_depF0 = outF0[[23]]
  fs_depF0 = outF0[[24]]
  depSADF0 = cs_depF0[,(policystart-startyear+1):(y2-startyear+1)]*(death_csF[,(policystart-startyear+1):(y2-startyear+1)]-death_nsF[,(policystart-startyear+1):(y2-startyear+1)])+
    fs_depF0[,(policystart-startyear+1):(y2-startyear+1)]*(death_fsF[,(policystart-startyear+1):(y2-startyear+1)]-death_nsF[,(policystart-startyear+1):(y2-startyear+1)])
  cSADdep = sum(depSADF0,na.rm=TRUE)
  
  deppopF0 = outF0[[13]]
  prev_dep = round(getmodelprevs(cs_depF0,deppopF0)[6, paste(c(policystart, seq(y1,y2,20)))]*100,1)
    
  return(list(cSADdep, depSADF0,prev_dep))
}

outF0 = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, mpc = 0, txeffcess=1.0, util="1.0") 
outM0 = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, mpc = 0, txeffcess=1.0, util="1.0")

depSADF0 <- eachscen(outF0, 2040,endyear, death_nsF, death_csF, death_fsF)[[2]] # current MDE
depSADM0 <- eachscen(outM0, 2040,endyear, death_nsM, death_csM, death_fsM)[[2]]

depSADtable = NULL
depLYGtable = NULL
depprevtable =  NULL

for (t in txeffset){
  for (u in utilset){
    for (m in c(0:1)){
      if(m>0&t>1.0) next
      outF = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, mpc=m, txeffcess = t, util=u)
      outM = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, mpc=m, txeffcess = t, util=u)
      
      scenresultsF = eachscen(outF, 2040,endyear, death_nsF, death_csF, death_fsF)
      scenresultsM = eachscen(outM, 2040,endyear, death_nsM, death_csM, death_fsM)
      
      depSADtable = rbind(depSADtable,c("females",m,t,u,scenresultsF[[1]]),c("males",m,t,u,scenresultsM[[1]]))
      depLYGtable = rbind(depLYGtable,
                          c("females",m,t,u,round(sum((LE_fsF[,(policystart-startyear+1):(endyear-startyear+1)]-LE_csF[,(policystart-startyear+1):(endyear-startyear+1)])*(depSADF0-scenresultsF[[2]]),0))),
                          c("males",m,t,u,round(sum((LE_fsM[,(policystart-startyear+1):(endyear-startyear+1)]-LE_csM[,(policystart-startyear+1):(endyear-startyear+1)])*(depSADM0-scenresultsM[[2]]),0))))

      depprevtable = rbind(depprevtable,c("females",m,t,u,scenresultsF[[3]]),c("males",m,t,u,scenresultsM[[3]]))
    }
  }
}

dfs = list(depSADtable,depLYGtable,depprevtable)
dfs = lapply(dfs, function(x) {as.data.frame(x)})
colnames(dfs[[1]]) <- c("sex","mpc","txeff","util","SAD averted") 
colnames(dfs[[2]]) <- c("sex","mpc","txeff","util","LYG")
colnames(dfs[[3]]) <- c("sex","mpc","txeff","util",paste0("depprev",seq(policystart,endyear,20)))
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
  melted = melt(prevdf, id.vars = c("txeff","util","sex"), measure.vars = c(paste0(whichprev,seq(policystart,endyear,20))))
  melted$value = as.numeric(melted$value)
  eTable1F = cast(subset(melted,sex=="females"),txeff+variable~util)
  eTable1M = cast(subset(melted,sex=="males"),txeff+variable~util)
  eTable1 = cbind(eTable1F,eTable1M[3:6])
  colnames(eTable1) = c("TxF","yearF",paste0(utilset,"F"), paste0(utilset,"M"))
  return(eTable1)
}

Table2 <- createTable2(dfs[[1]],dfs[[2]])
eTable1 <- createPrevTable(dfs[[3]],"depprev") # Smoking prevalence among Current MDE


wb <- createWorkbook()
addWorksheet(wb, "Table2")
writeDataTable(wb, "Table2", Table2, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb, "eTable 1")
writeDataTable(wb, "eTable 1", eTable1, colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, paste0(folder,"Tables",namethisrun,".xlsx"), overwrite = TRUE)