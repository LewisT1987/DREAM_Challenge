#### Preparing Data

## Load required packages
require(Biobase)
require(caret)
require(ggplot2)
require(pls)
library(doMC)
require(pls)
require(glmnet)
require(pamr)
require(plyr)
require(caret)
require(pROC)
registerDoMC(cores = 6)

## load and process P1 data ----
load('data/Eset.Training.P1.sva.adj.2016-06-10.Rda')
eset.P1 <- sva.adj.eset.P1

#eset.INF <- eset.INF[, eset.INF$trt == "Tofacitinib"]

# Create -24 , 0 and 2 hr data
hr.neg24<-which(eset.P1$TIMEHOURS==-24)
hr.0<-which(eset.P1$TIMEHOURS==0)
hr.2<-which(eset.P1$TIMEHOURS==2)

X.neg24<-exprs(eset.P1)[,hr.neg24]
X.0<-exprs(eset.P1)[,hr.0]
X.2<-exprs(eset.P1)[,hr.2]

rownames(X.neg24)<-paste(rownames(X.neg24), "_Hneg24", sep = "")
rownames(X.0)<-paste(rownames(X.0), "_H0", sep = "")

important.variables<-c('SUBJECTID', 'SHEDDING_SC1', 'Virus', 'Study')

DataHneg24<-cbind.data.frame(pData(eset.P1[,hr.neg24])[,important.variables], t(X.neg24))
save(DataHneg24, file='data/Classification/Data_Hneg24.rdata')

DataH0<-cbind.data.frame(pData(eset.P1[,hr.0])[,important.variables], t(X.0))
save(DataH0, file='data/Classification/Data_H0.rdata')

DataH2<-cbind.data.frame(pData(eset.P1[,hr.2])[,important.variables], t(X.2))
save(DataH2, file='data/Classification/Data_H2.rdata')

Both<-intersect(sort(DataH0$SUBJECTID),sort(DataHneg24$SUBJECTID))

#DataH0<-DataH0[which(DataH0$SUBJECTID%in%Both),]
#DataHneg24<-DataHneg24[which(DataHneg24$SUBJECTID%in%Both),]

DataHneg24_H0<-merge(DataHneg24, DataH0, by=c('SUBJECTID','SHEDDING_SC1','Virus','Study'), all.x=FALSE, all.y = FALSE)
save(DataHneg24_H0, file='data/Classification/Data_Hneg24_H0.rdata')

