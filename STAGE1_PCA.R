## Performs an initial PCA on the expression data

require(magrittr)
require(affy)
require(affycoretools)
require(arrayQualityMetrics)
require(gcrma)
require(preprocessCore)
require(lattice)
require(plyr)
require(zoo)
require(ggplot2)
require(annotate)
require(hgu133a2.db)
require(caret)
require(sva)

load("data/Eset.Training.Collapsed.2016-06-10.rda")

## Extract Phase1 data

#eset.PhaseI<-eset[,(which(pData(eset)$TIMEHOURS<=2 & pData(eset)$STUDYID==c('Rhinovirus Duke','Rhinovirus UVA')))]  # 239 samples
eset.PhaseI<-eset[,(which(pData(eset)$TIMEHOURS<=2))]  # 239 samples

eset.PhaseI$Study<-mapvalues(eset.PhaseI$STUDYID,
                        c('DEE1 RSV','DEE2 H3N2','DEE3 H1N1','DEE4X H1N1','DEE5 H3N2','Rhinovirus Duke','Rhinovirus UVA'),
                        c('DEE1','DEE2','DEE3','DEE4X','DEE5','Duke','UVA'))



# PCA -----
pca.str<-prcomp(t(exprs(eset.PhaseI)))
db<-pca.str$x[,1:7]
var.ex<-round(100*summary(pca.str)$importance[2,1:7],0); var.ex  # proportion of explained variability
db<-cbind.data.frame(db,pData(eset.PhaseI))

ggplot(db,aes(x=PC1,y=PC2, color=STUDYID, group=SUBJECTID, shape=Virus))+
  geom_point(size=4)+
  scale_shape_manual(values=c(19,15,9,8))+
  scale_color_manual(values=c('forestgreen','red', 'blue','darkslategray3','pink3','gold','magenta'))+
  labs(x=paste('PC-1 (',var.ex[1],'%)',sep=''),y=paste('PC-2 (',var.ex[2],'%)',sep=''), title='PCA-Plot Phase I')+
  theme_bw()

# Batch Adj ComBat -----
sva.adj.exprs<-sva::ComBat(dat=exprs(eset.PhaseI), batch=eset.PhaseI$STUDYID, mod=model.matrix(~1, data=pData(eset.PhaseI)))  # using only intercept in the model

sva.adj.eset.P1<-eset.PhaseI
exprs(sva.adj.eset.P1)<-sva.adj.exprs
save(file=paste0("data/Eset.Training.P1.sva.adj.",Sys.Date(),".rda"),sva.adj.eset.P1)


# PCA Batch Adjust -----
pca.str<-prcomp(t(exprs(sva.adj.eset.P1)))
db<-pca.str$x[,1:7]
var.ex<-round(100*summary(pca.str)$importance[2,1:7],0); var.ex  # proportion of explained variability
db<-cbind.data.frame(db,pData(sva.adj.eset.P1))

ggplot(db,aes(x=PC1,y=PC2, color=STUDYID, group=SUBJECTID, shape=Virus))+
  geom_point(size=4)+
  scale_shape_manual(values=c(19,15,9,8))+
  scale_color_manual(values=c('forestgreen','red', 'blue','darkslategray3','pink3','gold','magenta'))+
  labs(x=paste('PC-1 (',var.ex[1],'%)',sep=''),y=paste('PC-2 (',var.ex[2],'%)',sep=''), title='PCA-Plot Phase I')+
  theme_bw()




# Phase I: Training & Test data 70:30 -----
# set.seed(123)
# inTrain<-createDataPartition(y = X.P1$SHEDDING_SC1, ## the outcome data are needed
#                              p = .70, list = FALSE)
# X.P1.train<-X.P1[inTrain,]
# X.P1.test<-X.P1[-inTrain,]
# table(X.P1.test$TIMEHOURS)
#
# eset.P1.train<-eset.PhaseI[,sampleNames(eset.PhaseI)%in%X.P1.train$CEL]
# eset.P1.test<-eset.PhaseI[,sampleNames(eset.PhaseI)%in%X.P1.test$CEL]  # will not use this data until we build models








