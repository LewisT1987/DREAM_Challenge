## load required packages

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
require(WGCNA)
allowWGCNAThreads()


## Run this script to Download and save DREAM Challenge Data
source("http://depot.sagebase.org/CRAN.R")
pkgInstall("synapseClient")

require(synapseClient)
## You can use my login detail if you want
## but, I recomend you visit
## https://www.synapse.org/#!Synapse:syn5647810/wiki/ to make your own
synapseLogin('l.tomalin1987@gmail.com','Dewgon11Jewwej##')

DREAM_Clinical<-synGet('syn6043449')
localFilePath<-DREAM_Clinical@filePath
clin_data<-read.delim(localFilePath,header=T)
rownames(clin_data)<-clin_data$CEL
save(clin_data, file='data/Clinical_Data.Rda')

DREAM_Symptoms<-synGet('syn6043450')
localFilePath<-DREAM_Symptoms@filePath
symp_data<-read.delim(localFilePath,header=T)
save(symp_data, file='data/Symptom_Data.Rda')

Merged_Data<-merge(x=clin_data,y=symp_data)
Merged_Data<-t(Merged_Data)
save(Merged_Data, file='data/Merged_Data.Rda')

DREAM_Expression_norm<-synGet('syn6043448')
localFilePath<-DREAM_Expression_norm@filePath
Exp_norm<-read.delim(localFilePath,header=T)
save(Exp_norm, file='data/Expression_Data_norm.Rda')

DREAM_Expression_raw<-synGet('syn6043347')
localFilePath<-DREAM_Expression_raw@filePath
untar(localFilePath, exdir = "data/Training/CEL")

flsFULL<-list.files("data/Training/CEL",full.names=T)
abatch<-ReadAffy(filenames=flsFULL)
save(abatch, file='data/abatch.Rda')

abatch<-abatch[,rownames(clin_data)]
phenoData(abatch)<-new("AnnotatedDataFrame",clin_data[sampleNames(abatch),])
annotation(abatch)<-'hgu133a2.db'
save(file=paste0("data/Abatch.Training.All.",Sys.Date(),".rda"), abatch)
load("data/Abatch.Training.All.2016-06-08.rda")

DREAM_Annotations<-synGet('syn5684262')
localFilePath<-DREAM_Annotations@filePath
untar(localFilePath, exdir = "data/Training/Annotation")

##### 3. gcRMA eset -----
eset<-gcrma(abatch)
phenoData(eset) <- new("AnnotatedDataFrame", clin_data) # extracting date from CEL files info
eset$Batch.Date<-sapply(pData(protocolData(eset)[sampleNames(eset),])$ScanDate,function(x){substr(x,1,10)},simplify=TRUE,USE.NAMES=FALSE)
eset$Batch.Date.Formatted<-as.Date(ifelse(grepl("/",eset$Batch.Date),as.Date(substr(eset$Batch.Date,1,8),"%m/%d/%y"), as.Date(eset$Batch.Date,"%Y-%m-%d")))
eset$Batch.Year<-format(eset$Batch.Date.Formatted,"%Y")
eset$Batch.Month<-format(eset$Batch.Date.Formatted,"%B")
eset$Virus<-mapvalues(eset$STUDYID,c('DEE1 RSV','DEE2 H3N2','DEE3 H1N1','DEE4X H1N1','DEE5 H3N2','Rhinovirus Duke','Rhinovirus UVA'),c('RSV','H3N2','H1N1','H1N1','H3N2','Rhinovirus','Rhinovirus'))
PD<-pData(eset)
save(file=paste0("data/Eset.Training.All.",Sys.Date(),".rda"),eset,PD)

### Annotation
ann<-read.csv("data/Training/Annotation/HG-U133A_2.na35.annot.csv",stringsAsFactors=F,skip=25)
ann<-subset(ann, select=c(Probe.Set.ID,Gene.Title,Gene.Symbol,Chromosomal.Location,Entrez.Gene,Pathway))
rownames(ann)<-ann$Probe.Set.ID
fData(eset)<-as.data.frame(ann)[featureNames(eset),]

### Colapse Probes
load("data/Eset.Training.All.2016-06-09.rda")
new_exprs<-collapseRows(datET=exprs(eset),
                        rowGroup=fData(eset)$Gene.Symbol,
                        rowID=fData(eset)$Probe.Set.ID)

datETcollap1<-new_exprs[[1]]
datETcollap<-datETcollap1[!grepl('///',rownames(datETcollap1)),]
exprs(eset)<-datETcollap
fData(eset)<-fData(eset)[which(fData(eset)$Gene.Symbol%in%rownames(exprs(eset))),]
save(file=paste0("data/Eset.Training.Collapsed.",Sys.Date(),".rda"),eset,PD)




