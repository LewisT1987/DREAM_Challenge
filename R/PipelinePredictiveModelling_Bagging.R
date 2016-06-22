#' PipelinePredictiveModelling_Bagging
#'
#' @export


PipelinePredictiveModelling_Bagging<-function(bagging.cutoff = NULL,
                                      feat.num=300,
                                      Time.point=NULL,
                                      n.bal.samp = 50,
                                      repl = 100,
                                      ncores = 6,
                                      number.folds = 10,
                                      strategy=NULL,
                                      bagging.method=NULL,
                                      Freq_or_Acc=NULL,
                                      folder.results='results/',
                                      folder.outputs='outputs/',
                                      folder.figs='figs/',
                                      folder.tables='tables/',
                                      myseed=myseed){

########## Parameters to help Debugging ##############
  # bagging.cutoff = 0.5
  # feat.num=300
  # Time.point=c('Hneg24')
  # n.bal.samp = 4
  # repl = 4
  # ncores = 6
  # number.folds = 10
  # bagging.method= 'quant'
  # strategy='no_bootstrap'
  # Freq_or_Acc='Acc'
  # folder.results='results/'
  # folder.outputs='outputs/'
  # folder.figs='figs/'
  # folder.tables='tables/'
  # myseed=6969
######################################################

  ### 1 Preparing Data ### ----

  basic.fname<-paste('Data',paste(Time.point, collapse='_'),sep = "_")
  file.data=paste0('data/Classification/',basic.fname,'.rdata')
  Data<-get(load(file.data))

  dependent.variables<-lapply(Time.point, function(x) grep(x, colnames(Data)))
  dependent.variables<-unlist(dependent.variables)
  save.fname<-paste0(folder.results,folder.outputs ,basic.fname)
  save.folder<-paste0(folder.results)
  Ya=Data[,'SHEDDING_SC1']
  names(Ya)<-as.character(Data[,'SUBJECTID'])
  Ya=mapvalues(Ya,from=c(1,0),to=c('Sh','NSh'))
  Xa=t(Data[,dependent.variables])
  colnames(Xa)<-as.character(Data[,'SUBJECTID'])

 ### 2 Run Ensemble Bagging ### ----
 ### 2.1 Selected genes for Bagging ### ----

  all.genes<-rownames(Xa)

 if(Freq_or_Acc=='Freq'){
 bagging.genes <- BagGenes_for_Method_Freq(fileInput=paste0(save.fname, '_train','_',strategy,'.Rda'), bagging.method=bagging.method, bagging.cutoff=bagging.cutoff,
                                     save.folder=save.folder,basic.fname=paste0(basic.fname,strategy),all.genes=all.genes,ncores=ncores)
 }else if(Freq_or_Acc=='Acc'){
   bagging.genes <- BagGenes_for_Method_TestAcc(fileInput=paste0(save.fname, '_train','_',strategy,'.Rda'), bagging.method=bagging.method, bagging.cutoff=bagging.cutoff,
                                             save.folder=save.folder,basic.fname=paste0(basic.fname,strategy),ncores=ncores)
 }


 ### IMPORTANT pls.genes selection algorithm is subpar; we found that choosing ensemble genes, gives better results.
 bagging.genes[['pls.genes']]<-bagging.genes[['ensemble.genes']]
 #########################

 ### 2.2 Run Algorithms with ensemble_bagging genes ### ----
 # RES.Ensemble.bagging<-Run_Ensemble_Bagging(Ya, Xa, replic=repl,  number.folds = number.folds, n.balanced.samples=n.bal.samp, myseed=myseed,
 #                                                     fileout =paste0(save.fname,'_',bagging.method,'_all_bag.Rda'), ncores=ncores, run.bagging = TRUE,
 #                                                     bagged.genes = ensemble.genes)

 ### 2.3 Run Algorithms with Method-specific bagging ### ----
 RES.Method.bagging<-Run_Method_Bagging(Ya, Xa, replic=repl,  number.folds = number.folds, n.balanced.samples=n.bal.samp, myseed=myseed,
                                    fileout =paste0(save.fname, '_', bagging.method,'_', as.character(bagging.cutoff),'_bag','_', strategy,Freq_or_Acc,'.Rda'),
                                    ncores=ncores, strategy=strategy,bagging.genes=bagging.genes)

}
