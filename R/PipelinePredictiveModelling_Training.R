#' PipelinePredictiveModelling_Training
#'
#' @description Processes the data and runs training algorithms
#' @param Biomarker.Set: 'INF' or 'CVD.INF'
#' @param Drug: 'Tofa' or 'Etan'
#' @export

PipelinePredictiveModelling_Training<-function(Time.point='Hneg24',
                                      n.bal.samp=50,
                                      repl=100,
                                      ncores=6,
                                      number.folds=10,
                                      strategy=NULL,
                                      folder.results='results/',
                                      folder.outputs='outputs/',
                                      folder.figs='figs/',
                                      folder.tables='tables/',
                                      feat.num=200,
                                      myseed=myseed){

  ########## Parameters to help Debugging ##############
  # Time.point=c('Hneg24','H0')
  # n.bal.samp=3
  # repl=3
  # ncores=6
  # number.folds=10
  # strategy='no_bootstrap'
  # folder.results='results/'
  # folder.outputs='outputs/'
  # folder.figs='figs/'
  # folder.tables='tables/'
  # myseed=1234
  ######################################################

  ### 1 Preparing Data ### ----

  basic.fname<-paste('Data',paste(Time.point, collapse='_'),sep = "_")
  file.data=paste0('data/Classification/',basic.fname,'.rdata')
  Data<-get(load(file.data))

  dependent.variables<-lapply(Time.point, function(x) grep(x, colnames(Data)))
  dependent.variables<-unlist(dependent.variables)
  save.fname<-paste0(folder.results,folder.outputs ,basic.fname)
  Ya=Data[,'SHEDDING_SC1']
  names(Ya)<-as.character(Data[,'SUBJECTID'])
  Ya=mapvalues(Ya,from=c(1,0),to=c('Sh','NSh'))
  Xa=t(Data[,dependent.variables])
  colnames(Xa)<-as.character(Data[,'SUBJECTID'])

  ### 2 Run Training Algorithms ### ----

  RES.training<-Run_Training_Algorithms(Ya,Xa, replic=repl,  number.folds=number.folds, n.balanced.samples=n.bal.samp, myseed=myseed,
                                        fileout=paste0(save.fname,'_train','_',strategy,'.Rda',sep = ""), ncores=ncores, strategy=strategy, run.bagging = FALSE, bagged.genes = NULL, feat.num=feat.num)

}

