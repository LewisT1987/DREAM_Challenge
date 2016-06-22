#' @title BagGenes_for_Method_Freq
#'
#' @description <insert desc>
#' @export
#' @return NULL


BagGenes_for_Method_Freq<-function(fileInput=NULL, bagging.method=NULL, bagging.cutoff=0.75, save.folder=NULL, basic.fname=NULL,ncores=6){

   #Functions
  create.Indicator.Vector<-function(genes,all.genes){
    out<-matrix(0,  nrow=length(all.genes), ncol=1, dimnames=list(all.genes,''));
    out[match(genes, all.genes),1]<-1
    return(out)}

  aux<-function(genes.per.method,all.genes){
    create.Indicator.Vector(genes.per.method,all.genes)[,1]
  }

  mysubstr<-function(x){substr(x,5,nchar(x))}
  getSetType<-function(x){return(unlist(strsplit(x,'_'))[3])}
  getSymbolName<-function(x){y<-unlist(strsplit(x,'_'))[2:3];
                             return(paste(y[1],y[2], sep='_'))}
  getKit<-function(x){factor(grep('*',x), labels=c('INF','CVD'))}
  ###

  save.fname.fig<-paste0(save.folder, 'figs/' ,basic.fname, '_',as.character(bagging.cutoff),'_')
  save.fname.table<-paste0(save.folder, 'tables/', basic.fname,'_', as.character(bagging.cutoff),'_')
  save.fname.outputs<-paste0(save.folder, 'outputs/', basic.fname,'_', as.character(bagging.cutoff),'_')

  DB.Train<-get(load(fileInput))

  all.genes<-setdiff(colnames(DB.Train[[1]]$data.train),'ytrain')
  D.genes<-mclapply(DB.Train, function(l,all.genes){
                                sapply(l$genes, aux, all.genes=all.genes, simplify=TRUE, USE.NAMES=TRUE)},
    all.genes, mc.cores=ncores)

  bagging.freq.genes<-sapply(colnames(D.genes[[1]]), function(method,L){
    Ind.matrix<-sapply(L, function(l,method){l[,method]}, method,simplify=TRUE, USE.NAMES=TRUE) ## Here!
    rowMeans(Ind.matrix)
  },
  D.genes, simplify=TRUE, USE.NAMES=TRUE)

  feature.name<-as.character(rownames(bagging.freq.genes))
  A<-cbind.data.frame(bagging.freq.genes,
                      feature.names=feature.name,
                      Symbol=gsub('_baseline','_W0', sapply(feature.name,getSymbolName)),
                      type=sapply(feature.name,getSetType),
                      kit=ifelse(grepl('*',feature.name,fixed=TRUE),'CVD','INF'))

  ##########################################
  ### PLEASE DO this in list form as to not be dependent on number of methods; CHECK DoSelectionSelection ######
  getBaggedGenes<-function(m, bagging.freq.genes,bagging.method,bagging.cutoff){
    b.freq<-bagging.freq.genes[,m]
    if (bagging.method=='perc'){ wich.genes<-which[b.freq>bagging.cutoff]
    } else if (bagging.method=='quant'){
      wich.genes<-which(b.freq > quantile(b.freq, bagging.cutoff))
    }
    return(rownames(bagging.freq.genes)[wich.genes])
  }

  bagged.genes<-sapply(colnames(bagging.freq.genes), getBaggedGenes, bagging.freq.genes, bagging.method,bagging.cutoff, USE.NAMES=TRUE, simplify=FALSE)

## Run plots for each methods
  # for (m in names(bagged.genes)){
  #   temp<-subset(cbind(A, V1=A[,m]), feature.names%in%bagged.genes[[m]])
  #   Num.Markers = nrow(temp)
  #   type.names<-table(temp$type)
  #   title.name = paste(m,  paste(names(type.names),'=',type.names, sep=' ', collapse = ', '))
  #   pred.num=length(table)
  #
  #   GX<-ggplot(temp,
  #            aes(x=reorder(Symbol,as.numeric(as.character(V1))),y=as.numeric(as.character(V1)),fill=type))+
  #   geom_bar(stat = 'identity',position = position_dodge(width=0.25)) +
  #   scale_fill_manual(values=c('skyblue','darkred','darkgreen')) +
  #   labs(title = title.name, y= 'Frequency of Selection %',x='Proteins') +
  #   coord_flip(ylim=c(0,1))
  # ggsave(GX,file= paste0(save.fname.fig,'_', m ,'_Selection_Freq.pdf'), width=6, height= (0.162*Num.Markers), limitsize = FALSE)
  # }

  write.csv(bagging.freq.genes, file=paste0(save.fname.table,'selection_frequencies.csv'),row.names = TRUE)
  save(bagging.genes=bagged.genes, file=paste0(save.fname.outputs,bagging.method,bagging.cutoff, '_BaggedGenes','.Rda'))

  return(bagged.genes)
}
