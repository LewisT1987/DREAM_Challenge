#' extractPerformance_bag
#'
#' @export

extractPerformance_bag<-function(L){


  if(length(L)>1){
    obs = factor(L$data.test[,'ytest'],levels=c(0,1))
    pam = factor(ifelse(L$pred.test$pam>0.5,1,0),levels=c(0,1))
    pls = factor(ifelse(L$pred.test$pls>0.5,1,0),levels=c(0,1))
    glm = factor(ifelse(L$pred.test$glm>0.5,1,0),levels=c(0,1))
    ensemble = factor(ifelse(L$pred.test$ensemble>0.5,1,0),levels=c(0,1))
    ensemble.log = factor(ifelse(L$pred.test$ensemble.log>0.5,1,0),levels=c(0,1))

    pam.acc=caret::confusionMatrix(table(obs,pam))$overall
    pls.acc=caret::confusionMatrix(table(obs,pls))$overall
    glm.acc=caret::confusionMatrix(table(obs,glm))$overall
    ensemble.acc=caret::confusionMatrix(table(obs,ensemble))$overall
    ensemble.log.acc=caret::confusionMatrix(table(obs,ensemble.log))$overall

    glm.genes=L$glm.genes
    pls.genes=L$pls.genes
    pam.genes=L$pam.genes
    ensemble.genes=L$ensemble.genes

    rx <-cbind.data.frame(rbind(pam.acc,pls.acc,glm.acc,ensemble.acc,ensemble.log.acc),method=c('pam','pls','glm','ensemble','ensemble.log'))
    rownames(rx)<-NULL

    return(rx)}

  else{rx<-NULL}


}
