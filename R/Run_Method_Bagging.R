#' Run_Method_Bagging
#'
#' @export

Run_Method_Bagging<-function(Ya, Xa,
                             replic = 10,
                             n.balanced.samples = 5,
                             myseed = 1234,
                             fileout = "Name.not.specified.Rda",
                             bagging.genes,
                             strategy=NULL,
                             number.folds = NULL,
                             feat.num=300,
                             ncores = 1)
{


  ### testing rows ###
  # load('../data/olink_baseline_delta_week4.rdata')
  # Ya = Y
  # Xa = X
  # replic = 3
  # n.balanced.samples=3
  # myseed=1523
  # fileout = 'Name.not.specified.Rda'
  # run.bagging = FALSE
  # bagged.genes = bagging.genes
  # ncores = 1
  # indRow=1
  # number.folds = 10
  ##############################

  s = n.balanced.samples
  m = replic

  JDX = seq(1, n.balanced.samples, 1) # creates a vector of length number of balanced samples
  IDX = seq(1, replic, 1) # creates vector length number of replicates in each resample (number of 20:80 splits))

  mIdx <-expand.grid(IDX,JDX);names(mIdx)<-c('i','j') # creates a matrix so that each resample 'j' has 'i' replicates (20:80 splits)
  mIdx$i <- c(1:length(mIdx$j))
  rownames(mIdx)<-1:dim(mIdx)[1] # adds numerical rownames to mIdx


  modelEval<-function(indRow){ # creates function modelEval with variable 'indRow' (indRow is a row name on mIdx)
    # vector of probability of selection
    probs = ifelse(Ya == 'NSh', 0.9, 0.1) ## ensures all NRs are selected 0.9*109 > 42

    # balanced samples
    set.seed(myseed+mIdx$j[indRow])
    sub.samples <- sample(names(Ya),2*length(which(Ya=='NR')), replace = FALSE, prob = probs) ### creates a vector of patientIDs, containing all NR and an equal number of R

    Y.bal=Ya[sub.samples] #Y<-plyr::mapvalues(Y,from=c('NR','R'),to=c(0,1)) ### creates table of R and NR for patient IDs in resamples
    X.bal=as.matrix(scale(Xa[, sub.samples]))

    X.bal.samp<-sample(rownames(X.bal),feat.num, replace=FALSE) ### creates a vector of patientIDs, containing all NR and an equal number of R
    X.bal<-X.bal[which(rownames(X.bal)%in%X.bal.samp),]

    # tuning the model
    ctrl <- caret::trainControl(method = "repeatedcv",
                                number = number.folds, # default number of folds is 10
                                repeats = 1, #default for repeats is 1
                                selectionFunction = "best",
                                savePredictions = TRUE,
                                classProbs = TRUE,
                                allowParallel = TRUE) # caret, criteria to train algorithm 10-fold cross validation

    if(strategy=='no_bootstrap'){
      X.bal.resample<-X.bal
      Y.bal.resample<-Y.bal
    }else {
      set.seed(myseed+mIdx$i[indRow])
      ind.resamp<-createResample(1:length(Y.bal), times = 1, list = TRUE)
      ind.resamp=unlist(ind.resamp)
      X.bal.resample<-X.bal[,ind.resamp]
      Y.bal.resample<-Y.bal[ind.resamp]}

    PatientID.Train.index<-caret::createDataPartition(Y.bal.resample, times = 1, p = 0.8) ## creates 1 resample, indexes for the 80% of the samples for training, from ys (R or NR)
    PatientID.Train.index<-unlist(PatientID.Train.index); names(PatientID.Train.index)=NULL

    PatientID.Test.index<-c(1:length(Y.bal.resample))[-PatientID.Train.index]
    PatientID.Train.names<-names(Y.bal.resample[PatientID.Train.index])
    PatientID.Test.names<-names(Y.bal.resample[PatientID.Test.index])

    if(strategy=='bootstrap'){
      ##### get dupplicates and reasign.
      dup.patients<-intersect(PatientID.Train.names, PatientID.Test.names)
      r<-rbinom(length(dup.patients),size=1,prob=0.5)
      dup.names.train<-dup.patients[(r==1)]
      dup.names.test<-dup.patients[(r==0)]
      PatientID.Train.names<-PatientID.Train.names[PatientID.Train.names%in%setdiff(PatientID.Train.names,dup.names.test)]
      PatientID.Test.names<-PatientID.Test.names[PatientID.Test.names%in%setdiff(PatientID.Test.names,dup.names.train)]
    }

    X.bal.part = X.bal.resample[,PatientID.Train.names]
    Y.bal.part = Y.bal.resample[PatientID.Train.names]
    xtest      = X.bal.resample[,PatientID.Test.names]

    plsFit <- caret::train(t(X.bal.part[bagging.genes$pls.genes,]), factor(Y.bal.part), method = "pls", tuneLength = 5, metric = "Accuracy", trControl = ctrl)
    glmFit <- caret::train(t(X.bal.part[bagging.genes$glm.genes,]), factor(Y.bal.part), method = "glmnet", tuneLength = 5, metric = "Accuracy", trControl = ctrl)
    pamFit <- caret::train(t(X.bal.part[bagging.genes$pam.genes,]), factor(Y.bal.part), method = "pam", tuneLength = 30, metric = "Accuracy", trControl = ctrl)

    tune.PLS = plsFit$bestTune
    tune.GLM = glmFit$bestTune
    tune.PAM = pamFit$bestTune

    ytrain = as.numeric(Y.bal.part == "R")
    ytest = as.numeric(Y.bal[PatientID.Test.names] == "R")

    # predict does not work with pam so we need to do this call.
    data.pam<-list(x = X.bal.part[bagging.genes$pam.genes,],
                   y = ytrain, genenames = rownames(X.bal.part[bagging.genes$pam.genes,]),
                   geneid =rownames(X.bal.part[bagging.genes$pam.genes,]))
    pam.fit = pamr::pamr.train(data=data.pam, threshold = tune.PAM[, 1])

    # prediction
    # training sample
    pred.train.pls = predict(plsFit$finalModel,t(X.bal.part[bagging.genes$pls.genes, ]), type = "prob")[,'R',1] #prob of response
    pred.train.glm = predict(glmFit$finalModel, t(X.bal.part[bagging.genes$glm.genes, ]), type = "response",s=glmFit$bestTune$lambda)[,1] # in glm case, type=response gives you the probability
    pred.train.pam = pamr.predict(pam.fit, newx = X.bal.part[bagging.genes$pam.genes, ], threshold = tune.PAM[, 1], type = "posterior")[, 2]
    pred.train.ensemble <- (pred.train.glm + pred.train.pam + pred.train.pls)/3

    train.ensemble.preds.db<-data.frame(Y=ytrain,  p.glm=pred.train.glm, p.pam=pred.train.pam, p.pls=pred.train.pls)
    start.values.ensemble=rep(1/(ncol(train.ensemble.preds.db)-1), (ncol(train.ensemble.preds.db)-1))
    train.logistic.ensemble<-glm(Y~0+., data=train.ensemble.preds.db, family='binomial',
                                 control = glm.control(maxit=500),method='glm.fit2',start=start.values.ensemble)
    datapred.train.logistic.ensemble <- predict(train.logistic.ensemble, train.ensemble.preds.db, type='response')

    # testing sample
    pred.test.pls = predict(plsFit$finalModel,  t(xtest[bagging.genes$pls.genes,]), type = "prob")[,'R', 1]
    pred.test.glm = predict(glmFit$finalModel,  t(xtest[bagging.genes$glm.genes,]), type = "response",s=glmFit$bestTune$lambda)[,1]
    pred.test.pam = pamr.predict(pam.fit, newx = (xtest[bagging.genes$pam.genes,]), threshold = tune.PAM[, 1], type = "posterior")[, 2]
    pred.test.ensemble <- (pred.test.glm + pred.test.pam + pred.test.pls)/3

    test.ensemble.preds.db<-data.frame(Y=ytest, p.glm=pred.test.glm,  p.pam=pred.test.pam, p.pls=pred.test.pls)
    pred.test.logistic.ensemble <- predict(train.logistic.ensemble, test.ensemble.preds.db, type='response')

    # betas coefficients
    beta.pls = plsFit$finalModel$coefficients[,'R',1]
    beta.glm = as.matrix(predict(glmFit$finalModel,type='coefficient', s=glmFit$bestTune$lambda))[-1,]
    #betas.all <- rbind(beta.pls, beta.glm)

    # beta coefficinets for plotting
    beta.plot.glm = subset(data.frame(Beta=as.matrix(coef(glmFit$finalModel, s=glmFit$bestTune$lambda))[,1]),Beta!=0);
    beta.plot.pls = subset(data.frame(plsFit$finalModel$coefficients[,'R',1]))
    betas.plot.all <- list(pls = beta.plot.pls, glm = beta.plot.glm)

    # proportion of betas different from zero
    ####  prop.0 = length(which(beta.glm!=0))/length(beta.glm) #??? REALLY BAD

    # gene selection by loading weights
    wp = plsFit$finalModel$loadings[,1]
    pls.genes = names(which(abs(wp) > quantile(abs(wp), 0.7)))
    pam.genes = as.character(pamr.listgenes(pam.fit, data = data.pam, threshold = tune.PAM[, 1])[, 1])
    glm.genes = names(which(beta.glm != 0))
    ensemble.genes = unique(c(pls.genes, pam.genes, glm.genes))

    data.train=cbind(colnames(X.bal.part),ytrain)
    data.test=cbind(colnames(xtest),ytest)

    pred.train = cbind.data.frame(pls = pred.train.pls, glm = pred.train.glm, pam = pred.train.pam, ensemble = pred.train.ensemble, ensemble.log = datapred.train.logistic.ensemble, PatientID=names(pred.train.pls))
    pred.test = cbind.data.frame(pls = pred.test.pls, glm = pred.test.glm, pam = pred.test.pam, ensemble = pred.test.ensemble, ensemble.log= pred.test.logistic.ensemble, PatientID=names(pred.test.pls))
    genes = list(pls.genes = pls.genes, pam.genes = pam.genes, glm.genes = glm.genes, ensemble.genes = ensemble.genes)

    print(indRow)

    return(list(data.train=data.train, data.test=data.test, pred.train = pred.train, pred.test = pred.test, genes = genes, betas.plot = betas.plot.all))
  }

  KDX<-seq(1,dim(mIdx)[1],1)
  RES<-parallel::mclapply(KDX, modelEval, mc.cores = ncores)
  save(RES, bagging.genes, file = fileout)
  return(RES)
}
