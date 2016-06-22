#' Run_Training_Algorithms
#' @description Run ALgorithms with bootstrap, patients can occur in both test and train
#' @export

Run_Training_Algorithms<-function(Ya,
                                  Xa,
                                  replic=10,
                                  n.balanced.samples=5,
                                  myseed=1234,
                                  fileout="Name.not.specified.Rda",
                                  run.bagging=FALSE,
                                  bagged.genes=NULL,
                                  number.folds=10,
                                  strategy=NULL,
                                  feat.num=200,
                                  ncores = 1)
{


  ### testing rows ###
  # Ya = Ya
  # Xa = Xa
  # replic = 5
  # n.balanced.samples=5
  # myseed=1523
  # fileout = 'Name.not.specified.Rda'
  # run.bagging = FALSE
  # shuffle = TRUE
  # bagged.genes = NULL
  # ncores = 6
  # indRow=1
  # number.folds=10
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
    sub.samples <- sample(names(Ya),2*length(which(Ya=='NSh')), replace = FALSE, prob = probs) ### creates a vector of patientIDs, containing all NR and an equal number of R

    Y.bal = Ya[sub.samples] ### creates table of R and NR for patient IDs in resamples

    if(run.bagging==FALSE){
      X.bal = Xa[,sub.samples] # X is a data frame for the expression data for each patient in resample
    }else if(run.bagging==TRUE){
      X.bal = Xa[bagged.genes,sub.samples] # this will pick out only the gene data selected for bagging
    }

    X.bal<-as.matrix(scale(X.bal))  ### watch out; we are scaling the Variables

    X.bal.samp<-sample(rownames(X.bal),feat.num, replace=FALSE) ### creates a vector of patientIDs, containing all NR and an equal number of R

    X.bal<-X.bal[which(rownames(X.bal)%in%X.bal.samp),]

    # tuning the model
    ctrl <- caret::trainControl(method = "repeatedcv",
                                number = number.folds, #number = number of folds ie 10.
                                repeats = 1,
                                selectionFunction = "best",
                                savePredictions = TRUE,
                                classProbs = TRUE,
                                allowParallel = TRUE) # best train assesed using 10-fold cross validation


    if(strategy=='no_bootstrap'){
      X.bal.resample<-X.bal
      Y.bal.resample<-Y.bal
      }else if(strategy=='bootstrap'){
      ## CODE FOR BOOTSTRAPPING
      set.seed(myseed+mIdx$i[indRow])
      ind.resamp<-createResample(1:length(Y.bal), times = 1, list = TRUE)
      ind.resamp=unlist(ind.resamp)

      X.bal.resample<-X.bal[,ind.resamp]
      Y.bal.resample<-Y.bal[ind.resamp]
      }

    PatientID.Train.index = caret::createDataPartition(Y.bal.resample, times = 1, p = 0.8) ## creates 1 resample, indexes for the 80% of the samples for training
    PatientID.Train.index = as.vector(unlist(PatientID.Train.index))

    PatientID.Test.index = c(1:length(Y.bal.resample))[-PatientID.Train.index]
    PatientID.Train.names <- names(Y.bal.resample[PatientID.Train.index])
    PatientID.Test.names <- names(Y.bal.resample[PatientID.Test.index])

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
    X.bal.part.t = t(X.bal.part)

    Y.bal.part = Y.bal.resample[PatientID.Train.names]

    plsFit <- caret::train(X.bal.part.t, factor(Y.bal.part), method = "pls", tuneLength = 5, metric = "Accuracy", trControl = ctrl) # train using pls
    glmFit <- caret::train(X.bal.part.t, factor(Y.bal.part), method = "glmnet", tuneLength = 5, metric = "Accuracy", trControl = ctrl) # train using glm
    pamFit <- caret::train(X.bal.part.t, factor(Y.bal.part), method = "pam", tuneLength = 30, metric = "Accuracy", trControl = ctrl) # train using pam

    tune.PLS = plsFit$bestTune # best tuning parameters for PLS
    tune.GLM = glmFit$bestTune # best tuning parameters for GLM
    tune.PAM = pamFit$bestTune # best tuning parameters for PAM

    ytrain = as.numeric(Y.bal.part == "Sh")
    ytest = as.numeric(Y.bal[PatientID.Test.names] == "Sh")

    xtrain = X.bal.part.t
    xtest = t(X.bal[,PatientID.Test.names])

    PatientID.train<-rownames(xtrain)
    data.train = cbind.data.frame(xtrain, ytrain = ytrain, PatientID.train=PatientID.train)
    data.test = cbind.data.frame(xtest, ytest = ytest, PatientID.test=PatientID.Test.names)

    data = list(x = t(xtrain), y = ytrain, genenames = colnames(xtrain), geneid = colnames(xtrain))
    pam.fit = pamr::pamr.train(data = data, threshold = tune.PAM[, 1])

    # prediction
    ### Joel used this before to obtain predictions
    #pls.fit = cppls(ytrain ~ xtrain, ncomp = tune.PLS[, 1], scale = FALSE)
    #pred.train.pls.old<-predict(pls.fit, xtrain, type = "scores")[, , 1]
    #pred.train.pls.old[which(pred.train.pls.old > 1)] <- 1
    #pred.train.pls.old[which(pred.train.pls.old < 0)] <- 0

    # training sample
    pred.train.pls = predict(plsFit$finalModel, X.bal.part.t, type = "prob")[,'Sh',1]
    pred.train.glm = predict(glmFit$finalModel, X.bal.part.t, type = "response",s=glmFit$bestTune$lambda)
    pred.train.pam = pamr.predict(pam.fit, newx = X.bal.part, threshold = tune.PAM[, 1], type = "posterior")[, 2]
    pred.train.ensemble <- (pred.train.glm + pred.train.pam + pred.train.pls)/3

    train.ensemble.preds.db<-data.frame(Y=ytrain,  p.glm=pred.train.glm[,1], p.pam=pred.train.pam, p.pls=pred.train.pls)
    start.values.ensemble=rep(1/(ncol(train.ensemble.preds.db)-1), (ncol(train.ensemble.preds.db)-1))
    train.logistic.ensemble<-glm(Y~0+., data=train.ensemble.preds.db, family='binomial',
                                 control = glm.control(maxit=500),method='glm.fit2',start=start.values.ensemble)
    datapred.train.logistic.ensemble <- predict(train.logistic.ensemble, train.ensemble.preds.db, type='response')

    # testing sample
    pred.test.pls = predict(plsFit$finalModel, xtest, type = "prob")[,'Sh', 1]
    pred.test.glm = predict(glmFit$finalModel, xtest, type = "response",s=glmFit$bestTune$lambda)
    pred.test.pam = pamr.predict(pam.fit, newx = t(xtest), threshold = tune.PAM[, 1], type = "posterior")[, 2]
    pred.test.ensemble <- (pred.test.glm + pred.test.pam + pred.test.pls)/3
    test.ensemble.preds.db<-data.frame(Y=ytest,p.glm=pred.test.glm,
                                         p.pam=pred.test.pam,
                                         p.pls=pred.test.pls)
    test.ensemble.preds.db<-data.frame(Y=ytest, p.glm=pred.test.glm[,1],  p.pam=pred.test.pam, p.pls=pred.test.pls)
    pred.test.logistic.ensemble <- predict(train.logistic.ensemble, test.ensemble.preds.db, type='response')

    # betas coefficients
    beta.pls = plsFit$finalModel$coefficients[,'Sh',1]
    beta.glm = as.matrix(predict(glmFit$finalModel,type='coefficient',s=glmFit$bestTune$lambda))[-1,]
    betas.all <- rbind(beta.pls, beta.glm)

    # beta coefficinets for plotting
    beta.plot.glm = subset(data.frame(Beta=as.matrix(coef(glmFit$finalModel, s=glmFit$bestTune$lambda))[,1]),Beta!=0);
    beta.plot.pls = subset(data.frame(plsFit$finalModel$coefficients[,'Sh',1]))
    betas.plot.all <- list(best.pls.betas = beta.plot.pls, best.glm.betas = beta.plot.glm)

    # gene selection by loading weights
    wp = plsFit$finalModel$loadings[,1] # loadins for each feature in the first principle component of pls results
    pls.genes = names(which(abs(wp) > quantile(abs(wp), 0.7))) # selects features with loadings greater than 70% quantile
    pam.genes = as.character(pamr.listgenes(pam.fit, data = data, threshold = tune.PAM[, 1])[, 1])
    glm.genes = names(which(beta.glm != 0))
    ensemble.genes = unique(c(pls.genes, pam.genes, glm.genes))

    data.train = cbind(xtrain, ytrain)
    data.test = cbind(xtest, ytest)

    pred.train = cbind.data.frame(pls = pred.train.pls, glm = pred.train.glm[,1], pam = pred.train.pam, ensemble = pred.train.ensemble[,1], ensemble.log =  datapred.train.logistic.ensemble, PatientID=names(pred.train.pls))
    pred.test = cbind.data.frame(pls = pred.test.pls, glm = pred.test.glm[,1], pam = pred.test.pam, ensemble = pred.test.ensemble[,1], ensemble.log = pred.test.logistic.ensemble, PatientID=names(pred.test.pls))
    genes = list(pls.genes = pls.genes, pam.genes = pam.genes, glm.genes = glm.genes, ensemble.genes = ensemble.genes)

    print(indRow)

    return(list(data.train = data.train, data.test = data.test, pred.train = pred.train, pred.test = pred.test, genes = genes, betas.plot = betas.plot.all))
}

  KDX<-seq(1,dim(mIdx)[1],1)
  RES<-parallel::mclapply(KDX,modelEval, mc.cores = ncores)
  save(RES, file=fileout)
}
