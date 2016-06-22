# Load required packages
require(Biobase)
require(ggplot2)
require(pls)
library(doMC)
require(pls)
require(glmnet)
require(pamr)
require(plyr)
require(caret)
require(pROC)
require(ROCR)
library(glm2)
registerDoMC(cores = 6)

require(DREAM.Package)

# source('R/PipelinePredictiveModelling_Bagging.R')
# source('R/Run_Method_Bagging.R')
# source('R/BagGenes_for_Method_Freq.R')
# source('R/BagGenes_for_Method_TestAcc.R')
# source('R/extractPerformance_bag.R')

#### Run Bagging Pipeline #### ----

## set training parameters ##

n.balanced.samples=20 # number of patient resamples
n.replicates=40 # number of 80/20 splits
ncores=6
myseed=7243
bagging.method= 'quant' ## IMPORTANT!  arguments 'quant' or 'perc'
strategy='no_bootstrap'
Freq_or_Acc='Acc'
feat.num=300

bagging.cutoff=c(0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95)

### Hneg24 ###

for(i in bagging.cutoff){
  PipelinePredictiveModelling_Bagging(bagging.cutoff=i,
                                      Freq_or_Acc=Freq_or_Acc,
                                      Time.point= c('Hneg24'),
                                      n.bal.samp=n.balanced.samples,
                                      repl=n.replicates,
                                      strategy=strategy,
                                      ncores=ncores,
                                      number.folds=10,
                                      bagging.method=bagging.method,
                                      feat.num=feat.num,
                                      myseed=myseed)}
