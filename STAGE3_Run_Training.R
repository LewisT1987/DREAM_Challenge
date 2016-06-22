## This code will run the training algorithms for W0, W4 and W0&W4

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
library(glm2)
registerDoMC(cores = 6)

#require(DREAM.Package)
#source('R/PipelinePredictiveModelling_Training.R')
#source('R/Run_Training_Algorithms.R')

#### Run Training Pipeline #### ----

## set training parameters ##

n.balanced.samples=50 # number of balanced patient resamples
n.replicates=100 # number of 80/20 splits
ncores=6
myseed=7243
strategy='no_bootstrap' ## Args 'no_bootstrap' or 'bootstrap'
feat.num=300

######### Hneg24 ######### ----

### Run Piepline Function Tofa 'INF', 'W0' ###
PipelinePredictiveModelling_Training(Time.point='Hneg24',
                                     n.bal.samp=n.balanced.samples,
                                     repl=n.replicates,
                                     ncores=ncores,
                                     strategy=strategy,
                                     feat.num=feat.num,
                                     myseed=myseed)

######### H0 ######### ----

### Run Piepline Function Tofa 'INF', 'W0' ###
PipelinePredictiveModelling_Training(Time.point='H0',
                                     n.bal.samp=n.balanced.samples,
                                     repl=n.replicates,
                                     ncores=ncores,
                                     strategy=strategy,
                                     feat.num=feat.num,
                                     myseed=myseed)

######### Hneg24H2 ######### ----

### Run Piepline Function Tofa 'INF', 'W0' ###
PipelinePredictiveModelling_Training(Time.point=c('Hneg24','H0'),
                                     n.bal.samp=n.balanced.samples,
                                     repl=n.replicates,
                                     ncores=ncores,
                                     strategy=strategy,
                                     feat.num=feat.num,
                                     myseed=myseed)

