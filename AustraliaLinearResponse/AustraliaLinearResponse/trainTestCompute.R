train.test=list()
#Percent that will be in training set (>0.5)
percTrainingSet = 0.6
#Getting the train set
train.test$trainSet = head(simLinRECIST, length(simLinRECIST) * percTrainingSet)
#Getting the test set before redefininig the variables
train.test$testSet = tail(simLinRECIST, length(simLinRECIST) * (1 - percTrainingSet))
#Redefining these variables so that I get the trainging set
simDataFullCols <- head(simDataFullCols, (dim(simDataFullCols)[1]) * percTrainingSet)
simLinRECIST <- head(simLinRECIST, length(simLinRECIST) * percTrainingSet)


tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
tmpProjName = 'AustraliaLinearResponse'
tmpDataDir = file.path(tmpBaseDir, tmpProjName)
tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

# do the variable cutpoint thing for the training set
source(file.path(tmpScriptDir, 'funcPerfMetrics.R'))
tmpOutput = perf.metrics(corr.Cuts, simDataFullCols, simLinRECIST, classRECIST)
corr.bestClass = tmpOutput$corr.bestClass
corr.firstClass = tmpOutput$corr.firstClass
corr.lastClass = tmpOutput$corr.lastClass
logrank.bestClass = tmpOutput$logrank.bestClass
logrank.firstClass = tmpOutput$logrank.firstClass
logrank.lastClass = tmpOutput$logrank.lastClass

#Getting the max C-stat, and PR, PD cutoffs at the max C-stat of the train set
train.test$trainSet.maxC = corr.bestClass$max.stats$max.C
train.test$trainSet.PR = corr.bestClass$max.stats$max.C.PR
train.test$trainSet.PD = corr.bestClass$max.stats$max.C.PD

train.test$testSet.TimeChangePair =
            cbind(
                do.call(cbind, list(lapply(train.test$testSet, `[[`, c('SurvTimeDays')))),
                do.call(cbind, list(lapply(lapply(train.test$testSet, `[[`, c('meas')), '[', 2))))

train.test$testSet.PatientsPR = list()
for (i in 1:nrow(train.test$testSet.TimeChangePair)) {
    if (train.test$testSet.TimeChangePair[, 2][[i]] <= (1 + (corr.bestClass$max.stats$max.C.PR))) {
        train.test$testSet.PatientsPR <- c(train.test$testSet.PatientsPR, train.test$testSet.TimeChangePair[i,])
    }
}
rm(i)
train.test$testSet.PatientsPR <- t(matrix(train.test$testSet.PatientsPR, nrow = 2))

train.test$testSet.PatientsPD = list()
for (i in 1:nrow(train.test$testSet.TimeChangePair)) {
    if (train.test$testSet.TimeChangePair[, 2][[i]] >= (1 + (corr.bestClass$max.stats$max.C.PD))) {
        train.test$testSet.PatientsPD <- c(train.test$testSet.PatientsPD, train.test$testSet.TimeChangePair[i,])
    }
}
rm(i)
train.test$testSet.PatientsPD <- t(matrix(train.test$testSet.PatientsPD, nrow = 2))

#Test Set PPV: 
train.test$testSet.ppv = sum(train.test$testSet.PatientsPR[, 1] > 365) / nrow(train.test$testSet.PatientsPR)
#Test Set NPV:
train.test$testSet.npv = sum(train.test$testSet.PatientsPD[, 1] < 365) / nrow(train.test$testSet.PatientsPD)
#Test Set Sensitivity
train.test$testSet.sens = sum(train.test$testSet.PatientsPR[, 1] > 365) / (sum(train.test$testSet.PatientsPR[, 1] > 365) + sum(train.test$testSet.PatientsPD[, 1] > 365))
#Test Set Specificity
train.test$testSet.spec = sum(train.test$testSet.PatientsPD[, 1] < 365) / (sum(train.test$testSet.PatientsPD[, 1] < 365) + sum(train.test$testSet.PatientsPR[, 1] < 365))