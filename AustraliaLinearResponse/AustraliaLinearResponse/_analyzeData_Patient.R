# Analyze real patietn data.
#JL 2017

###((Run 'importData_Simulation.R', followed by 'generateSimulatedData.R'###

library(Hmisc)

# vary the RECIST cutoffs and see what happens.

# set the file names


    tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
    tmpProjName = 'AustraliaLinearResponse'
    tmpDataDir = file.path(tmpBaseDir, tmpProjName)
    tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

    source(file.path(tmpScriptDir, '_importData_Patient.R'))
#Adding another column to the data that has the minimum change. This will be used for all bestClass computations
dataFullColsAnonymized$minRecistChange = sapply(dataLinRECIST, `[[`, "changeMin")
#Which cutoffs are being tested
corr.Cuts = seq(-1, 1, by = 0.05)

        tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
        tmpProjName = 'AustraliaLinearResponse'
        tmpDataDir = file.path(tmpBaseDir, tmpProjName)
        tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

# do the variable cutpoint thing for the whole dataset

## changing the tThresh for calculating the other metrics
tThreshArray = seq(
    from = floor(quantile(dataFullColsAnonymized$SurvTimeDays, c(0.05, 0.95)))[[1]],
    to = quantile(dataFullColsAnonymized$SurvTimeDays, c(0.05, 0.95))[[2]],
    by = 50)
# empty list to fill with data later when varying tThresh
variedtThreshData=list()
for (i in c(1:length(tThreshArray))) {
    tThresh <- tThreshArray[[i]]

        source(file.path(tmpScriptDir, '_funcPerfMetrics_Patient.R'))
        tmpOutput = perf.metrics(corr.Cuts, dataFullColsAnonymized, dataLinRECIST, classRECIST,tThresh)
        corr.bestClass = tmpOutput$corr.bestClass
        corr.firstClass = tmpOutput$corr.firstClass
        corr.lastClass = tmpOutput$corr.lastClass
        logrank.bestClass = tmpOutput$logrank.bestClass
        logrank.firstClass = tmpOutput$logrank.firstClass
        logrank.lastClass = tmpOutput$logrank.lastClass

    variedtThreshData[[toString(tThresh)]][["npv"]] <- corr.bestClass$npv
    variedtThreshData[[toString(tThresh)]][["ppv"]] <- corr.bestClass$ppv
    variedtThreshData[[toString(tThresh)]][["spec"]] <- corr.bestClass$spec
    variedtThreshData[[toString(tThresh)]][["sens"]] <- corr.bestClass$sens
    variedtThreshData[[toString(tThresh)]][["npvNonResp"]] <- corr.bestClass$npvNonResp
    variedtThreshData[[toString(tThresh)]][["ppvProgFree"]] <- corr.bestClass$ppvProgFree
    variedtThreshData[[toString(tThresh)]][["specNonResp"]] <- corr.bestClass$specNonResp
    variedtThreshData[[toString(tThresh)]][["sensProgFree"]] <- corr.bestClass$sensProgFree

    npv0 = corr.bestClass$npv
    npv0[is.na(npv0)] <-0
    ppv0 = corr.bestClass$ppv
    ppv0[is.na(ppv0)] <- 0
    spec0 = corr.bestClass$spec
    spec0[is.na(spec0)] <- 0
    sens0 = corr.bestClass$sens
    sens0[is.na(sens0)] <- 0
    npvNonResp0 = corr.bestClass$npvNonResp
    npvNonResp0[is.na(npvNonResp0)] <- 0
    ppvProgFree0 = corr.bestClass$ppvProgFree
    ppvProgFree0[is.na(ppvProgFree0)] <- 0
    specNonResp0 = corr.bestClass$specNonResp
    specNonResp0[is.na(specNonResp0)] <- 0
    sensProgFree0 = corr.bestClass$sensProgFree
    sensProgFree0[is.na(sensProgFree0)] <- 0
    variedtThreshData[[toString(tThresh)]][["total"]] <- npv0 + ppv0 + spec0 + sens0 + npvNonResp0 + ppvProgFree0 + specNonResp0 + sensProgFree0
}
rm(i)
rm(npv0, ppv0, spec0, sens0, npvNonResp0, ppvProgFree0, specNonResp0, sensProgFree0)

c1 = 15
c2 = 24

temp = matrix(nrow = length(tThreshArray), ncol = 2)
for (i in c(1:length(tThreshArray))) {  
    temp[[i,1]]=lapply(lapply(lapply(variedtThreshData, `[[`, "sensProgFree"), '[[', c1, c2), '[[', 1)[[i]]
    temp[[i, 2]] = lapply(lapply(lapply(variedtThreshData, `[[`, "specNonResp"), '[[', c1, c2), '[[', 1)[[i]]   
}
temp <- na.omit(temp)
plot(1 - temp[, 1], temp[, 2], pch = 20, xlab = "1 - SensitivityProgFree", ylab = "SpecificityNonResp", main = "SpecificityNonResp vs 1-SensitivityProgFree \nfor various tTHresh", ylim = c(0, 1), xlim = c(0, 1))

temp = matrix(nrow = length(tThreshArray), ncol = 2)
for (i in c(1:length(tThreshArray))) {
    temp[[i, 1]] = lapply(lapply(lapply(variedtThreshData, `[[`, "ppvProgFree"), '[[', c1, c2), '[[', 1)[[i]]
    temp[[i, 2]] = lapply(lapply(lapply(variedtThreshData, `[[`, "npvNonResp"), '[[', c1, c2), '[[', 1)[[i]]
}
temp <- na.omit(temp)
plot(temp[, 1], temp[, 2], pch = 20, xlab = "PPVProgFree", ylab = "NPVNonResp", main = "PPVProgFree vs NPVNonResp \nfor various tTHresh", ylim = c(0, 1), xlim = c(0, 1))

temp = matrix(nrow = length(tThreshArray), ncol = 2)
for (i in c(1:length(tThreshArray))) {
    temp[[i, 1]] = tThreshArray[[i]]
    temp[[i, 2]] = lapply(lapply(lapply(variedtThreshData, `[[`, "specNonResp"), '[[', c1, c2), '[[', 1)[[i]]
}
temp <- na.omit(temp)
plot(temp[, 1], temp[, 2], pch = 20, xlab = "tThresh", ylab = "SpecificityNonResp", main = "SpecificityNonResp for various tTHresh")

temp = matrix(nrow = length(tThreshArray), ncol = 2)
for (i in c(1:length(tThreshArray))) {
    temp[[i, 1]] = tThreshArray[[i]]
    temp[[i, 2]] = lapply(lapply(lapply(variedtThreshData, `[[`, "sensProgFree"), '[[', c1, c2), '[[', 1)[[i]]
}
temp <- na.omit(temp)
plot(temp[, 1], temp[, 2], pch = 20, xlab = "tThresh", ylab = "SensitivityProgFree", main = "SensitivityProgFree for various tTHresh")

temp = matrix(nrow = length(tThreshArray), ncol = 2)
for (i in c(1:length(tThreshArray))) {
    temp[[i, 1]] = tThreshArray[[i]]
    temp[[i, 2]] = lapply(lapply(lapply(variedtThreshData, `[[`, "npvNonResp"), '[[', c1, c2), '[[', 1)[[i]]
}
temp <- na.omit(temp)
plot(temp[, 1], temp[, 2], pch = 20, xlab = "tThresh", ylab = "NPVNonResp", main = "NPVNonResp for various tTHresh")

temp = matrix(nrow = length(tThreshArray), ncol = 2)
for (i in c(1:length(tThreshArray))) {
    temp[[i, 1]] = tThreshArray[[i]]
    temp[[i, 2]] = lapply(lapply(lapply(variedtThreshData, `[[`, "ppvProgFree"), '[[', c1, c2), '[[', 1)[[i]]
}
temp <- na.omit(temp)
plot(temp[, 1], temp[, 2], pch = 20, xlab = "tThresh", ylab = "PPVProgFree", main = "PPVProgFree for various tTHresh")



save(list = c(
            ls(pattern = 'corr.+'),
            ls(pattern = 'logrank.+')),
            file = file.path(tmpDataDir, 'SimulatedDataSetPerformance.rdata'))
rm(list = ls(pattern = 'tmp.*'))
rm(list = ls(pattern = 'i[[:upper:]]'))
