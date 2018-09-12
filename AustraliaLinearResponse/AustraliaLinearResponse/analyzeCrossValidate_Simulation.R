# Run cross validation on simulated data.
#
# Run the importData.R and analyzeData_Simulation.R scripts first.
#
# Written by E Gudmundsson, June 2017

tmpBaseDir = file.path('F:', 'MesoResponse')
tmpProjName = 'AustraliaLinearResponse'
tmpDataDir = file.path(tmpBaseDir, tmpProjName)
tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)
library(Hmisc)

# do the variable cutpoint thing each validation subset
corr.Cuts = seq(-1, 1, by = 0.05)
tmpCurPR = match(-0.3, round(corr.Cuts, 5))
tmpCurPD = match(0.2, round(corr.Cuts, 5))
source(file.path(tmpScriptDir, 'funcPerfMetrics.R'))

# the idea with this cross validation is that I will leave out the patients,
# one at a time, and get the "best" model from the other (N-1) patients.
# I will use that model to classify the Nth patient, thereby leading to a
#classification for all N patients by the end of the run.
#So leave each patient out once and do the model fitting.
tmpPatList = names(simLinRECIST)
cv.Full = vector("list", length(tmpPatList))
cv.corr.bestClass = data.frame(patOut = NA,
    max = NA, max.PR = NA, max.PD = NA, max.p = NA,
    std = NA, minP = NA, minP.PR = NA, minP.PD = NA)
cv.corr.firstClass = data.frame(patOut = NA,
    max = NA, max.PR = NA, max.PD = NA, max.p = NA,
    std = NA, minP = NA, minP.PR = NA, minP.PD = NA)
cv.logrank.bestClass = data.frame(patOut = NA, min = NA, min.PR = NA, min.PD = NA)
cv.logrank.firstClass = data.frame(patOut = NA, min = NA, min.PR = NA, min.PD = NA)

for (iPat in tmpPatList) {
    iCount = match(iPat, tmpPatList)

    # make the cross validation patient cohort (leaving out iPat)
    tmpNewPatList = subset(tmpPatList, !(tmpPatList %in% iPat))
    tmp.dataFullCols = subset(simDataFullCols, (URN %in% tmpNewPatList))
    tmp.dataLinRECIST = simLinRECIST[match(tmpNewPatList, tmpPatList)]

    # find the performance for each damn run
    tmpOutput = perf.metrics(corr.Cuts, tmp.dataFullCols, tmp.dataLinRECIST, classRECIST)
    cv.Full[[iCount]] = tmpOutput

    # find the things that matter
    cv.corr.bestClass[iCount, 'patOut'] = iPat
    cv.corr.bestClass[iCount, 'max'] = tmpOutput$corr.bestClass$max
    cv.corr.bestClass[iCount, 'max.PR'] = tmpOutput$corr.bestClass$max.PR
    cv.corr.bestClass[iCount, 'max.PD'] = tmpOutput$corr.bestClass$max.PD
    cv.corr.bestClass[iCount, 'max.p'] = with(tmpOutput$corr.bestClass, perf.p[max.indPR, max.indPD])
    cv.corr.bestClass[iCount, 'std'] = with(tmpOutput$corr.bestClass, perf[tmpCurPR, tmpCurPD])
    cv.corr.bestClass[iCount, 'minP'] = tmpOutput$corr.bestClass$minP
    cv.corr.bestClass[iCount, 'minP.PR'] = tmpOutput$corr.bestClass$minP.PR
    cv.corr.bestClass[iCount, 'minP.PD'] = tmpOutput$corr.bestClass$minP.PD

    cv.corr.firstClass[iCount, 'patOut'] = iPat
    cv.corr.firstClass[iCount, 'max'] = tmpOutput$corr.firstClass$max
    cv.corr.firstClass[iCount, 'max.PR'] = tmpOutput$corr.firstClass$max.PR
    cv.corr.firstClass[iCount, 'max.PD'] = tmpOutput$corr.firstClass$max.PD
    cv.corr.firstClass[iCount, 'max.p'] = with(tmpOutput$corr.firstClass, perf.p[max.indPR, max.indPD])
    cv.corr.firstClass[iCount, 'std'] = with(tmpOutput$corr.firstClass, perf[tmpCurPR, tmpCurPD])
    cv.corr.firstClass[iCount, 'minP'] = tmpOutput$corr.firstClass$minP
    cv.corr.firstClass[iCount, 'minP.PR'] = tmpOutput$corr.firstClass$minP.PR
    cv.corr.firstClass[iCount, 'minP.PD'] = tmpOutput$corr.firstClass$minP.PD

    cv.logrank.bestClass[iCount, 'patOut'] = iPat
    cv.logrank.bestClass[iCount, 'min'] = tmpOutput$logrank.bestClass$min
    cv.logrank.bestClass[iCount, 'min.PR'] = tmpOutput$logrank.bestClass$min.PR
    cv.logrank.bestClass[iCount, 'min.PD'] = tmpOutput$logrank.bestClass$min.PD

    cv.logrank.firstClass[iCount, 'patOut'] = iPat
    cv.logrank.firstClass[iCount, 'min'] = tmpOutput$logrank.firstClass$min
    cv.logrank.firstClass[iCount, 'min.PR'] = tmpOutput$logrank.firstClass$min.PR
    cv.logrank.firstClass[iCount, 'min.PD'] = tmpOutput$logrank.firstClass$min.PD

}

# now, apply those cut points to the left-out cases and see how that does.
cv.Class = simDataFullCols[, c(1, 30, 81, 82, 128)]
cv.Class$bestCuts.bestClass = NA
cv.Class$bestCuts.firstClass = NA
cv.Class$firstCuts.bestClass = NA
cv.Class$firstCuts.firstClass = NA

cv.Class$minPbestCuts.bestClass = NA
cv.Class$minPbestCuts.firstClass = NA
cv.Class$minPfirstCuts.bestClass = NA
cv.Class$minPfirstCuts.firstClass = NA

cv.Class$logBestCuts.bestClass = NA
cv.Class$logBestCuts.firstClass = NA
cv.Class$logFirstCuts.bestClass = NA
cv.Class$logFirstCuts.firstClass = NA
for (iPat in tmpPatList) {
    iCount = match(iPat, tmpPatList)
    iRow = match(as.character(iPat), cv.corr.bestClass$patOut)

    # there are 6 ways to get these cutpoints
    tmp.bestClass.PD = cv.corr.bestClass$max.PD[iRow]
    tmp.bestClass.PR = cv.corr.bestClass$max.PR[iRow]

    tmp.firstClass.PD = cv.corr.firstClass$max.PD[iRow]
    tmp.firstClass.PR = cv.corr.firstClass$max.PR[iRow]

    tmp.minPbestClass.PD = cv.corr.bestClass$minP.PD[iRow]
    tmp.minPbestClass.PR = cv.corr.bestClass$minP.PR[iRow]

    tmp.minPfirstClass.PD = cv.corr.firstClass$minP.PD[iRow]
    tmp.minPfirstClass.PR = cv.corr.firstClass$minP.PR[iRow]

    tmp.logBestClass.PD = cv.logrank.bestClass$min.PD[iRow]
    tmp.logBestClass.PR = cv.logrank.bestClass$min.PR[iRow]

    tmp.logFirstClass.PD = cv.logrank.firstClass$min.PD[iRow]
    tmp.logFirstClass.PR = cv.logrank.firstClass$min.PR[iRow]

    # this is the "find the RECIST score" stuff from funcPerfMetrics.
    tmpPatRecord = simLinRECIST[[as.character(iPat)]]
    tmpNumUpdates = length(tmpPatRecord$changeBaseline)
    tmpPatRecord$bestCuts.classRECIST = NA
    tmpPatRecord$firstCuts.classRECIST = NA
    tmpPatRecord$minPbestCuts.classRECIST = NA
    tmpPatRecord$minPfirstCuts.classRECIST = NA
    tmpPatRecord$logBestCuts.classRECIST = NA
    tmpPatRecord$logFirstCuts.classRECIST = NA

    # Generate the RECIST classifications
    for (iUpdate in 1:tmpNumUpdates) {
        # (1) use the cutpoints from the best class output
        tmpPatRecord$bestCuts.classRECIST[iUpdate] = classRECIST(
            tmpPatRecord$changeBaseline[iUpdate],
            tmpPatRecord$changeNadir[iUpdate],
            c(tmp.bestClass.PR, tmp.bestClass.PD))
        # (2) use the cutpoints from the first class output
        tmpPatRecord$firstCuts.classRECIST[iUpdate] = classRECIST(
            tmpPatRecord$changeBaseline[iUpdate],
            tmpPatRecord$changeNadir[iUpdate],
            c(tmp.firstClass.PR, tmp.firstClass.PD))
        # (3) use the cutpoints from the best class minP output
        tmpPatRecord$minPbestCuts.classRECIST[iUpdate] = classRECIST(
            tmpPatRecord$changeBaseline[iUpdate],
            tmpPatRecord$changeNadir[iUpdate],
            c(tmp.minPbestClass.PR, tmp.minPbestClass.PD))
        # (4) use the cutpoints from the first class minP output
        tmpPatRecord$minPfirstCuts.classRECIST[iUpdate] = classRECIST(
            tmpPatRecord$changeBaseline[iUpdate],
            tmpPatRecord$changeNadir[iUpdate],
            c(tmp.minPfirstClass.PR, tmp.minPfirstClass.PD))
        # (5) use the cutpoints from the logrank best class output
        tmpPatRecord$logBestCuts.classRECIST[iUpdate] = classRECIST(
            tmpPatRecord$changeBaseline[iUpdate],
            tmpPatRecord$changeNadir[iUpdate],
            c(tmp.logBestClass.PR, tmp.logBestClass.PD))
        # (6) use the cutpoints from the logrank first class output
        tmpPatRecord$logFirstCuts.classRECIST[iUpdate] = classRECIST(
            tmpPatRecord$changeBaseline[iUpdate],
            tmpPatRecord$changeNadir[iUpdate],
            c(tmp.logFirstClass.PR, tmp.logFirstClass.PD))
    }

    # Determine the classifications that "stick"
    tmpPatRecord$bestCuts.classRECIST = factor(tmpPatRecord$bestCuts.classRECIST,
        levels = c('PD', 'SD', 'PR'), labels = c(1, 2, 3))
    tmpPatRecord$firstCuts.classRECIST = factor(tmpPatRecord$firstCuts.classRECIST,
        levels = c('PD', 'SD', 'PR'), labels = c(1, 2, 3))

    tmpPatRecord$minPbestCuts.classRECIST = factor(tmpPatRecord$minPbestCuts.classRECIST,
        levels = c('PD', 'SD', 'PR'), labels = c(1, 2, 3))
    tmpPatRecord$minPfirstCuts.classRECIST = factor(tmpPatRecord$minPfirstCuts.classRECIST,
        levels = c('PD', 'SD', 'PR'), labels = c(1, 2, 3))

    tmpPatRecord$logBestCuts.classRECIST = factor(tmpPatRecord$logBestCuts.classRECIST,
        levels = c('PD', 'SD', 'PR'), labels = c(1, 2, 3))
    tmpPatRecord$logFirstCuts.classRECIST = factor(tmpPatRecord$logFirstCuts.classRECIST,
        levels = c('PD', 'SD', 'PR'), labels = c(1, 2, 3))

    tmpRow = match(as.character(iPat), cv.Class$URN)
    cv.Class$bestCuts.bestClass[tmpRow] = max(as.numeric(tmpPatRecord$bestCuts.classRECIST))
    cv.Class$bestCuts.firstClass[tmpRow] = tmpPatRecord$bestCuts.classRECIST[1]
    cv.Class$firstCuts.bestClass[tmpRow] = max(as.numeric(tmpPatRecord$firstCuts.classRECIST))
    cv.Class$firstCuts.firstClass[tmpRow] = tmpPatRecord$firstCuts.classRECIST[1]

    cv.Class$minPbestCuts.bestClass[tmpRow] = max(as.numeric(tmpPatRecord$minPbestCuts.classRECIST))
    cv.Class$minPbestCuts.firstClass[tmpRow] = tmpPatRecord$minPbestCuts.classRECIST[1]
    cv.Class$minPfirstCuts.bestClass[tmpRow] = max(as.numeric(tmpPatRecord$minPfirstCuts.classRECIST))
    cv.Class$minPfirstCuts.firstClass[tmpRow] = tmpPatRecord$minPfirstCuts.classRECIST[1]

    cv.Class$logBestCuts.bestClass[tmpRow] = max(as.numeric(tmpPatRecord$logBestCuts.classRECIST))
    cv.Class$logBestCuts.firstClass[tmpRow] = tmpPatRecord$logBestCuts.classRECIST[1]
    cv.Class$logFirstCuts.bestClass[tmpRow] = max(as.numeric(tmpPatRecord$logFirstCuts.classRECIST))
    cv.Class$logFirstCuts.firstClass[tmpRow] = tmpPatRecord$logFirstCuts.classRECIST[1]
}

# now, to see how they all do.
cv.perf.bestCuts.bestClass = rcorr.cens(cv.Class$bestCuts.bestClass,
    Surv(cv.Class$SurvTimeDays, cv.Class$Status), outx = TRUE)
cv.perf.bestCuts.firstClass = rcorr.cens(cv.Class$bestCuts.firstClass,
    Surv(cv.Class$SurvTimeDays, cv.Class$Status), outx = TRUE)

cv.perf.firstCuts.bestClass = rcorr.cens(cv.Class$firstCuts.bestClass,
    Surv(cv.Class$SurvTimeDays, cv.Class$Status), outx = TRUE)
cv.perf.firstCuts.firstClass = rcorr.cens(cv.Class$firstCuts.firstClass,
    Surv(cv.Class$SurvTimeDays, cv.Class$Status), outx = TRUE)

cv.perf.minPbestCuts.bestClass = rcorr.cens(cv.Class$minPbestCuts.bestClass,
    Surv(cv.Class$SurvTimeDays, cv.Class$Status), outx = TRUE)
cv.perf.minPbestCuts.firstClass = rcorr.cens(cv.Class$minPbestCuts.firstClass,
    Surv(cv.Class$SurvTimeDays, cv.Class$Status), outx = TRUE)

cv.perf.minPfirstCuts.bestClass = rcorr.cens(cv.Class$minPfirstCuts.bestClass,
    Surv(cv.Class$SurvTimeDays, cv.Class$Status), outx = TRUE)
cv.perf.minPfirstCuts.firstClass = rcorr.cens(cv.Class$minPfirstCuts.firstClass,
    Surv(cv.Class$SurvTimeDays, cv.Class$Status), outx = TRUE)

cv.log.bestCuts.bestClass = survdiff(Surv(SurvTimeDays, Status) ~ logBestCuts.bestClass, data = cv.Class)
cv.log.bestCuts.firstClass = survdiff(Surv(SurvTimeDays, Status) ~ logBestCuts.firstClass, data = cv.Class)

cv.log.firstCuts.bestClass = survdiff(Surv(SurvTimeDays, Status) ~ logFirstCuts.bestClass, data = cv.Class)
cv.log.firstCuts.firstClass = survdiff(Surv(SurvTimeDays, Status) ~ logFirstCuts.firstClass, data = cv.Class)

save(list = c(
    ls(pattern = 'cv.+')),
    file = file.path(tmpDataDir, 'simCVDataSetPerformance.rdata'))
rm(list = ls(pattern = 'tmp.*'))
rm(list = ls(pattern = 'i[[:upper:]]'))