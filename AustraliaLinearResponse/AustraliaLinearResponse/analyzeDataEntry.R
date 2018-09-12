# vary the RECIST cutoffs and see what happens.
tmpBaseDir=file.path('F:','MesoResponse')
tmpProjName='AustraliaLinearResponse'
tmpDataDir=file.path(tmpBaseDir,tmpProjName)
tmpScriptDir=file.path(tmpBaseDir,tmpProjName,tmpProjName)

library(Hmisc)

# set up the continuous metrics
dataContinuous = dataFullCols[,c(1,30,81,82,128,135,139)]
dataContinuous$RecistBest=NA
for (iPat in names(dataLinRECIST)) {
	tmpRow=match(as.character(iPat),dataContinuous$URN)
	dataContinuous$RecistBest[tmpRow]=dataLinRECIST[[as.character(iPat)]]$changeMin
}
# calculate some of the continuous metrics, since those are independent of the loops
corr.Continuous=list()
corr.Continuous$bestChange=rcorr.cens(-dataContinuous$RecistBest,
	Surv(dataContinuous$EntryTimeDays,dataContinuous$Status),outx=TRUE)

corr.Continuous$lastChange=rcorr.cens(-dataContinuous$RecistChangeN,
	Surv(dataContinuous$EntryTimeDays,dataContinuous$Status),outx=TRUE)

corr.Continuous$firstChange=rcorr.cens(-dataContinuous$RecistChange2,
	Surv(dataContinuous$EntryTimeDays,dataContinuous$Status),outx=TRUE)

# do the variable cutpoint thing for the whole dataset

corr.Cuts=seq(-1,1,by=0.01)

source(file.path(tmpScriptDir,'funcPerfMetricsEntry.R'))
tmpOutput=perf.metrics.entry(corr.Cuts,dataFullCols,dataLinRECIST,classRECIST)
corr.bestClass=tmpOutput$corr.bestClass
corr.firstClass=tmpOutput$corr.firstClass
corr.lastClass=tmpOutput$corr.lastClass
logrank.bestClass=tmpOutput$logrank.bestClass
logrank.firstClass=tmpOutput$logrank.firstClass
logrank.lastClass=tmpOutput$logrank.lastClass

save(list=c(
	ls(pattern='corr.+'),
	ls(pattern='logrank.+')),
	file=file.path(tmpDataDir,'FullDataSetPerformanceEntry.rdata'))
rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))