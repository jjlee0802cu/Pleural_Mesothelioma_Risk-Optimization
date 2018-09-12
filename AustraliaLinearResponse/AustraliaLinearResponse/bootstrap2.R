tmpBaseDir=file.path('F:','MesoResponse')
tmpProjName='AustraliaLinearResponse'
tmpDataDir=file.path(tmpBaseDir,tmpProjName)
tmpScriptDir=file.path(tmpBaseDir,tmpProjName,tmpProjName)
library(Hmisc)

# do the variable cutpoint thing each validation subset
corr.Cuts=c(-0.64,-0.3,0.2,0.5) # changed from before
tmpCurPR=match(-0.3,round(corr.Cuts,5))
tmpCurPD=match(0.2,round(corr.Cuts,5))
tmpNewPR=match(-0.64,round(corr.Cuts,5))
tmpNewPD=match(0.50,round(corr.Cuts,5))
source(file.path(tmpScriptDir,'funcPerfMetrics_noLR.R'))

tmpPatList=names(dataLinRECIST)
boot.bestClass=data.frame(new=NA,std=NA)
boot.firstClass=data.frame(new=NA,std=NA)

iBootNum=500
iBootSampN=round(0.75*length(dataLinRECIST))

for (iCount in seq(1,iBootNum)) {
	# new patient fake cohort
	tmpBootSamp=sample(1:length(dataLinRECIST),iBootSampN,replace=TRUE)
	
	# make the cross validation patient cohort (leaving out iPat)
	tmpNewPatList=tmpPatList[tmpBootSamp]
	tmp.dataLinRECIST=dataLinRECIST[match(tmpNewPatList,tmpPatList)]
	# this one is trickier...URNs are not in the same order!
	tmp.dataFullCols=dataFullCols[match(tmpNewPatList,dataFullCols$URN),]
	
	# find the performance for each damn run
	try(expr=(tmpOutput=perf.metrics(corr.Cuts,tmp.dataFullCols,tmp.dataLinRECIST,classRECIST)))
	
	# find the things that matter
	boot.bestClass[iCount,'new']=with(tmpOutput$corr.bestClass,perf[tmpNewPR,tmpNewPD])
	boot.bestClass[iCount,'std']=with(tmpOutput$corr.bestClass,perf[tmpCurPR,tmpCurPD])
	
	boot.firstClass[iCount,'new']=with(tmpOutput$corr.firstClass,perf[tmpNewPR,tmpNewPD])
	boot.firstClass[iCount,'std']=with(tmpOutput$corr.firstClass,perf[tmpCurPR,tmpCurPD])
}

save(list=c(
	ls(pattern='boot.+')),
	file=file.path(tmpDataDir,'Boot2DataSetPerformance.rdata'))
rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))