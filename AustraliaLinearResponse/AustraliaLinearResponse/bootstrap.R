tmpBaseDir=file.path('F:','MesoResponse')
tmpProjName='AustraliaLinearResponse'
tmpDataDir=file.path(tmpBaseDir,tmpProjName)
tmpScriptDir=file.path(tmpBaseDir,tmpProjName,tmpProjName)
library(Hmisc)

# do the variable cutpoint thing each validation subset
corr.Cuts=seq(-1,1,by=0.01)
tmpCurPR=match(-0.3,round(corr.Cuts,5))
tmpCurPD=match(0.2,round(corr.Cuts,5))
source(file.path(tmpScriptDir,'funcPerfMetrics_noLR.R'))

tmpPatList=names(dataLinRECIST)
boot.bestClass=data.frame(max=NA,max.PR=NA,max.PD=NA,std=NA)
boot.firstClass=data.frame(max=NA,max.PR=NA,max.PD=NA,std=NA)

iBootNum=200

for (iCount in seq(1,iBootNum)) {
	# new patient fake cohort
	tmpBootSamp=sample(1:length(dataLinRECIST),length(dataLinRECIST),replace=TRUE)
	
	# make the cross validation patient cohort (leaving out iPat)
	tmpNewPatList=tmpPatList[tmpBootSamp]
	tmp.dataLinRECIST=dataLinRECIST[match(tmpNewPatList,tmpPatList)]
	# this one is trickier...URNs are not in the same order!
	tmp.dataFullCols=dataFullCols[match(tmpNewPatList,dataFullCols$URN),]
	
	# find the performance for each damn run
	try(expr=(tmpOutput=perf.metrics(corr.Cuts,tmp.dataFullCols,tmp.dataLinRECIST,classRECIST)))
	
	# find the things that matter
	boot.bestClass[iCount,'max']=tmpOutput$corr.bestClass$max
	boot.bestClass[iCount,'max.PR']=tmpOutput$corr.bestClass$max.PR
	boot.bestClass[iCount,'max.PD']=tmpOutput$corr.bestClass$max.PD
	boot.bestClass[iCount,'std']=with(tmpOutput$corr.bestClass,perf[tmpCurPR,tmpCurPD])
	
	boot.firstClass[iCount,'max']=tmpOutput$corr.firstClass$max
	boot.firstClass[iCount,'max.PR']=tmpOutput$corr.firstClass$max.PR
	boot.firstClass[iCount,'max.PD']=tmpOutput$corr.firstClass$max.PD
	boot.firstClass[iCount,'std']=with(tmpOutput$corr.firstClass,perf[tmpCurPR,tmpCurPD])
}

save(list=c(
	ls(pattern='boot.+')),
	file=file.path(tmpDataDir,'BootDataSetPerformance.rdata'))
rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))