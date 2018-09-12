tmpBaseDir=file.path('F:','MesoResponse')
tmpProjName='AustraliaLinearResponse'
tmpDataDir=file.path(tmpBaseDir,tmpProjName)
tmpScriptDir=file.path(tmpBaseDir,tmpProjName,tmpProjName)
library(Hmisc)
source(file.path(tmpScriptDir,'funcPerfMetrics.R'))
## --------------------------- ##


jk.best.corrOptC=with(cv.corr.bestClass,cor(max,std,method='pearson'))
jk.first.corrOptC=with(cv.corr.firstClass,cor(max,std,method='pearson'))

tmp.new.PR=-0.64 #-0.68 happens too.
tmp.new.PD=0.5 #0.26 happens too.

newCut = dataFullCols[,c(1,30,81,82,83,128)]
newCut$bestClass=NA
newCut$firstClass=NA
newCut$stdBestClass=NA
newCut$stdFirstClass=NA
newCut$std.scanOfPD=NA
newCut$new.scanOfPD=NA
newCut$std.timeToPD=NA
newCut$new.timeToPD=NA
newCut$worstClass=NA
newCut$stdWorstClass=NA

for (iPat in names(dataLinRECIST)) {
	
	# this is the "find the RECIST score" stuff from funcPerfMetrics.
	tmpPatRecord = dataLinRECIST[[as.character(iPat)]]
	tmpNumUpdates=length(tmpPatRecord$changeBaseline)
	tmpPatRecord$newCuts.classRECIST=NA
	tmpPatRecord$stdCuts.classRECIST=NA
	
	# Generate the RECIST classifications
	for (iUpdate in 1:tmpNumUpdates) {
		# (1) use the cutpoints from the new (supplied) numbers
		tmpPatRecord$newCuts.classRECIST[iUpdate]=classRECIST(
			tmpPatRecord$changeBaseline[iUpdate],
			tmpPatRecord$changeNadir[iUpdate],
			c(tmp.new.PR,tmp.new.PD))
		# (2) use the cutpoints from the standard algorithm
		tmpPatRecord$stdCuts.classRECIST[iUpdate]=classRECIST(
			tmpPatRecord$changeBaseline[iUpdate],
			tmpPatRecord$changeNadir[iUpdate],
			c(-0.3,0.2))
	}
	
	# Determine the classifications that "stick"
	tmpPatRecord$newCuts.classRECIST=factor(tmpPatRecord$newCuts.classRECIST,
		levels=c('PD','SD','PR'),labels=c(1,2,3))
	tmpPatRecord$stdCuts.classRECIST=factor(tmpPatRecord$stdCuts.classRECIST,
		levels=c('PD','SD','PR'),labels=c(1,2,3))
	
	# find out some time-to-progression stuff
	tmpPatRecord$std.scanOfPD=match(1,tmpPatRecord$stdCuts.classRECIST)+1
	tmpPatRecord$new.scanOfPD=match(1,tmpPatRecord$newCuts.classRECIST)+1
	tmpPatRecord$std.timeToPD=tmpPatRecord$scanDates[tmpPatRecord$std.scanOfPD]
	tmpPatRecord$new.timeToPD=tmpPatRecord$scanDates[tmpPatRecord$new.scanOfPD]
	
	tmpRow=match(as.character(iPat),newCut$URN)
	newCut$bestClass[tmpRow]=max(as.numeric(tmpPatRecord$newCuts.classRECIST))
	newCut$firstClass[tmpRow]=tmpPatRecord$newCuts.classRECIST[1]
	newCut$stdBestClass[tmpRow]=max(as.numeric(tmpPatRecord$stdCuts.classRECIST))
	newCut$stdFirstClass[tmpRow]=tmpPatRecord$stdCuts.classRECIST[1]
	
	newCut$worstClass[tmpRow]=min(as.numeric(tmpPatRecord$newCuts.classRECIST))
	newCut$stdWorstClass[tmpRow]=min(as.numeric(tmpPatRecord$stdCuts.classRECIST))
	
	newCut$std.scanOfPD[tmpRow]=tmpPatRecord$std.scanOfPD
	newCut$new.scanOfPD[tmpRow]=tmpPatRecord$new.scanOfPD
	newCut$std.timeToPD[tmpRow]=
		as.numeric(tmpPatRecord$std.timeToPD-tmpPatRecord$scanDates[1])
	newCut$new.timeToPD[tmpRow]=
		as.numeric(tmpPatRecord$new.timeToPD-tmpPatRecord$scanDates[1])
	
}
rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))

cv.Class$fullStd.bestClass = newCut$stdBestClass
cv.Class$fullStd.firstClass = newCut$stdFirstClass

jk.best.corrOptCuts=with(cv.Class,
	cor(fullStd.bestClass,bestCuts.bestClass,method='kendall'))

jk.first.corrOptCuts=with(cv.Class,
	cor(fullStd.firstClass,bestCuts.firstClass,method='kendall'))

## ----- try some jackknife stuff, motivated by MRMC Dorfman, Metz 1992

jk.bestClass=cv.corr.bestClass[,c('patOut','std')]
jk.firstClass=cv.corr.firstClass[,c('patOut','std')]
jk.bestClass$opt=NA
jk.firstClass$opt=NA
tmpPatList=jk.bestClass$patOut
tmpCuts=c(-0.64,0.5)
for (iPat in tmpPatList) {
	iCount=match(iPat,tmpPatList)
	
	# make the cross validation patient cohort (leaving out iPat)
	tmpNewPatList=subset(tmpPatList,!(tmpPatList %in% iPat))
	tmp.dataFullCols=subset(dataFullCols,(URN %in% tmpNewPatList))
	tmp.dataLinRECIST=dataLinRECIST[match(tmpNewPatList,tmpPatList)]
	
	# find the performance for each damn run
	tmpOutput=perf.metrics(tmpCuts,tmp.dataFullCols,tmp.dataLinRECIST,classRECIST)
	jk.bestClass[iCount,'opt']=tmpOutput$corr.bestClass$max
	jk.firstClass[iCount,'opt']=tmpOutput$corr.firstClass$max
}
rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))

tmpN=length(jk.bestClass$patOut)
jk.bestClass$pseudo.opt=tmpN*corr.bestClass$max - (tmpN-1)*jk.bestClass$opt
jk.bestClass$pseudo.std=tmpN*corr.bestClass$perf[71,121] - (tmpN-1)*jk.bestClass$std

jk.firstClass$pseudo.opt=tmpN*corr.firstClass$max - (tmpN-1)*jk.firstClass$opt
jk.firstClass$pseudo.std=tmpN*corr.firstClass$perf[71,121] - (tmpN-1)*jk.firstClass$std

jk=data.frame(
	patOut=rep(jk.bestClass$patOut,2),
	pseudo.best=c(jk.bestClass$pseudo.std,jk.bestClass$pseudo.opt),
	type=as.factor(c(rep('std',tmpN),rep('opt',tmpN))),
	real.best=c(jk.bestClass$std,jk.bestClass$opt),
	pseudo.first=c(jk.firstClass$pseudo.std,jk.firstClass$pseudo.opt),
	real.first=c(jk.firstClass$std,jk.firstClass$opt))

library(lme4)
jk.bestModel=lmer(pseudo.best~0+type+(1|patOut),data=jk)
jk.firstModel=lmer(pseudo.first~0+type+(1|patOut),data=jk)

## ------------ calculate the sigmas for differences.
tmpInd.PR=match(-0.64,round(corr.Cuts,5))
tmpInd.PD=match(0.50,round(corr.Cuts,5))

sigma.jkModel.best=sqrt(2/tmpN)*attr(VarCorr(jk.bestModel),'sc')
sigma.jkModel.first=sqrt(2/tmpN)*attr(VarCorr(jk.firstModel),'sc')

tmpS.bestOpt=corr.bestClass$perf.stdev[tmpInd.PR,tmpInd.PD]
tmpS.firstOpt=corr.firstClass$perf.stdev[tmpInd.PR,tmpInd.PD]
tmpS.bestStd=corr.bestClass$perf.stdev[71,121]
tmpS.firstStd=corr.firstClass$perf.stdev[71,121]

sigma.indy.best=sqrt(tmpS.bestOpt^2+tmpS.bestStd^2)
sigma.indy.first=sqrt(tmpS.firstOpt^2+tmpS.firstStd^2)

sigma.corrC.best=sqrt(tmpS.bestOpt^2+tmpS.bestStd^2-2*tmpS.bestOpt*tmpS.bestStd*jk.best.corrOptC)
sigma.corrC.first=sqrt(tmpS.firstOpt^2+tmpS.firstStd^2-2*tmpS.firstOpt*tmpS.firstStd*jk.first.corrOptC)

sigma.corrCuts.best=sqrt(tmpS.bestOpt^2+tmpS.bestStd^2-2*tmpS.bestOpt*tmpS.bestStd*jk.best.corrOptCuts)
sigma.corrCuts.first=sqrt(tmpS.firstOpt^2+tmpS.firstStd^2-2*tmpS.firstOpt*tmpS.firstStd*jk.first.corrOptCuts)

rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))