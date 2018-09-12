tmpBaseDir=file.path('F:','MesoResponse')
tmpProjName='AustraliaLinearResponse'
tmpDataDir=file.path(tmpBaseDir,tmpProjName)
tmpScriptDir=file.path(tmpBaseDir,tmpProjName,tmpProjName)

library(Hmisc)

tmpN=length(corr.Cuts)
tmpArray=matrix(nrow=1,ncol=tmpN,data=0)

corr2.bestClass=list(perf=tmpArray,perf.sd=tmpArray)
corr2.firstClass=list(perf=tmpArray,perf.sd=tmpArray)

tmpStarterData=dataFullCols[,c(1,30,81,83)]
tmpListPats=names(dataLinRECIST)

for (iCol in 1:tmpN) {
	iCut=corr.Cuts[iCol]
	tmpPatClass = tmpStarterData
	tmpPatClass$bestClass=NA
	tmpPatClass$lastClass=NA
	tmpPatClass$firstClass=NA
	
	for (iPat in tmpListPats) {

		tmpPatRecord = dataLinRECIST[[as.character(iPat)]]
		tmpNumUpdates=length(tmpPatRecord$changeBaseline)
		tmpPatRecord$classRECIST=NA
		
		# Generate the RECIST classifications
		for (iUpdate in 1:tmpNumUpdates) {
			tmpPatRecord$classRECIST[iUpdate]=classRECIST(
				tmpPatRecord$changeBaseline[iUpdate],
				tmpPatRecord$changeNadir[iUpdate],
				c(iCut,100))
		}
		
		# Determine the classifications that "stick"
		tmpPatRecord$classRECIST=factor(tmpPatRecord$classRECIST,
			levels=c('PD','SD','PR'),labels=c(1,2,3))
		tmpRow=match(as.character(iPat),tmpPatClass$URN)
		tmpPatClass$bestClass[tmpRow]=max(as.numeric(tmpPatRecord$classRECIST))
		tmpPatClass$lastClass[tmpRow]=tmpPatRecord$classRECIST[tmpNumUpdates]
		tmpPatClass$firstClass[tmpRow]=tmpPatRecord$classRECIST[1]
	}
	
	# calculate the model performance at this cutpoint
	tmpCorrObj1=rcorr.cens(tmpPatClass$bestClass,
		Surv(tmpPatClass$SurvTimeDays,tmpPatClass$Status),outx=TRUE)
	
	tmpCorrObj3=rcorr.cens(tmpPatClass$firstClass,
		Surv(tmpPatClass$SurvTimeDays,tmpPatClass$Status),outx=TRUE)
	
	corr2.bestClass$perf[iCol]=tmpCorrObj1['C Index']
	corr2.bestClass$perf.sd[iCol]=tmpCorrObj1['S.D.']/2
	
	corr2.firstClass$perf[iCol]=tmpCorrObj3['C Index']
	corr2.firstClass$perf.sd[iCol]=tmpCorrObj3['S.D.']/2
}

corr2.bestClass$max=max(corr2.bestClass$perf,na.rm=TRUE)
corr2.bestClass$max.PR=corr.Cuts[arrayInd(which.max(corr2.bestClass$perf),dim(corr2.bestClass$perf))[2]]

corr2.firstClass$max=max(corr2.firstClass$perf,na.rm=TRUE)
corr2.firstClass$max.PR=corr.Cuts[arrayInd(which.max(corr2.firstClass$perf),dim(corr2.firstClass$perf))[2]]

rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))

## ------- Plotting
plot(corr.Cuts,corr2.bestClass$perf,
	type='l',lty='solid',lwd=2,
	xlab='Classification Cut-point',
	ylab='Classification Performance, C',
	xlim=c(-1,0),ylim=c(0,1))

lines(corr.Cuts,corr2.bestClass$perf+1.96*corr2.bestClass$perf.sd,
	lty='dashed',lwd=1)

lines(corr.Cuts,corr2.bestClass$perf-1.96*corr2.bestClass$perf.sd,
	lty='dashed',lwd=1)

## --------- same plot, first cuts
plot(corr.Cuts,corr2.firstClass$perf,
	type='l',lty='solid',lwd=2,
	xlab='Classification Cut-point',
	ylab='Classification Performance, C',
	xlim=c(-1,0),ylim=c(0,1))

lines(corr.Cuts,corr2.firstClass$perf+1.96*corr2.firstClass$perf.sd,
	lty='dashed',lwd=1)

lines(corr.Cuts,corr2.firstClass$perf-1.96*corr2.firstClass$perf.sd,
	lty='dashed',lwd=1)