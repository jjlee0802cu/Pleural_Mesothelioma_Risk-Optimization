perf.metrics.entry=function(corr.Cuts,dataFullCols,dataLinRECIST,classRECIST) {

library(Hmisc)
library(doSMP)
tmpWorkerList=startWorkers(4)
registerDoSMP(tmpWorkerList)

tmpN=length(corr.Cuts)
tmpMatrix=matrix(nrow=tmpN,ncol=tmpN,data=0)

# perf. is the matrix of P_k for varying RECIST cutoff points.  If the PD,PR
# point is invalid (i.e., PD<PR), then just set P_k = 0.  The up and lo.95
# matrices are the upper and lower 95% confidence intervals on the stat.
# That multiplicative factor is 1.95996
corr.bestClass=list(perf=tmpMatrix,up95=tmpMatrix,lo95=tmpMatrix)
# that first data frame is for doing the STANDARD recist algorithm.
# let's also see how this performs using the last continuous change per patient,
# and the first change per patient (from baseline to very first followup)
corr.lastClass=list(perf=tmpMatrix,up95=tmpMatrix,lo95=tmpMatrix)
corr.firstClass=list(perf=tmpMatrix,up95=tmpMatrix,lo95=tmpMatrix)
# prepare the log-rank matrices too.
logrank.bestClass=list(chisq=tmpMatrix,df=tmpMatrix,p=tmpMatrix+1)
logrank.firstClass=list(chisq=tmpMatrix,df=tmpMatrix,p=tmpMatrix+1)
logrank.lastClass=list(chisq=tmpMatrix,df=tmpMatrix,p=tmpMatrix+1)

# get the patients in this collection
tmpListPats=names(dataLinRECIST)

# make the starter data for the rest of these runs
tmpStarterData=dataFullCols[,c(1,30,81,82,128)]

# calculate the classification tasks for the array of cutpoints.
# this is the part that can be done in parallel.
# the net output of this section of code is the corr.bestClass,
# corr.firstClass, and corr.lastClass variables.  That's it.
tmpOpts=list(chunkSize=ceil(tmpN/getDoParWorkers()))

tmpParOutput=foreach(iRow=1:tmpN, .combine='rbind', .inorder=FALSE, .packages=c('surviva','Hmisc'), .options.smp=tmpOpts) %:%
	foreach(iCol=1:tmpN, .combine='rbind', .inorder=FALSE, .packages=c('survival','Hmisc')) %dopar% {
		tmpCutPR = corr.Cuts[iRow]	# row tied to PR
		tmpCutPD = corr.Cuts[iCol]	# column tied to PD
		
		if (tmpCutPD>=0 && tmpCutPR<=0) { # otherwise it doesn't make sense.
		
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
					c(tmpCutPR,tmpCutPD))
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
			Surv(tmpPatClass$EntryTimeDays,tmpPatClass$Status),outx=TRUE)
		
		tmpCorrObj2=rcorr.cens(tmpPatClass$lastClass,
			Surv(tmpPatClass$EntryTimeDays,tmpPatClass$Status),outx=TRUE)
		
		tmpCorrObj3=rcorr.cens(tmpPatClass$firstClass,
			Surv(tmpPatClass$EntryTimeDays,tmpPatClass$Status),outx=TRUE)
		
		# calculate the log-rank tests at this cutpoint
		if (length(unique(tmpPatClass$bestClass))>1) {
			tmpLogRank1=survdiff(Surv(EntryTimeDays,Status)~bestClass,data=tmpPatClass)
			tmpLRdf1=length(tmpLogRank1$obs)-1
			tmpLRcs1=tmpLogRank1$chisq
			tmpLRp1=pchisq(tmpLRcs1,tmpLRdf1,lower.tail=FALSE)
		}
		else {
			tmpLRdf1=0; tmpLRcs1=0; tmpLRp1=1;
		}
		
		if (length(unique(tmpPatClass$lastClass))>1) {
			tmpLogRank2=survdiff(Surv(EntryTimeDays,Status)~lastClass,data=tmpPatClass)
			tmpLRdf2=length(tmpLogRank2$obs)-1
			tmpLRcs2=tmpLogRank2$chisq
			tmpLRp2=pchisq(tmpLRcs2,tmpLRdf2,lower.tail=FALSE)
		}
		else {
			tmpLRdf2=0; tmpLRcs2=0; tmpLRp2=1;
		}
		
		if (length(unique(tmpPatClass$firstClass))>1) {
			tmpLogRank3=survdiff(Surv(EntryTimeDays,Status)~firstClass,data=tmpPatClass)
			tmpLRdf3=length(tmpLogRank3$obs)-1
			tmpLRcs3=tmpLogRank3$chisq
			tmpLRp3=pchisq(tmpLRcs3,tmpLRdf3,lower.tail=FALSE)
		}
		else {
			tmpLRdf3=0; tmpLRcs3=0; tmpLRp3=1;
		}

		data.frame(row=iRow,col=iCol,
			bcP=tmpCorrObj1['C Index'],
			bcU=tmpCorrObj1['C Index'] + 1.95996*tmpCorrObj1['S.D.']/2,
			bcL=tmpCorrObj1['C Index'] - 1.95996*tmpCorrObj1['S.D.']/2,
			lcP=tmpCorrObj2['C Index'],
			lcU=tmpCorrObj2['C Index'] + 1.95996*tmpCorrObj2['S.D.']/2,
			lcL=tmpCorrObj2['C Index'] - 1.95996*tmpCorrObj2['S.D.']/2,
			fcP=tmpCorrObj3['C Index'],
			fcU=tmpCorrObj3['C Index'] + 1.95996*tmpCorrObj3['S.D.']/2,
			fcL=tmpCorrObj3['C Index'] - 1.95996*tmpCorrObj3['S.D.']/2,
			lr.BCcs=tmpLRcs1, lr.BCdf=tmpLRdf1, lr.BCp=tmpLRp1,
			lr.LCcs=tmpLRcs2, lr.LCdf=tmpLRdf2, lr.LCp=tmpLRp2,
			lr.FCcs=tmpLRcs3, lr.FCdf=tmpLRdf3, lr.FCp=tmpLRp3)
		
		}
		else {
			data.frame(row=iRow,col=iCol,
			bcP=0,bcU=0,bcL=0,lcP=0,lcU=0,lcL=0,fcP=0,fcU=0,fcL=0,
			lr.BCcs=0,lr.BCdf=0,lr.BCp=1,
			lr.LCcs=0,lr.LCdf=0,lr.LCp=1,
			lr.FCcs=0,lr.FCdf=0,lr.FCp=1)
		}
}

stopWorkers(tmpWorkerList)

tmpParOutput=tmpParOutput[order(tmpParOutput$row,tmpParOutput$col),]
corr.bestClass$perf=matrix(tmpParOutput$bcP,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.bestClass$up95=matrix(tmpParOutput$bcU,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.bestClass$lo95=matrix(tmpParOutput$bcL,nrow=tmpN,ncol=tmpN,byrow=TRUE)

corr.lastClass$perf=matrix(tmpParOutput$lcP,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.lastClass$up95=matrix(tmpParOutput$lcU,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.lastClass$lo95=matrix(tmpParOutput$lcL,nrow=tmpN,ncol=tmpN,byrow=TRUE)

corr.firstClass$perf=matrix(tmpParOutput$fcP,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.firstClass$up95=matrix(tmpParOutput$fcU,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.firstClass$lo95=matrix(tmpParOutput$fcL,nrow=tmpN,ncol=tmpN,byrow=TRUE)

logrank.bestClass$chisq=matrix(tmpParOutput$lr.BCcs,nrow=tmpN,ncol=tmpN,byrow=TRUE)
logrank.bestClass$df=matrix(tmpParOutput$lr.BCdf,nrow=tmpN,ncol=tmpN,byrow=TRUE)
logrank.bestClass$p=matrix(tmpParOutput$lr.BCp,nrow=tmpN,ncol=tmpN,byrow=TRUE)

logrank.lastClass$chisq=matrix(tmpParOutput$lr.LCcs,nrow=tmpN,ncol=tmpN,byrow=TRUE)
logrank.lastClass$df=matrix(tmpParOutput$lr.LCdf,nrow=tmpN,ncol=tmpN,byrow=TRUE)
logrank.lastClass$p=matrix(tmpParOutput$lr.LCp,nrow=tmpN,ncol=tmpN,byrow=TRUE)

logrank.firstClass$chisq=matrix(tmpParOutput$lr.FCcs,nrow=tmpN,ncol=tmpN,byrow=TRUE)
logrank.firstClass$df=matrix(tmpParOutput$lr.FCdf,nrow=tmpN,ncol=tmpN,byrow=TRUE)
logrank.firstClass$p=matrix(tmpParOutput$lr.FCp,nrow=tmpN,ncol=tmpN,byrow=TRUE)

# find the maximal statistics
corr.bestClass$max=max(corr.bestClass$perf,na.rm=TRUE)
corr.bestClass$max.PR=corr.Cuts[arrayInd(which.max(corr.bestClass$perf),dim(corr.bestClass$perf))[1]]
corr.bestClass$max.PD=corr.Cuts[arrayInd(which.max(corr.bestClass$perf),dim(corr.bestClass$perf))[2]]
corr.bestClass$perf.stdev=with(corr.bestClass,(up95-perf)/1.95996)
corr.bestClass$max.indPR=arrayInd(which.max(corr.bestClass$perf),dim(corr.bestClass$perf))[1]
corr.bestClass$max.indPD=arrayInd(which.max(corr.bestClass$perf),dim(corr.bestClass$perf))[2]

corr.lastClass$max=max(corr.lastClass$perf,na.rm=TRUE)
corr.lastClass$max.PR=corr.Cuts[arrayInd(which.max(corr.lastClass$perf),dim(corr.lastClass$perf))[1]]
corr.lastClass$max.PD=corr.Cuts[arrayInd(which.max(corr.lastClass$perf),dim(corr.lastClass$perf))[2]]
corr.lastClass$perf.stdev=with(corr.lastClass,(up95-perf)/1.95996)
corr.lastClass$max.indPR=arrayInd(which.max(corr.lastClass$perf),dim(corr.lastClass$perf))[1]
corr.lastClass$max.indPD=arrayInd(which.max(corr.lastClass$perf),dim(corr.lastClass$perf))[2]

corr.firstClass$max=max(corr.firstClass$perf,na.rm=TRUE)
corr.firstClass$max.PR=corr.Cuts[arrayInd(which.max(corr.firstClass$perf),dim(corr.firstClass$perf))[1]]
corr.firstClass$max.PD=corr.Cuts[arrayInd(which.max(corr.firstClass$perf),dim(corr.firstClass$perf))[2]]
corr.firstClass$perf.stdev=with(corr.firstClass,(up95-perf)/1.95996)
corr.firstClass$max.indPR=arrayInd(which.max(corr.firstClass$perf),dim(corr.firstClass$perf))[1]
corr.firstClass$max.indPD=arrayInd(which.max(corr.firstClass$perf),dim(corr.firstClass$perf))[2]

logrank.bestClass$min=min(logrank.bestClass$p,na.rm=TRUE)
logrank.bestClass$min.PR=corr.Cuts[arrayInd(which.min(logrank.bestClass$p),dim(logrank.bestClass$p))[1]]
logrank.bestClass$min.PD=corr.Cuts[arrayInd(which.min(logrank.bestClass$p),dim(logrank.bestClass$p))[2]]

logrank.lastClass$min=min(logrank.lastClass$p,na.rm=TRUE)
logrank.lastClass$min.PR=corr.Cuts[arrayInd(which.min(logrank.lastClass$p),dim(logrank.lastClass$p))[1]]
logrank.lastClass$min.PD=corr.Cuts[arrayInd(which.min(logrank.lastClass$p),dim(logrank.lastClass$p))[2]]

logrank.firstClass$min=min(logrank.firstClass$p,na.rm=TRUE)
logrank.firstClass$min.PR=corr.Cuts[arrayInd(which.min(logrank.firstClass$p),dim(logrank.firstClass$p))[1]]
logrank.firstClass$min.PD=corr.Cuts[arrayInd(which.min(logrank.firstClass$p),dim(logrank.firstClass$p))[2]]

# calculate significance matrices for each result against the std model.
tmpCurPR=match(-0.3,round(corr.Cuts,5))
tmpCurPD=match(0.2,round(corr.Cuts,5))

tmpData=corr.bestClass
tmpStdPerf=with(tmpData,(perf[tmpCurPR,tmpCurPD]))
tmpStdStdev=with(tmpData,(perf.stdev[tmpCurPR,tmpCurPD]))
tmpDiff=tmpData$perf-tmpStdPerf
tmpStdev=sqrt(tmpData$perf.stdev^2+tmpStdStdev^2)
tmpZ=tmpDiff/tmpStdev
tmpZ[is.nan(tmpZ)]=0
tmpZ[tmpData$perf==0]=0
corr.bestClass$perf.p=pnorm(tmpZ,lower.tail=FALSE)

tmpData=corr.firstClass
tmpStdPerf=with(tmpData,(perf[tmpCurPR,tmpCurPD]))
tmpStdStdev=with(tmpData,(perf.stdev[tmpCurPR,tmpCurPD]))
tmpDiff=tmpData$perf-tmpStdPerf
tmpStdev=sqrt(tmpData$perf.stdev^2+tmpStdStdev^2)
tmpZ=tmpDiff/tmpStdev
tmpZ[is.nan(tmpZ)]=0
tmpZ[tmpData$perf==0]=0
corr.firstClass$perf.p=pnorm(tmpZ,lower.tail=FALSE)

tmpData=corr.lastClass
tmpStdPerf=with(tmpData,(perf[tmpCurPR,tmpCurPD]))
tmpStdStdev=with(tmpData,(perf.stdev[tmpCurPR,tmpCurPD]))
tmpDiff=tmpData$perf-tmpStdPerf
tmpStdev=sqrt(tmpData$perf.stdev^2+tmpStdStdev^2)
tmpZ=tmpDiff/tmpStdev
tmpZ[is.nan(tmpZ)]=0
tmpZ[tmpData$perf==0]=0
corr.lastClass$perf.p=pnorm(tmpZ,lower.tail=FALSE)

# find the minimum p-value
corr.bestClass$minP=min(corr.bestClass$perf.p,na.rm=TRUE)
corr.bestClass$minP.PR=corr.Cuts[arrayInd(which.min(corr.bestClass$perf.p),dim(corr.bestClass$perf.p))[1]]
corr.bestClass$minP.PD=corr.Cuts[arrayInd(which.min(corr.bestClass$perf.p),dim(corr.bestClass$perf.p))[2]]

corr.lastClass$minP=min(corr.lastClass$perf.p,na.rm=TRUE)
corr.lastClass$minP.PR=corr.Cuts[arrayInd(which.min(corr.lastClass$perf.p),dim(corr.lastClass$perf.p))[1]]
corr.lastClass$minP.PD=corr.Cuts[arrayInd(which.min(corr.lastClass$perf.p),dim(corr.lastClass$perf.p))[2]]

corr.firstClass$minP=min(corr.firstClass$perf.p,na.rm=TRUE)
corr.firstClass$minP.PR=corr.Cuts[arrayInd(which.min(corr.firstClass$perf.p),dim(corr.firstClass$perf.p))[1]]
corr.firstClass$minP.PD=corr.Cuts[arrayInd(which.min(corr.firstClass$perf.p),dim(corr.firstClass$perf.p))[2]]

return(list(
	corr.bestClass=corr.bestClass,
	corr.lastClass=corr.lastClass,
	corr.firstClass=corr.firstClass,
	logrank.bestClass=logrank.bestClass,
	logrank.lastClass=logrank.lastClass,
	logrank.firstClass=logrank.firstClass))

}