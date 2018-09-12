perf.metrics=function(corr.Cuts,dataFullCols,dataLinRECIST,classRECIST) {

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
corr.firstClass=list(perf=tmpMatrix,up95=tmpMatrix,lo95=tmpMatrix)

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
			tmpPatClass$firstClass[tmpRow]=tmpPatRecord$classRECIST[1]
		}
		
		# calculate the model performance at this cutpoint
		tmpCorrObj1=rcorr.cens(tmpPatClass$bestClass,
			Surv(tmpPatClass$SurvTimeDays,tmpPatClass$Status),outx=TRUE)
				
		tmpCorrObj3=rcorr.cens(tmpPatClass$firstClass,
			Surv(tmpPatClass$SurvTimeDays,tmpPatClass$Status),outx=TRUE)
		
		data.frame(row=iRow,col=iCol,
			bcP=tmpCorrObj1['C Index'],
			bcU=tmpCorrObj1['C Index'] + 1.95996*tmpCorrObj1['S.D.']/2,
			bcL=tmpCorrObj1['C Index'] - 1.95996*tmpCorrObj1['S.D.']/2,
			fcP=tmpCorrObj3['C Index'],
			fcU=tmpCorrObj3['C Index'] + 1.95996*tmpCorrObj3['S.D.']/2,
			fcL=tmpCorrObj3['C Index'] - 1.95996*tmpCorrObj3['S.D.']/2)
		
		}
		else {
			data.frame(row=iRow,col=iCol,
			bcP=0,bcU=0,bcL=0,fcP=0,fcU=0,fcL=0)
		}
}

stopWorkers(tmpWorkerList)

tmpParOutput=tmpParOutput[order(tmpParOutput$row,tmpParOutput$col),]
corr.bestClass$perf=matrix(tmpParOutput$bcP,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.bestClass$up95=matrix(tmpParOutput$bcU,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.bestClass$lo95=matrix(tmpParOutput$bcL,nrow=tmpN,ncol=tmpN,byrow=TRUE)

corr.firstClass$perf=matrix(tmpParOutput$fcP,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.firstClass$up95=matrix(tmpParOutput$fcU,nrow=tmpN,ncol=tmpN,byrow=TRUE)
corr.firstClass$lo95=matrix(tmpParOutput$fcL,nrow=tmpN,ncol=tmpN,byrow=TRUE)

# find the maximal statistics
corr.bestClass$max=max(corr.bestClass$perf,na.rm=TRUE)
corr.bestClass$max.PR=corr.Cuts[arrayInd(which.max(corr.bestClass$perf),dim(corr.bestClass$perf))[1]]
corr.bestClass$max.PD=corr.Cuts[arrayInd(which.max(corr.bestClass$perf),dim(corr.bestClass$perf))[2]]
corr.bestClass$perf.stdev=with(corr.bestClass,(up95-perf)/1.95996)
corr.bestClass$max.indPR=arrayInd(which.max(corr.bestClass$perf),dim(corr.bestClass$perf))[1]
corr.bestClass$max.indPD=arrayInd(which.max(corr.bestClass$perf),dim(corr.bestClass$perf))[2]

corr.firstClass$max=max(corr.firstClass$perf,na.rm=TRUE)
corr.firstClass$max.PR=corr.Cuts[arrayInd(which.max(corr.firstClass$perf),dim(corr.firstClass$perf))[1]]
corr.firstClass$max.PD=corr.Cuts[arrayInd(which.max(corr.firstClass$perf),dim(corr.firstClass$perf))[2]]
corr.firstClass$perf.stdev=with(corr.firstClass,(up95-perf)/1.95996)
corr.firstClass$max.indPR=arrayInd(which.max(corr.firstClass$perf),dim(corr.firstClass$perf))[1]
corr.firstClass$max.indPD=arrayInd(which.max(corr.firstClass$perf),dim(corr.firstClass$perf))[2]

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

# find the minimum p-value
corr.bestClass$minP=min(corr.bestClass$perf.p,na.rm=TRUE)
corr.bestClass$minP.PR=corr.Cuts[arrayInd(which.min(corr.bestClass$perf.p),dim(corr.bestClass$perf.p))[1]]
corr.bestClass$minP.PD=corr.Cuts[arrayInd(which.min(corr.bestClass$perf.p),dim(corr.bestClass$perf.p))[2]]

corr.firstClass$minP=min(corr.firstClass$perf.p,na.rm=TRUE)
corr.firstClass$minP.PR=corr.Cuts[arrayInd(which.min(corr.firstClass$perf.p),dim(corr.firstClass$perf.p))[1]]
corr.firstClass$minP.PD=corr.Cuts[arrayInd(which.min(corr.firstClass$perf.p),dim(corr.firstClass$perf.p))[2]]

return(list(
	corr.bestClass=corr.bestClass,
	corr.firstClass=corr.firstClass))

}