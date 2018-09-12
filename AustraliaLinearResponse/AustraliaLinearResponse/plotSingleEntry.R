# Script to analyze the RECIST/Survival correlation data.

## the first step is to generate the classifications according to the new
# set(s) of criteria.
tmp.new.PR=-0.64 #-0.68 happens too.
tmp.new.PD=0.5 #0.26 happens too.

newCut = dataFullCols[,c(1,30,81,82,128)]
newCut$bestClass=NA
newCut$firstClass=NA
newCut$stdBestClass=NA
newCut$stdFirstClass=NA

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
	
	tmpRow=match(as.character(iPat),newCut$URN)
	newCut$bestClass[tmpRow]=max(as.numeric(tmpPatRecord$newCuts.classRECIST))
	newCut$firstClass[tmpRow]=tmpPatRecord$newCuts.classRECIST[1]
	newCut$stdBestClass[tmpRow]=max(as.numeric(tmpPatRecord$stdCuts.classRECIST))
	newCut$stdFirstClass[tmpRow]=tmpPatRecord$stdCuts.classRECIST[1]
}

## Make a clean plot
tmp.op=par(no.readonly=TRUE)
tmpYlim=1500
par(mfcol=c(3,2))

## start with the "best response" approach first.
#std (best) RECIST, C Plot
plot(newCut$stdBestClass,newCut$EntryTimeDays,pch=3-2*newCut$Status,col=2-newCut$Status,
	xaxt='n',yaxt='n',xlim=c(1,3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (days)')
tmpCorr=with(newCut,rcorr.cens(stdBestClass,Surv(EntryTimeDays,Status),outx=TRUE))
axis(1,at=1:3,label=c('PD','SD','PR'))
axis(2,at=seq(0,tmpYlim,300))
title(main=bquote(paste('Standard RECIST, ',C==.(round(tmpCorr['C Index'],3))
%+-%.(round(tmpCorr['S.D.']/2,3)))))
box()
lines(x=c(2.3,2.3),y=c(-50,100),col='red')
arrows(2.3,100,2.7,100,col='red',length=0.15)
text(x=2.4,y=200,labels='-30%',col='red')
lines(x=c(1.7,1.7),y=c(-50,100),col='red')
arrows(1.7,100,1.3,100,col='red',length=0.15)
text(x=1.6,y=200,labels='+20%',col='red')

#optimized (best) RECIST, C plot
plot(newCut$bestClass,newCut$EntryTimeDays,pch=3-2*newCut$Status,col=2-newCut$Status,
	xaxt='n',yaxt='n',xlim=c(1,3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (days)')
tmpCorr=with(newCut,rcorr.cens(bestClass,Surv(EntryTimeDays,Status),outx=TRUE))
axis(1,at=1:3,label=c('PD','SD','PR'))
axis(2,at=seq(0,tmpYlim,300))
title(main=bquote(paste('Optimized RECIST, ',C==.(round(tmpCorr['C Index'],3))
%+-%.(round(tmpCorr['S.D.']/2,3)))))
box()
lines(x=c(2.3,2.3),y=c(-50,100),col='red')
arrows(2.3,100,2.7,100,col='red',length=0.15)
text(x=2.4,y=200,labels=bquote(paste(.(100*tmp.new.PR),'%')),col='red')
lines(x=c(1.7,1.7),y=c(-50,100),col='red')
arrows(1.7,100,1.3,100,col='red',length=0.15)
text(x=1.6,y=200,labels=bquote(paste('+',.(100*tmp.new.PD),'%')),col='red')

# cross-validated optimized (best) RECIST, C plot
plot(cv.Class$bestCuts.bestClass,cv.Class$EntryTimeDays,pch=3-2*cv.Class$Status,col=2-cv.Class$Status,
	xaxt='n',yaxt='n',xlim=c(1,3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (days)')
tmpCorr=cv.perf.bestCuts.bestClass
axis(1,at=1:3,label=c('PD','SD','PR'))
axis(2,at=seq(0,tmpYlim,300))
title(main=bquote(paste('Cross-Validated RECIST, ',C==.(round(tmpCorr['C Index'],3))
%+-%.(round(tmpCorr['S.D.']/2,3)))))
box()

# std (best) RECIST, surv. plot
survplot(survfit(Surv(EntryTimeDays,Status)~stdBestClass,data=newCut),
	xlab='Overall Survival from Diagnosis [Days]',xlim=c(0,tmpYlim),time.inc=300,
	conf='none',label.curves=list(keys='lines',labels=c('PD','SD','PR')),
	lty=c('dotted','dashed','solid'),lwd=1)
tmpLR=survdiff(Surv(EntryTimeDays,Status)~stdBestClass,data=newCut)
tmpLRdf=length(tmpLR$obs)-1
tmpLRcs=tmpLR$chisq
tmpLRp=pchisq(tmpLRcs,tmpLRdf,lower.tail=FALSE)
title(main=bquote(paste('Standard RECIST, log-rank ',p==.(round(tmpLRp,4)))))

# optimized (best) RECIST, surv. plot
survplot(survfit(Surv(EntryTimeDays,Status)~bestClass,data=newCut),
	xlab='Overall Survival from Diagnosis [Days]',xlim=c(0,tmpYlim),time.inc=300,
	conf='none',label.curves=list(keys='lines',labels=c('SD','PR')),
	lty=c('dashed','solid'),lwd=1)
tmpLR=survdiff(Surv(EntryTimeDays,Status)~bestClass,data=newCut)
tmpLRdf=length(tmpLR$obs)-1
tmpLRcs=tmpLR$chisq
tmpLRp=pchisq(tmpLRcs,tmpLRdf,lower.tail=FALSE)
title(main=bquote(paste('Optimized RECIST, log-rank ',p==.(round(tmpLRp,4)))))

# cross-validated optimized (best) RECIST, surv. plot
survplot(survfit(Surv(EntryTimeDays,Status)~bestCuts.bestClass,data=cv.Class),
	xlab='Overall Survival from Diagnosis [Days]',xlim=c(0,tmpYlim),time.inc=300,
	conf='none',label.curves=list(keys='lines',labels=c('PD','SD','PR')),
	lty=c('dotted','dashed','solid'),lwd=1)
tmpLR=survdiff(Surv(EntryTimeDays,Status)~bestCuts.bestClass,data=cv.Class)
tmpLRdf=length(tmpLR$obs)-1
tmpLRcs=tmpLR$chisq
tmpLRp=pchisq(tmpLRcs,tmpLRdf,lower.tail=FALSE)
title(main=bquote(paste('Cross-Validated RECIST, log-rank ',p==.(round(tmpLRp,4)))))




## do the whole thing again for the "first change" metric
dev.new()
tmpYlim=1500
par(mfcol=c(3,2))

#std (first) RECIST, C Plot
plot(newCut$stdFirstClass,newCut$EntryTimeDays,pch=3-2*newCut$Status,col=2-newCut$Status,
	xaxt='n',yaxt='n',xlim=c(1,3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (days)')
tmpCorr=with(newCut,rcorr.cens(stdFirstClass,Surv(EntryTimeDays,Status),outx=TRUE))
axis(1,at=1:3,label=c('PD','SD','PR'))
axis(2,at=seq(0,tmpYlim,300))
title(main=bquote(paste('Standard RECIST @ First Follow-up, ',C==.(round(tmpCorr['C Index'],3))
%+-%.(round(tmpCorr['S.D.']/2,3)))))
box()
lines(x=c(2.3,2.3),y=c(-50,100),col='red')
arrows(2.3,100,2.7,100,col='red',length=0.15)
text(x=2.4,y=200,labels='-30%',col='red')
lines(x=c(1.7,1.7),y=c(-50,100),col='red')
arrows(1.7,100,1.3,100,col='red',length=0.15)
text(x=1.6,y=200,labels='+20%',col='red')

#optimized (first) RECIST, C plot
plot(newCut$firstClass,newCut$EntryTimeDays,pch=3-2*newCut$Status,col=2-newCut$Status,
	xaxt='n',yaxt='n',xlim=c(1,3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (days)')
tmpCorr=with(newCut,rcorr.cens(firstClass,Surv(EntryTimeDays,Status),outx=TRUE))
axis(1,at=1:3,label=c('PD','SD','PR'))
axis(2,at=seq(0,tmpYlim,300))
title(main=bquote(paste('Optimized RECIST @ First Follow-up, ',C==.(round(tmpCorr['C Index'],3))
%+-%.(round(tmpCorr['S.D.']/2,3)))))
box()
lines(x=c(2.3,2.3),y=c(-50,100),col='red')
arrows(2.3,100,2.7,100,col='red',length=0.15)
text(x=2.4,y=200,labels=bquote(paste(.(100*tmp.new.PR),'%')),col='red')
lines(x=c(1.7,1.7),y=c(-50,100),col='red')
arrows(1.7,100,1.3,100,col='red',length=0.15)
text(x=1.6,y=200,labels=bquote(paste('+',.(100*tmp.new.PD),'%')),col='red')

# cross-validated optimized (first) RECIST, C plot
plot(cv.Class$bestCuts.firstClass,cv.Class$EntryTimeDays,pch=3-2*cv.Class$Status,col=2-cv.Class$Status,
	xaxt='n',yaxt='n',xlim=c(1,3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (days)')
tmpCorr=cv.perf.bestCuts.firstClass
axis(1,at=1:3,label=c('PD','SD','PR'))
axis(2,at=seq(0,tmpYlim,300))
title(main=bquote(paste('Cross-Validated RECIST @ First Follow-up, ',C==.(round(tmpCorr['C Index'],3))
%+-%.(round(tmpCorr['S.D.']/2,3)))))
box()

# std (first) RECIST, surv. plot
survplot(survfit(Surv(EntryTimeDays,Status)~stdFirstClass,data=newCut),
	xlab='Overall Survival from Diagnosis [Days]',xlim=c(0,tmpYlim),time.inc=300,
	conf='none',label.curves=list(keys='lines',labels=c('PD','SD','PR')),
	lty=c('dotted','dashed','solid'),lwd=1)
tmpLR=survdiff(Surv(EntryTimeDays,Status)~stdFirstClass,data=newCut)
tmpLRdf=length(tmpLR$obs)-1
tmpLRcs=tmpLR$chisq
tmpLRp=pchisq(tmpLRcs,tmpLRdf,lower.tail=FALSE)
title(main=bquote(paste('Standard RECIST @ First Follow-up, log-rank ',p==.(round(tmpLRp,4)))))

# optimized (first) RECIST, surv. plot
survplot(survfit(Surv(EntryTimeDays,Status)~firstClass,data=newCut),
	xlab='Overall Survival from Diagnosis [Days]',xlim=c(0,tmpYlim),time.inc=300,
	conf='none',label.curves=list(keys='lines',labels=c('SD','PR')),
	lty=c('dashed','solid'),lwd=1)
tmpLR=survdiff(Surv(EntryTimeDays,Status)~firstClass,data=newCut)
tmpLRdf=length(tmpLR$obs)-1
tmpLRcs=tmpLR$chisq
tmpLRp=pchisq(tmpLRcs,tmpLRdf,lower.tail=FALSE)
title(main=bquote(paste('Optimized RECIST @ First Follow-up, log-rank ',p==.(round(tmpLRp,4)))))

# cross-validated optimized (first) RECIST, surv. plot
survplot(survfit(Surv(EntryTimeDays,Status)~bestCuts.firstClass,data=cv.Class),
	xlab='Overall Survival from Diagnosis [Days]',xlim=c(0,tmpYlim),time.inc=300,
	conf='none',label.curves=list(keys='lines',labels=c('PD','SD','PR')),
	lty=c('dotted','dashed','solid'),lwd=1)
tmpLR=survdiff(Surv(EntryTimeDays,Status)~bestCuts.firstClass,data=cv.Class)
tmpLRdf=length(tmpLR$obs)-1
tmpLRcs=tmpLR$chisq
tmpLRp=pchisq(tmpLRcs,tmpLRdf,lower.tail=FALSE)
title(main=bquote(paste('Cross-Validated RECIST @ First Follow-up, log-rank ',p==.(round(tmpLRp,4)))))

# clean up the workspace
par(tmp.op)
rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))