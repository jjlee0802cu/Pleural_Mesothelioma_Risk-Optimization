tmpYlim=64
tmpXlim=64

dev.new(width=7.5,height=5.5)
survplot(survfit(Surv(SurvTimeMonths,Status)~1,data=dataFullCols),
	xlab='',
	ylab='Survival Probability',
	time.inc=12,
	conf='none',
	lwd=2,
	n.risk=TRUE,
	y.n.risk=-0.2,
	cex.n.risk=0.75,
	xlim=c(0,tmpXlim))

mtext('Months',side=1,line=0.6,at=c(-3),adj=1,cex=0.8)
mtext('Number at risk',side=1,line=2.05,at=c(-3),adj=1,cex=0.65)

text(18,0.8,
	labels='n = 78 (75 deaths)\nMedian Survival = 14.9 Months, 95% CI (12.5-17.0)',
	adj=0)

## ---------- some of the new cutpoint plots -------- ##

dev.new(width=6,height=5)
## ---------------- std (best) RECIST, C Plot
plot(jitter(newCut$stdBestClass,factor=0.25),newCut$SurvTimeMonths,
	pch=1+4*newCut$Status,col=2-newCut$Status,
	cex=2,lwd=2,cex.lab=1.2,
	xaxt='n',yaxt='n',xlim=c(0.7,3.3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (months)')
axis(1,at=1:3,label=c('PD\n(+20%)','SD\n','PR\n(-30%)'),tick=FALSE)
axis(2,at=seq(0,tmpYlim,12))
box()

dev.new(width=6,height=5)
## ------------------ optimized (best) RECIST, C plot
plot(jitter(newCut$bestClass,factor=0.25),newCut$SurvTimeMonths,
	pch=1+4*newCut$Status,col=2-newCut$Status,
	cex=2,lwd=2,cex.lab=1.2,
	xaxt='n',yaxt='n',xlim=c(0.7,3.3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (months)')
axis(1,at=1:3,label=c('PD\n(+50%)','SD\n','PR\n(-64%)'),tick=FALSE)
axis(2,at=seq(0,tmpYlim,12))
box()

dev.new(width=7.5,height=5.5)
## ------------------- std (best) RECIST, surv. plot
survplot(survfit(Surv(SurvTimeMonths,Status)~stdBestClass,data=newCut),
	xlab='',
	xlim=c(0,tmpXlim),time.inc=12,
	ylab='Survival Probability',
	conf='none',
	label.curves=list(keys='lines',labels=c('PD, n=2','SD, n=50','PR, n=26')),
	lty=c('dotted','dashed','solid'),lwd=2,
	col=c('red','black','goldenrod'))
mtext('Months',side=1,line=0.6,at=c(-3),adj=1,cex=0.8)

dev.new(width=7.5,height=5.5)
## ------------------- optimized (best) RECIST, surv. plot
survplot(survfit(Surv(SurvTimeMonths,Status)~bestClass,data=newCut),
	xlab='',
	xlim=c(0,tmpXlim),time.inc=12,
	ylab='Survival Probability',
	conf='none',
	label.curves=list(keys='lines',labels=c('SD, n=67','PR, n=11')),
	lty=c('dashed','solid'),lwd=2,
	col=c('black','goldenrod'))
mtext('Months',side=1,line=0.6,at=c(-3),adj=1,cex=0.8)

dev.new(width=6,height=5)
## ---------------- std (first) RECIST, C Plot
plot(jitter(newCut$stdFirstClass,factor=0.25),newCut$SurvTimeMonths,
	pch=1+4*newCut$Status,col=2-newCut$Status,
	cex=2,lwd=2,cex.lab=1.2,
	xaxt='n',yaxt='n',xlim=c(0.7,3.3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (months)')
axis(1,at=1:3,label=c('PD\n(+20%)','SD\n','PR\n(-30%)'),tick=FALSE)
axis(2,at=seq(0,tmpYlim,12))
box()

dev.new(width=6,height=5)
## ------------------ optimized (first) RECIST, C plot
plot(jitter(newCut$firstClass,factor=0.25),newCut$SurvTimeMonths,
	pch=1+4*newCut$Status,col=2-newCut$Status,
	cex=2,lwd=2,cex.lab=1.2,
	xaxt='n',yaxt='n',xlim=c(0.7,3.3),ylim=c(0,tmpYlim),
	xlab='Response Class',ylab='Survival (months)')
axis(1,at=1:3,label=c('PD\n(+50%)','SD\n','PR\n(-64%)'),tick=FALSE)
axis(2,at=seq(0,tmpYlim,12))
box()