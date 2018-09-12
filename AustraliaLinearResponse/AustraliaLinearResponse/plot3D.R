# Plot various crap for the RECIST correlation project
library(rgl)

tmpN=length(corr.Cuts)
tmpFifty = matrix(nrow=tmpN,ncol=tmpN,data=0.5)
tmpCurPR=match(-0.3,round(corr.Cuts,5))
tmpCurPD=match(0.2,round(corr.Cuts,5))
tmpShowLimits=FALSE
tmpLogFlag=FALSE

# three-D surface plots of recist across multiple cut points
open3d()
persp3d(corr.Cuts,corr.Cuts,corr.bestClass$perf,col="slategray",alpha=1,
	xlim=c(-1.0,0.0),ylim=c(0,1),zlim=c(0,1),aspect=c(1,1,0.5),
	add=FALSE,box=TRUE,axes=FALSE,xlab='',ylab='',zlab='')
if (tmpShowLimits) {
#surface3d(corr.Cuts,corr.Cuts,corr.bestClass$up95,col="goldenrod",alpha=0.2)
surface3d(corr.Cuts,corr.Cuts,corr.bestClass$lo95,col="goldenrod",alpha=0.2)
}
persp3d(corr.Cuts,corr.Cuts,tmpFifty,col="tomato",alpha=1,
	add=TRUE,box=FALSE,axes=FALSE,xlab='',ylab='',zlab='')
axes3d(c('x','y','z-+','z++'))
title3d(main="RECIST (Best Response) Model Performance, C metric",
	xlab="PR Cut",ylab="PD Cut")
spheres3d(x=-0.30,y=0.20,z=corr.bestClass$perf[tmpCurPR,tmpCurPD],color='red',rad=0.02)
spheres3d(x=corr.bestClass$max.PR,
	y=corr.bestClass$max.PD,
	z=corr.bestClass$max,color='green',rad=0.02)
#surface3d(corr.Cuts,corr.Cuts,logrank.bestClass$df/2,col='goldenrod',alpha=0.4,back='line')
#surface3d(corr.Cuts,corr.Cuts,1-corr.bestClass$perf.p,col='goldenrod',alpha=0.4)
#surface3d(corr.Cuts,corr.Cuts,tmpFifty+0.45,col="green2",alpha=0.2)

# three-D surface plots of first change across multiple cut points
open3d()
material3d(col='black')
surface3d(corr.Cuts,corr.Cuts,corr.firstClass$perf,col="slategray",alpha=1,back='line')
if (tmpShowLimits) {
#surface3d(corr.Cuts,corr.Cuts,corr.firstClass$up95,col="goldenrod",alpha=0.2)
surface3d(corr.Cuts,corr.Cuts,corr.firstClass$lo95,col="goldenrod",alpha=0.2)
}
surface3d(corr.Cuts,corr.Cuts,tmpFifty,col="tomato",alpha=1)
axes3d(c('x','y','z-+'))
title3d(main="First Change with variable Cut Points",xlab="PR Cut",ylab="PD Cut")
spheres3d(x=-0.30,y=0.20,z=corr.firstClass$perf[tmpCurPR,tmpCurPD],color='red',rad=0.02)
spheres3d(x=corr.firstClass$max.PR,
	y=corr.firstClass$max.PD,
	z=corr.firstClass$max,color='green',rad=0.02)
#surface3d(corr.Cuts,corr.Cuts,logrank.firstClass$df/2,col='goldenrod',alpha=0.4,back='line')
#surface3d(corr.Cuts,corr.Cuts,1-corr.firstClass$perf.p,col='goldenrod',alpha=0.4)
#surface3d(corr.Cuts,corr.Cuts,tmpFifty+0.45,col="green2",alpha=0.2)

if (tmpLogFlag) {
# make some of the logrank plots  - first, the standard algorithm
open3d()
material3d(col='black')
surface3d(corr.Cuts,corr.Cuts,logrank.bestClass$p,zlog=TRUE,col='slategray',back='line')
axes3d(c('x','y','z++'))
title3d(main="Std RECIST log rank p-value",xlab="PR Cut",ylab="PD Cut")
spheres3d(x=-0.30,y=0.20,z=logrank.bestClass$p[tmpCurPR,tmpCurPD],color='red',rad=0.02)
spheres3d(x=logrank.bestClass$min.PR,
	y=logrank.bestClass$min.PD,
	z=logrank.bestClass$min,color='green',rad=0.02)
surface3d(corr.Cuts,corr.Cuts,tmpFifty/10,col="tomato",alpha=0.5)

# then the first change
open3d()
material3d(col='black')
surface3d(corr.Cuts,corr.Cuts,logrank.firstClass$p,zlog=TRUE,col='slategray',back='line')
axes3d(c('x','y','z++'))
title3d(main="First Change log rank p-value",xlab="PR Cut",ylab="PD Cut")
spheres3d(x=-0.30,y=0.20,z=logrank.firstClass$p[tmpCurPR,tmpCurPD],color='red',rad=0.02)
spheres3d(x=logrank.firstClass$min.PR,
	y=logrank.firstClass$min.PD,
	z=logrank.firstClass$min,color='green',rad=0.02)
surface3d(corr.Cuts,corr.Cuts,tmpFifty/10,col="tomato",alpha=0.5)
}


# prettier plot
# three-D surface plots of recist across multiple cut points
open3d()
persp3d(corr.Cuts[1:101],corr.Cuts[101:201],corr.bestClass$perf[1:101,101:201],col="slategray",alpha=1,
	xlim=c(-1.0,0.0),ylim=c(0,1),zlim=c(0,1),aspect=c(1,1,0.5),
	add=FALSE,box=TRUE,axes=FALSE,xlab='',ylab='',zlab='')
axes3d(c('x','y','z-+','z++'))
title3d(main="RECIST (Best Response) Model Performance, C metric",
	xlab="PR Cut",ylab="PD Cut")
spheres3d(x=-0.30,y=0.20,z=corr.bestClass$perf[tmpCurPR,tmpCurPD],color='red',rad=0.02)
spheres3d(x=corr.bestClass$max.PR,
	y=corr.bestClass$max.PD,
	z=corr.bestClass$max,color='green',rad=0.02)


rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))