mdlBest=coxph(Surv(BaseTimeDays,Status)~bestClass,data=newCut)
mdlStd=coxph(Surv(BaseTimeDays,Status)~stdBestClass,data=newCut)
cindex(list('best'=mdlBest,'std'=mdlStd),
	formula=Surv(BaseTimeDays,Status)~1,data=newCut,
	tiedPredictionsIn=TRUE)

with(newCut,rcorr.cens(bestClass,Surv(BaseTimeDays,Status),outx=TRUE))