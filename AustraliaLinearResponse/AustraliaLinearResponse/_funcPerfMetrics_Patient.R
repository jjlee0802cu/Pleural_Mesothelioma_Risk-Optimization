perf.metrics = function(corr.Cuts, dataFullColsAnonymized, dataLinRECIST, classRECIST,tThresh) {

    library(Hmisc)
    library(doSNOW)
    cl <- makeCluster(2, type = "SOCK")
    registerDoSNOW(cl)


    tmpN = length(corr.Cuts)
    tmpMatrix = matrix(nrow = tmpN, ncol = tmpN, data = 0)

    # perf. is the matrix of P_k for varying RECIST cutoff points.  If the PD,PR
    # point is invalid (i.e., PD<PR), then just set P_k = 0.  The up and lo.95
    # matrices are the upper and lower 95% confidence intervals on the stat.
    # That multiplicative factor is 1.95996
    corr.bestClass = list(perf = tmpMatrix, up95 = tmpMatrix, lo95 = tmpMatrix, npv = tmpMatrix, ppv = tmpMatrix, spec = tmpMatrix, sens = tmpMatrix)
    # that first data frame is for doing the STANDARD recist algorithm.
    # let's also see how this performs using the last continuous change per patient,
    # and the first change per patient (from baseline to very first followup)
    corr.lastClass = list(perf = tmpMatrix, up95 = tmpMatrix, lo95 = tmpMatrix)
    corr.firstClass = list(perf = tmpMatrix, up95 = tmpMatrix, lo95 = tmpMatrix)
    # prepare the log-rank matrices too.
    logrank.bestClass = list(chisq = tmpMatrix, df = tmpMatrix, p = tmpMatrix + 1)
    logrank.firstClass = list(chisq = tmpMatrix, df = tmpMatrix, p = tmpMatrix + 1)
    logrank.lastClass = list(chisq = tmpMatrix, df = tmpMatrix, p = tmpMatrix + 1)

    # get the patients in this collection
    tmpListPats = names(dataLinRECIST)

    # make the starter data for the rest of these runs
    tmpStarterColsOfInterest = sort(
    c(
        match("URN", names(dataFullColsAnonymized)),
        match("Status", names(dataFullColsAnonymized)),
        match("SurvTimeDays", names(dataFullColsAnonymized)),
        match("EntryTimeDays", names(dataFullColsAnonymized)),
        match("BaseTimeDays", names(dataFullColsAnonymized))
        )
    )
    tmpStarterData = dataFullColsAnonymized[, tmpStarterColsOfInterest]

    # calculate the classification tasks for the array of cutpoints.
    # this is the part that can be done in parallel.
    # the net output of this section of code is the corr.bestClass,
    # corr.firstClass, and corr.lastClass variables.  That's it.
    tmpOpts = list(chunkSize = ceil(tmpN / getDoParWorkers()))



tmpParOutput=foreach(iRow=1:tmpN, .combine='rbind', .inorder=FALSE, .packages=c('survival','Hmisc'), .options.smp=tmpOpts) %:%
	foreach(iCol=1:tmpN, .combine='rbind', .inorder=FALSE, .packages=c('survival','Hmisc')) %dopar% {
		tmpCutPR = corr.Cuts[iRow]	# row tied to PR
        tmpCutPD = corr.Cuts[iCol] # column tied to PD


        if (tmpCutPD >= 0 && tmpCutPR <= 0) { # otherwise it doesn't make sense.

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
                levels = c('PD', 'SD', 'PR'), labels = c(1, 2, 3))
			tmpRow=match(as.character(iPat),tmpPatClass$URN)
			tmpPatClass$bestClass[tmpRow]=max(as.numeric(tmpPatRecord$classRECIST))
			tmpPatClass$lastClass[tmpRow]=tmpPatRecord$classRECIST[tmpNumUpdates]
			tmpPatClass$firstClass[tmpRow]=tmpPatRecord$classRECIST[1]
        }


        ### Built in c-stat ###
		# calculate the model performance at this cutpoint
		tmpCorrObj1=rcorr.cens(tmpPatClass$bestClass,
			Surv(tmpPatClass$SurvTimeDays,tmpPatClass$Status),outx=TRUE)
		
		tmpCorrObj2=rcorr.cens(tmpPatClass$lastClass,
			Surv(tmpPatClass$SurvTimeDays,tmpPatClass$Status),outx=TRUE)
		
		tmpCorrObj3=rcorr.cens(tmpPatClass$firstClass,
			Surv(tmpPatClass$SurvTimeDays,tmpPatClass$Status),outx=TRUE)
   
            
		# calculate the log-rank tests at this cutpoint
		if (length(unique(tmpPatClass$bestClass))>1) {
			tmpLogRank1=survdiff(Surv(SurvTimeDays,Status)~bestClass,data=tmpPatClass)
			tmpLRdf1=length(tmpLogRank1$obs)-1
			tmpLRcs1=tmpLogRank1$chisq
			tmpLRp1=pchisq(tmpLRcs1,tmpLRdf1,lower.tail=FALSE)
		}
		else {
			tmpLRdf1=0; tmpLRcs1=0; tmpLRp1=1;
		}
		
		if (length(unique(tmpPatClass$lastClass))>1) {
			tmpLogRank2=survdiff(Surv(SurvTimeDays,Status)~lastClass,data=tmpPatClass)
			tmpLRdf2=length(tmpLogRank2$obs)-1
			tmpLRcs2=tmpLogRank2$chisq
			tmpLRp2=pchisq(tmpLRcs2,tmpLRdf2,lower.tail=FALSE)
		}
		else {
			tmpLRdf2=0; tmpLRcs2=0; tmpLRp2=1;
		}
		
		if (length(unique(tmpPatClass$firstClass))>1) {
			tmpLogRank3=survdiff(Surv(SurvTimeDays,Status)~firstClass,data=tmpPatClass)
			tmpLRdf3=length(tmpLogRank3$obs)-1
			tmpLRcs3=tmpLogRank3$chisq
			tmpLRp3=pchisq(tmpLRcs3,tmpLRdf3,lower.tail=FALSE)
		}
		else {
			tmpLRdf3=0; tmpLRcs3=0; tmpLRp3=1;
		}
### First splitting patients only as PR vs PD
pPatientsPR <- tmpPatClass[tmpPatClass$bestClass == 3, c("SurvTimeDays", "bestClass")]
pPatientsPD <- tmpPatClass[tmpPatClass$bestClass == 1, c("SurvTimeDays", "bestClass")]


### Calculating NPV, PPV, Sensitivity, Specificity
#PPV: 
ppv=sum(pPatientsPR[, 1] > tThresh) / nrow(pPatientsPR)
#NPV:
npv=sum(pPatientsPD[, 1] < tThresh) / nrow(pPatientsPD)
#Sensitivity
sens=sum(pPatientsPR[, 1] > tThresh) / (sum(pPatientsPR[, 1] > tThresh) + sum(pPatientsPD[, 1] > tThresh))
#Specificity
spec=sum(pPatientsPD[, 1] < tThresh) / (sum(pPatientsPD[, 1] < tThresh) + sum(pPatientsPR[, 1] < tThresh))

### Next splitting patients as (PR+SD) vs PD or PR vs (SD+PD)
pPatientsSD <- tmpPatClass[tmpPatClass$bestClass == 2, c("SurvTimeDays", "bestClass")]

# PPVprogfree: (PR+SD) vs PD
ppvProgFree=(sum(pPatientsPR[, 1] > tThresh)+sum(pPatientsSD[, 1] > tThresh))/(nrow(pPatientsPR)+nrow(pPatientsSD))
# NPVnonresp: PR vs (SD+PD)
npvNonResp=(sum(pPatientsSD[, 1] < tThresh)+sum(pPatientsPD[, 1] < tThresh))/(nrow(pPatientsSD)+nrow(pPatientsPD))
# Sensprogfree: (PR+SD) vs PD
sensProgFree=(sum(pPatientsPR[, 1] > tThresh)+sum(pPatientsSD[, 1] > tThresh))/(sum(pPatientsPR[, 1] > tThresh)+sum(pPatientsSD[,1]>tThresh) + sum(pPatientsPD[, 1] > tThresh))                
# Specnonresp: PR vs (SD+PD)
specNonResp=(sum(pPatientsSD[, 1] < tThresh)+sum(pPatientsPD[, 1] < tThresh))/(sum(pPatientsPR[, 1] < tThresh)+sum(pPatientsSD[,1]<tThresh) + sum(pPatientsPD[, 1] < tThresh))

### Weighted c-statistic
flagWeightedCStat=TRUE
            pPatientsPRStatus <- tmpPatClass[tmpPatClass$bestClass == 3, c("SurvTimeDays", "bestClass","Status")]
            pPatientsPDStatus <- tmpPatClass[tmpPatClass$bestClass == 1, c("SurvTimeDays", "bestClass", "Status")]
            pPatientsSDStatus <- tmpPatClass[tmpPatClass$bestClass == 2, c("SurvTimeDays", "bestClass", "Status")]

            sum = 0
            comparisons = 0

            if (nrow(pPatientsPRStatus) > 0) {
                if (nrow(pPatientsSDStatus) > 0) {


                    for (i in c(1:nrow(pPatientsPRStatus))) {
                        for (j in c(1:nrow(pPatientsSDStatus))) {

                            if (pPatientsPRStatus[[i, 3]] == 0 && pPatientsSDStatus[[j,3]]==1) {
                                if (pPatientsPRStatus[[i, 1]] >= pPatientsSDStatus[[j, 1]]) {
                                    comparisons<-comparisons+1
                                    sum <- sum+1
                                } else {      
                                }
                                } else {
                                    if (pPatientsPRStatus[[i, 3]] == 1 && pPatientsSDStatus[[j, 3]] == 0) {
                                        if (pPatientsPRStatus[[i, 1]] <= pPatientsSDStatus[[j, 1]]) {
                                            comparisons<-comparisons+1
                                        } else {            
                                        }
                                    } else {
                                        if (pPatientsPRStatus[[i, 3]] == 0 && pPatientsSDStatus[[j, 3]] == 0) {         
                                        } else {
                                            if (pPatientsPRStatus[[i, 3]] == 1 && pPatientsSDStatus[[j, 3]] == 1) {
                                                if (pPatientsPRStatus[[i, 1]] > pPatientsSDStatus[[j, 1]]) {
                                                    comparisons <- comparisons + 1
                                                    sum <- sum + 1
                                                } else if (pPatientsPRStatus[[i, 1]] == pPatientsSDStatus[[j, 1]]) { }
                                                else { comparisons <- comparisons + 1 }
                                            }
                                        }
                                            }
                                }
                        }
                    }

                } else { }
                } else { }
            rm(i)
            rm(j)

            if (nrow(pPatientsSDStatus) > 0) {
                if (nrow(pPatientsPDStatus) > 0) {


                    for (i in c(1:nrow(pPatientsSDStatus))) {
                        for (j in c(1:nrow(pPatientsPDStatus))) {

                            if (pPatientsSDStatus[[i, 3]] == 0 && pPatientsPDStatus[[j, 3]] == 1) {
                                if (pPatientsSDStatus[[i, 1]] >= pPatientsPDStatus[[j, 1]]) {
                                    comparisons <- comparisons + 1
                                    sum <- sum + 1
                                } else {
                                }
                            } else {
                                if (pPatientsSDStatus[[i, 3]] == 1 && pPatientsPDStatus[[j, 3]] == 0) {
                                    if (pPatientsSDStatus[[i, 1]] <= pPatientsPDStatus[[j, 1]]) {
                                        comparisons <- comparisons + 1
                                    } else {
                                    }
                                } else {
                                    if (pPatientsSDStatus[[i, 3]] == 0 && pPatientsPDStatus[[j, 3]] == 0) {
                                    } else {
                                        if (pPatientsSDStatus[[i, 3]] == 1 && pPatientsPDStatus[[j, 3]] == 1) {
                                            if (pPatientsSDStatus[[i, 1]] > pPatientsPDStatus[[j, 1]]) {
                                                comparisons <- comparisons + 1
                                                sum <- sum + 1
                                            } else if (pPatientsSDStatus[[i, 1]] == pPatientsPDStatus[[j, 1]]) { } else { comparisons <- comparisons + 1 }
                                            }
                                    }
                                }
                            }
                        }
                    }

                } else { }
                } else { }
            rm(i)
            rm(j)



            if (nrow(pPatientsPRStatus) > 0) {
                if (nrow(pPatientsPDStatus) > 0) {


                    for (i in c(1:nrow(pPatientsPRStatus))) {
                        for (j in c(1:nrow(pPatientsPDStatus))) {

                            if (pPatientsPRStatus[[i, 3]] == 0 && pPatientsPDStatus[[j, 3]] == 1) {
                                if (pPatientsPRStatus[[i, 1]] >= pPatientsPDStatus[[j, 1]]) {
                                    comparisons <- comparisons + 1
                                    sum <- sum + 1
                                } else {
                                }
                            } else {
                                if (pPatientsPRStatus[[i, 3]] == 1 && pPatientsPDStatus[[j, 3]] == 0) {
                                    if (pPatientsPRStatus[[i, 1]] <= pPatientsPDStatus[[j, 1]]) {
                                        comparisons <- comparisons + 1
                                        if (flagWeightedCStat == TRUE) {
                                            comparisons <- comparisons + 1
                                            #Here is the weight: if wrong about PR vs PD, add another one to comparisons
                                        }
                                    } else {
                                    }
                                } else {
                                    if (pPatientsPRStatus[[i, 3]] == 0 && pPatientsPDStatus[[j, 3]] == 0) {
                                    } else {
                                        if (pPatientsPRStatus[[i, 3]] == 1 && pPatientsPDStatus[[j, 3]] == 1) {
                                            if (pPatientsPRStatus[[i, 1]] > pPatientsPDStatus[[j, 1]]) {
                                                comparisons <- comparisons + 1
                                                sum <- sum + 1
                                            } else if (pPatientsPRStatus[[i, 1]] == pPatientsPDStatus[[j, 1]]) {
                                            } else {
                                                comparisons <- comparisons + 1
                                                if (flagWeightedCStat == TRUE) {
                                                    comparisons <- comparisons + 1
                                                    #Here is the weight: if wrong about PR vs PD, add another one to comparisons
                                                }
                                                }
                                            }
                                    }
                                }
                            }
                        }
                    }

                } else { }
                } else { }
            rm(i)
            rm(j)
c.weighted=sum/comparisons

            #Find the number of ppl in each group
            no.pPatientsPR = nrow(pPatientsPR)
            no.pPatientsSD = nrow(pPatientsSD)
            no.pPatientsPD = nrow(pPatientsPD)

#Putting all data into data.frame
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
            lr.FCcs = tmpLRcs3, lr.FCdf = tmpLRdf3, lr.FCp = tmpLRp3,
            ppv = ppv,
            npv = npv,
            sens = sens,
            spec = spec,
            ppvProgFree = ppvProgFree,
            npvNonResp = npvNonResp,
            sensProgFree = sensProgFree,
            specNonResp = specNonResp,
            c.weighted = c.weighted,
            no.pPatientsPR = no.pPatientsPR,
            no.pPatientsSD = no.pPatientsSD,
            no.pPatientsPD = no.pPatientsPD
            )
		
		}
		else {
			data.frame(row=iRow,col=iCol,
			bcP=0,bcU=0,bcL=0,lcP=0,lcU=0,lcL=0,fcP=0,fcU=0,fcL=0,
			lr.BCcs=0,lr.BCdf=0,lr.BCp=1,
			lr.LCcs=0,lr.LCdf=0,lr.LCp=1,
            lr.FCcs = 0, lr.FCdf = 0, lr.FCp = 1,
            ppv = 0, npv = 0, sens = 0, spec = 0,
            ppvProgFree = 0,
            npvNonResp = 0,
            sensProgFree = 0,
            specNonResp = 0,
            c.weighted = 0,
            no.pPatientsPR = 0,
            no.pPatientsSD = 0,
            no.pPatientsPD = 0
            )
        }


}
stopCluster(cl)

#Populating corr.bestClass with new the variables
tmpParOutput = tmpParOutput[order(tmpParOutput$row, tmpParOutput$col),]
corr.bestClass$no.pPatientsPR = matrix(tmpParOutput$no.pPatientsPR, nrow = tmpN, ncol = tmpN, byrow = TRUE)
corr.bestClass$no.pPatientsSD = matrix(tmpParOutput$no.pPatientsSD, nrow = tmpN, ncol = tmpN, byrow = TRUE)
corr.bestClass$no.pPatientsPD = matrix(tmpParOutput$no.pPatientsPD, nrow = tmpN, ncol = tmpN, byrow = TRUE)

corr.bestClass$c.weighted = matrix(tmpParOutput$c.weighted, nrow = tmpN, ncol = tmpN, byrow = TRUE)

corr.bestClass$npv = matrix(tmpParOutput$npv, nrow = tmpN, ncol = tmpN, byrow = TRUE)
corr.bestClass$ppv = matrix(tmpParOutput$ppv, nrow = tmpN, ncol = tmpN, byrow = TRUE)
corr.bestClass$spec = matrix(tmpParOutput$spec, nrow = tmpN, ncol = tmpN, byrow = TRUE)
corr.bestClass$sens = matrix(tmpParOutput$sens, nrow = tmpN, ncol = tmpN, byrow = TRUE)

corr.bestClass$npvNonResp=matrix(tmpParOutput$npvNonResp, nrow = tmpN, ncol = tmpN, byrow = TRUE)
corr.bestClass$ppvProgFree = matrix(tmpParOutput$ppvProgFree, nrow = tmpN, ncol = tmpN, byrow = TRUE)
corr.bestClass$specNonResp = matrix(tmpParOutput$specNonResp, nrow = tmpN, ncol = tmpN, byrow = TRUE)
corr.bestClass$sensProgFree = matrix(tmpParOutput$sensProgFree, nrow = tmpN, ncol = tmpN, byrow = TRUE)

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



##### Nearest Metric Computation ####

    ### Since there can be many thresholds that can produce the same max c-statistic, 
    ### this finds the c-statistic closest to the clinical standard (-0.3,0.2)
    a = list()
    for (i in c(1:nrow(arrayInd(which(corr.bestClass$perf == max(corr.bestClass$perf, na.rm = TRUE)), dim(corr.bestClass$perf))))) {
        a <- c(a,
        ((((corr.Cuts[arrayInd(which(corr.bestClass$perf == max(corr.bestClass$perf, na.rm = TRUE)), dim(corr.bestClass$perf))[i, 1]] - (-0.3)) ^ 2) +
        ((corr.Cuts[arrayInd(which(corr.bestClass$perf == max(corr.bestClass$perf, na.rm = TRUE)), dim(corr.bestClass$perf))[i, 2]] - (0.2)) ^ 2)) ^ 0.5)
    )
    }
    ## Cutoffs that are closest to clinical standard
    corr.bestClass$max.stats$max.C.PR = corr.Cuts[arrayInd(which(corr.bestClass$perf == max(corr.bestClass$perf, na.rm = TRUE)), dim(corr.bestClass$perf))[which.min(a), 1]]
    corr.bestClass$max.stats$max.C.PD = corr.Cuts[arrayInd(which(corr.bestClass$perf == max(corr.bestClass$perf, na.rm = TRUE)), dim(corr.bestClass$perf))[which.min(a), 2]]
    ## Their respective indicies
    corr.bestClass$max.stats$max.indPR = arrayInd(which(corr.bestClass$perf == max(corr.bestClass$perf, na.rm = TRUE)), dim(corr.bestClass$perf))[which.min(a), 1]
    corr.bestClass$max.stats$max.indPD = arrayInd(which(corr.bestClass$perf == max(corr.bestClass$perf, na.rm = TRUE)), dim(corr.bestClass$perf))[which.min(a), 2]
    # the value of the c-statistic    
    corr.bestClass$max.stats$max.C = corr.bestClass$perf[
    arrayInd(which(corr.bestClass$perf == max(corr.bestClass$perf, na.rm = TRUE)), dim(corr.bestClass$perf))[which.min(a), 1],
    arrayInd(which(corr.bestClass$perf == max(corr.bestClass$perf, na.rm = TRUE)), dim(corr.bestClass$perf))[which.min(a), 2]
    ]
    rm(i)

corr.bestClass$nearest.comp=list()
### Finding max NPV, PPV, Spec, Sens that are closest in cutoff to the optimal c's cutoffs
    a = list()
    for (i in c(1:nrow(arrayInd(which(corr.bestClass$npv == max(corr.bestClass$npv, na.rm = TRUE)), dim(corr.bestClass$npv))))) {
        a <- c(a,
        ((((corr.Cuts[arrayInd(which(corr.bestClass$npv == max(corr.bestClass$npv, na.rm = TRUE)), dim(corr.bestClass$npv))[i, 1]] - (corr.bestClass$max.stats$max.C.PR)) ^ 2) +
        ((corr.Cuts[arrayInd(which(corr.bestClass$npv == max(corr.bestClass$npv, na.rm = TRUE)), dim(corr.bestClass$npv))[i, 2]] - (corr.bestClass$max.stats$max.C.PD)) ^ 2)) ^ 0.5)
    )
    }
    # Max npv
    corr.bestClass$nearest.comp$max.npv=corr.bestClass$npv[
        arrayInd(which(corr.bestClass$npv == max(corr.bestClass$npv, na.rm = TRUE)), dim(corr.bestClass$npv))[which.min(a), 1],
        arrayInd(which(corr.bestClass$npv == max(corr.bestClass$npv, na.rm = TRUE)), dim(corr.bestClass$npv))[which.min(a), 2]
    ]
    # Cutoffs that are closest to optimal c's cutoffs
    corr.bestClass$nearest.comp$nearest.max.npv=corr.Cuts[arrayInd(which(corr.bestClass$npv == max(corr.bestClass$npv, na.rm = TRUE)), dim(corr.bestClass$npv))[which.min(a),]]
    # Distance to optimal c cutoff
    corr.bestClass$nearest.comp$nearest.max.npv.dist=a[[which.min(a)]]
    ## Their respective indicies
    #arrayInd(which(corr.bestClass$npv == max(corr.bestClass$npv, na.rm = TRUE)), dim(corr.bestClass$npv))[which.min(a),]
    rm(i)


    a = list()
    for (i in c(1:nrow(arrayInd(which(corr.bestClass$ppv == max(corr.bestClass$ppv, na.rm = TRUE)), dim(corr.bestClass$ppv))))) {
        a <- c(a,
        ((((corr.Cuts[arrayInd(which(corr.bestClass$ppv == max(corr.bestClass$ppv, na.rm = TRUE)), dim(corr.bestClass$ppv))[i, 1]] - (corr.bestClass$max.stats$max.C.PR)) ^ 2) +
        ((corr.Cuts[arrayInd(which(corr.bestClass$ppv == max(corr.bestClass$ppv, na.rm = TRUE)), dim(corr.bestClass$ppv))[i, 2]] - (corr.bestClass$max.stats$max.C.PD)) ^ 2)) ^ 0.5)
    )
    }
    # Max ppv
    corr.bestClass$nearest.comp$max.ppv = corr.bestClass$ppv[
        arrayInd(which(corr.bestClass$ppv == max(corr.bestClass$ppv, na.rm = TRUE)), dim(corr.bestClass$ppv))[which.min(a), 1],
        arrayInd(which(corr.bestClass$ppv == max(corr.bestClass$ppv, na.rm = TRUE)), dim(corr.bestClass$ppv))[which.min(a), 2]
    ]
    # Cutoffs that are closest to optimal c's cutoffs
    corr.bestClass$nearest.comp$nearest.max.ppv = corr.Cuts[arrayInd(which(corr.bestClass$ppv == max(corr.bestClass$ppv, na.rm = TRUE)), dim(corr.bestClass$ppv))[which.min(a),]]
    # Distance to optimal c cutoff
    corr.bestClass$nearest.comp$nearest.max.ppv.dist = a[[which.min(a)]]
    ## Their respective indicies
    #arrayInd(which(corr.bestClass$ppv == max(corr.bestClass$ppv, na.rm = TRUE)), dim(corr.bestClass$ppv))[which.min(a),]
    rm(i)


    a = list()
    for (i in c(1:nrow(arrayInd(which(corr.bestClass$spec == max(corr.bestClass$spec, na.rm = TRUE)), dim(corr.bestClass$spec))))) {
        a <- c(a,
        ((((corr.Cuts[arrayInd(which(corr.bestClass$spec == max(corr.bestClass$spec, na.rm = TRUE)), dim(corr.bestClass$spec))[i, 1]] - (corr.bestClass$max.stats$max.C.PR)) ^ 2) +
        ((corr.Cuts[arrayInd(which(corr.bestClass$spec == max(corr.bestClass$spec, na.rm = TRUE)), dim(corr.bestClass$spec))[i, 2]] - (corr.bestClass$max.stats$max.C.PD)) ^ 2)) ^ 0.5)
    )
    }
    # Max spec
    corr.bestClass$nearest.comp$max.spec = corr.bestClass$spec[
        arrayInd(which(corr.bestClass$spec == max(corr.bestClass$spec, na.rm = TRUE)), dim(corr.bestClass$spec))[which.min(a), 1],
        arrayInd(which(corr.bestClass$spec == max(corr.bestClass$spec, na.rm = TRUE)), dim(corr.bestClass$spec))[which.min(a), 2]
    ]
    # Cutoffs that are closest to optimal c's cutoffs
    corr.bestClass$nearest.comp$nearest.max.spec = corr.Cuts[arrayInd(which(corr.bestClass$spec == max(corr.bestClass$spec, na.rm = TRUE)), dim(corr.bestClass$spec))[which.min(a),]]
    # Distance to optimal c cutoff
    corr.bestClass$nearest.comp$nearest.max.spec.dist = a[[which.min(a)]]
    ## Their respective indicies
    #arrayInd(which(corr.bestClass$spec == max(corr.bestClass$spec, na.rm = TRUE)), dim(corr.bestClass$spec))[which.min(a),]
    rm(i)


    a = list()
    for (i in c(1:nrow(arrayInd(which(corr.bestClass$sens == max(corr.bestClass$sens, na.rm = TRUE)), dim(corr.bestClass$sens))))) {
        a <- c(a,
        ((((corr.Cuts[arrayInd(which(corr.bestClass$sens == max(corr.bestClass$sens, na.rm = TRUE)), dim(corr.bestClass$sens))[i, 1]] - (corr.bestClass$max.stats$max.C.PR)) ^ 2) +
        ((corr.Cuts[arrayInd(which(corr.bestClass$sens == max(corr.bestClass$sens, na.rm = TRUE)), dim(corr.bestClass$sens))[i, 2]] - (corr.bestClass$max.stats$max.C.PD)) ^ 2)) ^ 0.5)
    )
    }
    # Max sens
    corr.bestClass$nearest.comp$max.sens = corr.bestClass$sens[
        arrayInd(which(corr.bestClass$sens == max(corr.bestClass$sens, na.rm = TRUE)), dim(corr.bestClass$sens))[which.min(a), 1],
        arrayInd(which(corr.bestClass$sens == max(corr.bestClass$sens, na.rm = TRUE)), dim(corr.bestClass$sens))[which.min(a), 2]
    ]
    # Cutoffs that are closest to optimal c's cutoffs
    corr.bestClass$nearest.comp$nearest.max.sens = corr.Cuts[arrayInd(which(corr.bestClass$sens == max(corr.bestClass$sens, na.rm = TRUE)), dim(corr.bestClass$sens))[which.min(a),]]
    # Distance to optimal c cutoff
    corr.bestClass$nearest.comp$nearest.max.sens.dist = a[[which.min(a)]]
    ## Their respective indicies
    #arrayInd(which(corr.bestClass$sens == max(corr.bestClass$sens, na.rm = TRUE)), dim(corr.bestClass$sens))[which.min(a),]
    rm(i)
    rm(a)

    a = list()
    for (i in c(1:nrow(arrayInd(which(corr.bestClass$npvNonResp == max(corr.bestClass$npvNonResp, na.rm = TRUE)), dim(corr.bestClass$npvNonResp))))) {
        a <- c(a,
        ((((corr.Cuts[arrayInd(which(corr.bestClass$npvNonResp == max(corr.bestClass$npvNonResp, na.rm = TRUE)), dim(corr.bestClass$npvNonResp))[i, 1]] - (corr.bestClass$max.stats$max.C.PR)) ^ 2) +
        ((corr.Cuts[arrayInd(which(corr.bestClass$npvNonResp == max(corr.bestClass$npvNonResp, na.rm = TRUE)), dim(corr.bestClass$npvNonResp))[i, 2]] - (corr.bestClass$max.stats$max.C.PD)) ^ 2)) ^ 0.5)
    )
    }
    # Max npvNonResp
    corr.bestClass$nearest.comp$max.npvNonResp = corr.bestClass$npvNonResp[
        arrayInd(which(corr.bestClass$npvNonResp == max(corr.bestClass$npvNonResp, na.rm = TRUE)), dim(corr.bestClass$npvNonResp))[which.min(a), 1],
        arrayInd(which(corr.bestClass$npvNonResp == max(corr.bestClass$npvNonResp, na.rm = TRUE)), dim(corr.bestClass$npvNonResp))[which.min(a), 2]
    ]
    # Cutoffs that are closest to optimal c's cutoffs
    corr.bestClass$nearest.comp$nearest.max.npvNonResp = corr.Cuts[arrayInd(which(corr.bestClass$npvNonResp == max(corr.bestClass$npvNonResp, na.rm = TRUE)), dim(corr.bestClass$npvNonResp))[which.min(a),]]
    # Distance to optimal c cutoff
    corr.bestClass$nearest.comp$nearest.max.npvNonResp.dist = a[[which.min(a)]]
    ## Their respective indicies
    #arrayInd(which(corr.bestClass$npvNonResp == max(corr.bestClass$npvNonResp, na.rm = TRUE)), dim(corr.bestClass$npvNonResp))[which.min(a),]
    rm(i)
    rm(a)

    a = list()
    for (i in c(1:nrow(arrayInd(which(corr.bestClass$ppvProgFree == max(corr.bestClass$ppvProgFree, na.rm = TRUE)), dim(corr.bestClass$ppvProgFree))))) {
        a <- c(a,
        ((((corr.Cuts[arrayInd(which(corr.bestClass$ppvProgFree == max(corr.bestClass$ppvProgFree, na.rm = TRUE)), dim(corr.bestClass$ppvProgFree))[i, 1]] - (corr.bestClass$max.stats$max.C.PR)) ^ 2) +
        ((corr.Cuts[arrayInd(which(corr.bestClass$ppvProgFree == max(corr.bestClass$ppvProgFree, na.rm = TRUE)), dim(corr.bestClass$ppvProgFree))[i, 2]] - (corr.bestClass$max.stats$max.C.PD)) ^ 2)) ^ 0.5)
    )
    }
    # Max ppvProgFree
    corr.bestClass$nearest.comp$max.ppvProgFree = corr.bestClass$ppvProgFree[
        arrayInd(which(corr.bestClass$ppvProgFree == max(corr.bestClass$ppvProgFree, na.rm = TRUE)), dim(corr.bestClass$ppvProgFree))[which.min(a), 1],
        arrayInd(which(corr.bestClass$ppvProgFree == max(corr.bestClass$ppvProgFree, na.rm = TRUE)), dim(corr.bestClass$ppvProgFree))[which.min(a), 2]
    ]
    # Cutoffs that are closest to optimal c's cutoffs
    corr.bestClass$nearest.comp$nearest.max.ppvProgFree = corr.Cuts[arrayInd(which(corr.bestClass$ppvProgFree == max(corr.bestClass$ppvProgFree, na.rm = TRUE)), dim(corr.bestClass$ppvProgFree))[which.min(a),]]
    # Distance to optimal c cutoff
    corr.bestClass$nearest.comp$nearest.max.ppvProgFree.dist = a[[which.min(a)]]
    ## Their respective indicies
    #arrayInd(which(corr.bestClass$ppvProgFree == max(corr.bestClass$ppvProgFree, na.rm = TRUE)), dim(corr.bestClass$ppvProgFree))[which.min(a),]
    rm(i)
    rm(a)

    a = list()
    for (i in c(1:nrow(arrayInd(which(corr.bestClass$specNonResp == max(corr.bestClass$specNonResp, na.rm = TRUE)), dim(corr.bestClass$specNonResp))))) {
        a <- c(a,
        ((((corr.Cuts[arrayInd(which(corr.bestClass$specNonResp == max(corr.bestClass$specNonResp, na.rm = TRUE)), dim(corr.bestClass$specNonResp))[i, 1]] - (corr.bestClass$max.stats$max.C.PR)) ^ 2) +
        ((corr.Cuts[arrayInd(which(corr.bestClass$specNonResp == max(corr.bestClass$specNonResp, na.rm = TRUE)), dim(corr.bestClass$specNonResp))[i, 2]] - (corr.bestClass$max.stats$max.C.PD)) ^ 2)) ^ 0.5)
    )
    }
    # Max specNonResp
    corr.bestClass$nearest.comp$max.specNonResp = corr.bestClass$specNonResp[
        arrayInd(which(corr.bestClass$specNonResp == max(corr.bestClass$specNonResp, na.rm = TRUE)), dim(corr.bestClass$specNonResp))[which.min(a), 1],
        arrayInd(which(corr.bestClass$specNonResp == max(corr.bestClass$specNonResp, na.rm = TRUE)), dim(corr.bestClass$specNonResp))[which.min(a), 2]
    ]
    # Cutoffs that are closest to optimal c's cutoffs
    corr.bestClass$nearest.comp$nearest.max.specNonResp = corr.Cuts[arrayInd(which(corr.bestClass$specNonResp == max(corr.bestClass$specNonResp, na.rm = TRUE)), dim(corr.bestClass$specNonResp))[which.min(a),]]
    # Distance to optimal c cutoff
    corr.bestClass$nearest.comp$nearest.max.specNonResp.dist = a[[which.min(a)]]
    ## Their respective indicies
    #arrayInd(which(corr.bestClass$specNonResp == max(corr.bestClass$specNonResp, na.rm = TRUE)), dim(corr.bestClass$specNonResp))[which.min(a),]
    rm(i)
    rm(a)

    a = list()
    for (i in c(1:nrow(arrayInd(which(corr.bestClass$sensProgFree == max(corr.bestClass$sensProgFree, na.rm = TRUE)), dim(corr.bestClass$sensProgFree))))) {
        a <- c(a,
        ((((corr.Cuts[arrayInd(which(corr.bestClass$sensProgFree == max(corr.bestClass$sensProgFree, na.rm = TRUE)), dim(corr.bestClass$sensProgFree))[i, 1]] - (corr.bestClass$max.stats$max.C.PR)) ^ 2) +
        ((corr.Cuts[arrayInd(which(corr.bestClass$sensProgFree == max(corr.bestClass$sensProgFree, na.rm = TRUE)), dim(corr.bestClass$sensProgFree))[i, 2]] - (corr.bestClass$max.stats$max.C.PD)) ^ 2)) ^ 0.5)
    )
    }
    # Max sensProgFree
    corr.bestClass$nearest.comp$max.sensProgFree = corr.bestClass$sensProgFree[
        arrayInd(which(corr.bestClass$sensProgFree == max(corr.bestClass$sensProgFree, na.rm = TRUE)), dim(corr.bestClass$sensProgFree))[which.min(a), 1],
        arrayInd(which(corr.bestClass$sensProgFree == max(corr.bestClass$sensProgFree, na.rm = TRUE)), dim(corr.bestClass$sensProgFree))[which.min(a), 2]
    ]
    # Cutoffs that are closest to optimal c's cutoffs
    corr.bestClass$nearest.comp$nearest.max.sensProgFree = corr.Cuts[arrayInd(which(corr.bestClass$sensProgFree == max(corr.bestClass$sensProgFree, na.rm = TRUE)), dim(corr.bestClass$sensProgFree))[which.min(a),]]
    # Distance to optimal c cutoff
    corr.bestClass$nearest.comp$nearest.max.sensProgFree.dist = a[[which.min(a)]]
    ## Their respective indicies
    #arrayInd(which(corr.bestClass$sensProgFree == max(corr.bestClass$sensProgFree, na.rm = TRUE)), dim(corr.bestClass$sensProgFree))[which.min(a),]
    rm(i)
    rm(a)


# max npv, ppv, spec,sens
corr.bestClass$max.stats$max.npv = max(corr.bestClass$npv, na.rm = TRUE)
corr.bestClass$max.stats$max.ppv = max(corr.bestClass$ppv, na.rm = TRUE)
corr.bestClass$max.stats$max.spec = max(corr.bestClass$spec, na.rm = TRUE)
corr.bestClass$max.stats$max.sens = max(corr.bestClass$sens, na.rm = TRUE)
# max npvNonResp, ppvProgFree,sensProgFree,specNonResp
corr.bestClass$max.stats$max.npvNonResp = max(corr.bestClass$npvNonResp, na.rm = TRUE)
corr.bestClass$max.stats$max.ppvProgFree = max(corr.bestClass$ppvProgFree, na.rm = TRUE)
corr.bestClass$max.stats$max.specNonResp = max(corr.bestClass$specNonResp, na.rm = TRUE)
corr.bestClass$max.stats$max.sensProgFree = max(corr.bestClass$sensProgFree, na.rm = TRUE)
# max C
#corr.bestClass$max.stats$max.C = max(corr.bestClass$perf, na.rm = TRUE)
# these just show where in the list of corr.Cuts, the max c was found (which row/column to look at when looking at perf matrix)
#corr.bestClass$max.stats$max.indPR = arrayInd(which.max(corr.bestClass$perf), dim(corr.bestClass$perf))[1]
#corr.bestClass$max.stats$max.indPD = arrayInd(which.max(corr.bestClass$perf), dim(corr.bestClass$perf))[2]
# max PR and PD when max C
#corr.bestClass$max.stats$max.C.PR = corr.Cuts[corr.bestClass$max.stats$max.indPR]
#corr.bestClass$max.stats$max.C.PD = corr.Cuts[corr.bestClass$max.stats$max.indPD]

# max npv,ppv,spec,sens when max C
corr.bestClass$max.stats$max.C.npv = corr.bestClass$npv[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
corr.bestClass$max.stats$max.C.ppv = corr.bestClass$ppv[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
corr.bestClass$max.stats$max.C.spec = corr.bestClass$spec[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
corr.bestClass$max.stats$max.C.sens = corr.bestClass$sens[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
## max npvNonResp, ppvProgFree,sensProgFree,specNonResp when max C
corr.bestClass$max.stats$max.C.npvNonResp = corr.bestClass$npvNonResp[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
corr.bestClass$max.stats$max.C.ppvProgFree = corr.bestClass$ppvProgFree[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
corr.bestClass$max.stats$max.C.specNonResp = corr.bestClass$specNonResp[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
corr.bestClass$max.stats$max.C.sensProgFree = corr.bestClass$sensProgFree[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
# standev of all the c-stats
corr.bestClass$perf.stdev = with(corr.bestClass, (up95 - perf) / 1.95996)

#Number of patients in each PR, SD, PD when using the cutoffs of the maximum c-statistic
corr.bestClass$no.Patients.comp=list()
    corr.bestClass$no.Patients.comp$PR = corr.bestClass$no.pPatientsPR[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
    corr.bestClass$no.Patients.comp$SD = corr.bestClass$no.pPatientsSD[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]
    corr.bestClass$no.Patients.comp$PD = corr.bestClass$no.pPatientsPD[corr.bestClass$max.stats$max.indPR, corr.bestClass$max.stats$max.indPD]


    
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
    logrank.firstClass = logrank.firstClass
    ))

}

