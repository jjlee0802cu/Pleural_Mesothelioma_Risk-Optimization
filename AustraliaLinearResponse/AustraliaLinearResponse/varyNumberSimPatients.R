source(file.path(tmpScriptDir, 'importData_Simulation.R'))

tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
tmpProjName = 'AustraliaLinearResponse'
tmpDataDir = file.path(tmpBaseDir, tmpProjName)
tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

# Empty matrix to fill up with data later on
variedPatientNumbersData <- data.frame(matrix(, nrow = 0, ncol = 18))
# no. of Patients I will be testing
nVarySimulatedPatients =c(50, 60, 70,80,100,120,140,160,180,200,225,250,275,300,350,400,450,500)
for (i in 1:length(nVarySimulatedPatients)) {
    nSimulatedPatients <- nVarySimulatedPatients[i]

    tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
    tmpProjName = 'AustraliaLinearResponse'
    tmpDataDir = file.path(tmpBaseDir, tmpProjName)
    tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

    source(file.path(tmpScriptDir, 'generateSimulatedData.R'))

    # calculate some of the continuous metrics, since those are independent of the loops
    corr.Continuous = list()
    corr.Continuous$bestChange = rcorr.cens(-simDataContinuous$RecistBest,
        Surv(simDataContinuous$SurvTimeDays, simDataContinuous$Status), outx = TRUE)
    corr.Continuous$lastChange = rcorr.cens(-simDataContinuous$RecistChangeN,
        Surv(simDataContinuous$SurvTimeDays, simDataContinuous$Status), outx = TRUE)
    corr.Continuous$firstChange = rcorr.cens(-simDataContinuous$RecistChange2,
        Surv(simDataContinuous$SurvTimeDays, simDataContinuous$Status), outx = TRUE)

    tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
    tmpProjName = 'AustraliaLinearResponse'
    tmpDataDir = file.path(tmpBaseDir, tmpProjName)
    tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

    # do the variable cutpoint thing for the whole dataset
    source(file.path(tmpScriptDir, 'funcPerfMetrics.R'))
    tmpOutput = perf.metrics(corr.Cuts, simDataFullCols, simLinRECIST, classRECIST)
    corr.bestClass = tmpOutput$corr.bestClass
    corr.firstClass = tmpOutput$corr.firstClass
    corr.lastClass = tmpOutput$corr.lastClass
    logrank.bestClass = tmpOutput$logrank.bestClass
    logrank.firstClass = tmpOutput$logrank.firstClass
    logrank.lastClass = tmpOutput$logrank.lastClass


    #Populate the variedPatientNumbersData matrix with data for each number of patients
    variedPatientNumbersData = rbind(
    variedPatientNumbersData,
    c(
        nSimulatedPatients,
        corr.bestClass$max.stats$max.C,
        corr.bestClass$max.stats$max.C.npv,
        corr.bestClass$max.stats$max.C.ppv,
        corr.bestClass$max.stats$max.C.spec,
        corr.bestClass$max.stats$max.C.sens,
        corr.bestClass$max.stats$max.npv,
        corr.bestClass$max.stats$max.ppv,
        corr.bestClass$max.stats$max.spec,
        corr.bestClass$max.stats$max.sens,
        corr.bestClass$max.stats$max.C.npvNonResp,
        corr.bestClass$max.stats$max.C.ppvProgFree,
        corr.bestClass$max.stats$max.C.specNonResp,
        corr.bestClass$max.stats$max.C.sensProgFree,
        corr.bestClass$max.stats$max.npvNonResp,
        corr.bestClass$max.stats$max.ppvProgFree,
        corr.bestClass$max.stats$max.specNonResp,
        corr.bestClass$max.stats$max.sensProgFree,
corr.bestClass$nearest.comp$nearest.max.npv.dist,
corr.bestClass$nearest.comp$nearest.max.ppv.dist,
corr.bestClass$nearest.comp$nearest.max.spec.dist,
corr.bestClass$nearest.comp$nearest.max.sens.dist,

corr.bestClass$nearest.comp$nearest.max.npvNonResp.dist,
corr.bestClass$nearest.comp$nearest.max.ppvProgFree.dist,
corr.bestClass$nearest.comp$nearest.max.specNonResp.dist,
corr.bestClass$nearest.comp$nearest.max.sensProgFree.dist,

corr.bestClass$no.Patients.comp$PR,
corr.bestClass$no.Patients.comp$SD,
corr.bestClass$no.Patients.comp$PD


     )
    )

}
rm(i)

colnames(variedPatientNumbersData) <-
    c("nSimulatedPatients",
        "max.C",
        "NPV @ max.C",
        "PPV @ max.C",
        "spec @ max.C",
        "sens @ max.C",
        "max NPV",
        "max PPV",
        "max spec",
        "max sens",
        "npvNonResp @ max.C",
        "ppvProgFree @ max.C",
        "specNonResp @ max.C",
        "sensProgFree @ max.C",
        "max npvNonResp",
        "max ppvProgFree",
        "max specNonResp",
        "max sensProgFree",
        "nearest max NPV dist",
        "nearest max PPV dist",
        "nearest max Spec dist",
        "nearest max Sens dist",
        "nearest max NPV NonResp dist",
        "nearest max PPV ProgFree",
        "nearest max Spec NonResp",
        "nearest max Sens ProgFree",
        "number of PR patients",
        "number of SD patients",
        "number of PD patients"
        )
# flag to plot metrics c, npv,ppv,spec,sens all together.
flagPlotMetrics = TRUE
if (flagPlotMetrics == TRUE) {
    tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
    tmpProjName = 'AustraliaLinearResponse'
    tmpDataDir = file.path(tmpBaseDir, tmpProjName)
    tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

    source(file.path(tmpScriptDir, 'plotMetrics.R'))
} else { }