# Analyze simulated patient data to see how the optimization compares to when we have real data.
#
# Before running this script, either:
# i) run the 'importData.R' script, followed by the 'generateSimulatedData.R' script.
# ii) Load pre-generated simulated patient data, using the 'importData_Simulation.R' script, or
#   manually loading it through the user interface.
#
# Once the simulated data is generated/loaded, this script runs the analysis on the simulated data.
#
# Written by E Gudmundsson, June 2017

###((Run 'importData_Simulation.R', followed by 'generateSimulatedData.R'###

library(Hmisc)

# vary the RECIST cutoffs and see what happens.

# set the file names
tmpBaseDir = file.path('C:','Users','Justin','Documents','_UC intern')
tmpProjName = 'AustraliaLinearResponse'
tmpDataDir = file.path(tmpBaseDir, tmpProjName)
tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

# Flags indicating how to optimize
flagMidPercentile = FALSE
midPercentile = 80
flagCloseToClinical = FALSE
flagCloseToClinicalNo0 = FALSE
source(file.path(tmpScriptDir, 'createThresholds.R'))

# Flag to indicate whether I want to vary the number of Sim patients 
# or just generate the data myself
flagVaryPatientNumber = FALSE
if (flagVaryPatientNumber == TRUE) {
    source(file.path(tmpScriptDir, 'varyNumberSimPatients.R')) 
    } else {
    source(file.path(tmpScriptDir, 'importData_Simulation.R'))

    tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
    tmpProjName = 'AustraliaLinearResponse'
    tmpDataDir = file.path(tmpBaseDir, tmpProjName)
    tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

        nSimulatedPatients <- 50
    source(file.path(tmpScriptDir, 'generateSimulatedData.R'))

        # calculate some of the continuous metrics, since those are independent of the loops
        corr.Continuous = list()
        corr.Continuous$bestChange = rcorr.cens(-simDataContinuous$RecistBest,
            Surv(simDataContinuous$SurvTimeDays, simDataContinuous$Status), outx = TRUE)
        corr.Continuous$lastChange = rcorr.cens(-simDataContinuous$RecistChangeN,
            Surv(simDataContinuous$SurvTimeDays, simDataContinuous$Status), outx = TRUE)
        corr.Continuous$firstChange = rcorr.cens(-simDataContinuous$RecistChange2,
            Surv(simDataContinuous$SurvTimeDays, simDataContinuous$Status), outx = TRUE)

    # flag to do train/test split
    flagTrainTest=FALSE
        if (flagTrainTest == TRUE) {

            tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
            tmpProjName = 'AustraliaLinearResponse'
            tmpDataDir = file.path(tmpBaseDir, tmpProjName)
            tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

            source(file.path(tmpScriptDir, 'trainTestCompute.R'))

        } else {

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
            



        }



}

save(list = c(
            ls(pattern = 'corr.+'),
            ls(pattern = 'logrank.+')),
            file = file.path(tmpDataDir, 'SimulatedDataSetPerformance.rdata'))
rm(list = ls(pattern = 'tmp.*'))
rm(list = ls(pattern = 'i[[:upper:]]'))



