# Generate simulated patient data to see how the optimization of classification thresholds compares 
# to when we have real data.
#
# Before running this script, run the 'importData.R' script. This script then modifies the dataLinRECIST
# data array. Alternatively, run the 'importData_Simulation.R' script to import pre-generated simulated data
# (in that case, this script doesn't need to be run before analyzing the simulated data).
#
#
# Written by E Gudmundsson, June 2017

# vary the RECIST cutoffs and see what happens.
tmpBaseDir = file.path('F:', 'MesoResponse')
tmpProjName = 'AustraliaLinearResponse'
tmpDataDir = file.path(tmpBaseDir, tmpProjName)
tmpScriptDir = file.path(tmpBaseDir, tmpProjName, tmpProjName)

library(Hmisc)

### Simulate patient survival time and Recist change ###
# Should survival time and measurement change be uncorrelated?
flagCorrelatedSurvival = FALSE

# Simulated patient cohort statistics
#nSimulatedPatients = 200
meanSimSurvivalTime = 400 # days
meanSimMeasurementChange = 1 # fraction of baseline
stdSimSurvivalTime = 50
stdSimMeasurementChange = 0.2

# Arrays to store simulated survival times and measurement changes
simSurvivalTimes = c()
simMeasChange = c()

# Set the random number generator seed for reproducibility
set.seed(0)

# We simulate patient data that should have no correlation between survival and change in RECIST 
# by drawing from two normal distributions, one for determining the survival of any given patient
# and another for determining the measurement change. All patients will have only two measurements,
# one at baseline and one follow-up, to keep things simple. The measurement at follow up is always 
# 1 cm.

#dataLinRECIST has following attributes for each patient:
# meas: thickness measurements in numeric array
# baseline: baseline measurement in numeric array of same length as meas
# nadir: what measurement at any time point is the nadir, numeric array of same length as meas
# scanClass: factor with 4 levels "2", "3", "4", "6" (EG: not sure what this means)
# scanDates: Date array of format "yyyy-mm-dd", same length as meas
# changeBaseline: relative change since baseline, numeric array of length 1 less than meas
# changeNadir: relative change from nadir, numeric array of length 1 less than meas. If current measurement
#   is nadir, this is 0.
# changeMin: minimum change from baseline, single number
# SurvTimeDays: days of survival since diagnosis, single number
# SurvTimeMonths: months of survival since diagnosis (seems to be ~= SurvTimeDays / 30.4), single number
# EntryTimeDays: days of entry on clinical trial? (EG: Not sure) Single number
# EntryTimeMonths: months of entry on clinical trial? Single number

# Initialize the list containing the simulated RECIST measurements
simLinRECIST = list()

for (iPat in 1:nSimulatedPatients) {
    if (iPat == 1) {
        # Create the first row in tmpLinRECIST, which will contain the data for the simulated patients
        tmpSurvDays = rnorm(1, meanSimSurvivalTime, stdSimSurvivalTime)
        tmpMeasChange = rnorm(1, meanSimMeasurementChange, stdSimMeasurementChange)

        # Artificially insert a patient
        #tmpSurvDays = 650
        #tmpMeasChange = 0.5

        if (flagCorrelatedSurvival) {
            # Keep the two calls to rnorm above, so that random number generation is still reproducible
            # Here we make sure that tmpSurvDays and tmpMeasChange are correlated
            tmpMeasChange = 1 / (tmpSurvDays / meanSimSurvivalTime)
        }

        tmpMeasFollowUp = 1 * tmpMeasChange

        if (tmpMeasChange < 1) {
            # in this case, the follow-up is the nadir measurement
            tmpNadir = c(1, tmpMeasFollowUp)
            tmpChangeNadir = 0
        } else {
            tmpNadir = c(1, 1)
            tmpChangeNadir = tmpMeasChange - 1
        }

        tmpList = list(
            meas = c(1, tmpMeasFollowUp),
            baseline = c(1, 1),
            nadir = tmpNadir,
            scanClass = factor(c("2", "3"), c("2", "3", "4", "6")), #using "2" and "3" for no good reason
            scanDates = c("", ""),
            changeBaseline = tmpMeasChange - 1,
            changeNadir = tmpChangeNadir,
            changeMin = tmpMeasChange - 1,
            SurvTimeDays = tmpSurvDays,
            SurvTimeMonths = tmpSurvDays / 30.4,
            EntryTimeDays = tmpSurvDays,
            EntryTimeMonths = tmpSurvDays / 30.4
        )

        simLinRECIST[['1']] = tmpList
    } else {
        tmpPatIndex = toString(iPat)
        tmpSurvDays = rnorm(1, meanSimSurvivalTime, stdSimSurvivalTime)
        tmpMeasChange = rnorm(1, meanSimMeasurementChange, stdSimMeasurementChange)

        if (flagCorrelatedSurvival) {
            # Keep the two calls to rnorm above, so that random number generation is still reproducible
            # Here we make sure that tmpSurvDays and tmpMeasChange are correlated
            tmpMeasChange = 1 / (tmpSurvDays / meanSimSurvivalTime)
        }

        tmpMeasFollowUp = 1 * tmpMeasChange

        if (tmpMeasChange < 1) {
            # in this case, the follow-up is the nadir measurement
            tmpNadir = tmpMeasFollowUp
            tmpChangeNadir = 0
        } else {
            tmpNadir = 1
            tmpChangeNadir = tmpMeasChange - 1
        }

        tmpList = list(
            meas = c(1, tmpMeasFollowUp),
            baseline = c(1, 1),
            nadir = c(tmpNadir, tmpNadir),
            scanClass = factor(c("2", "3"), c("2", "3", "4", "6")), #using "2" and "3" for no good reason
            scanDates = c("", ""),
            changeBaseline = tmpMeasChange - 1,
            changeNadir = tmpChangeNadir,
            changeMin = tmpMeasChange - 1,
            SurvTimeDays = tmpSurvDays,
            SurvTimeMonths = tmpSurvDays / 30.4,
            EntryTimeDays = tmpSurvDays,
            EntryTimeMonths = tmpSurvDays / 30.4
        )

        simLinRECIST[[tmpPatIndex]] = tmpList
    }
    simSurvivalTimes[iPat] = tmpSurvDays
    simMeasChange[iPat] = tmpMeasChange
}

# set up the continuous metrics
tmpColsOfInterest = sort(
    c(
        match("URN", names(simDataFullCols)),
        match("Status", names(simDataFullCols)),
        match("SurvTimeDays", names(simDataFullCols)),
        match("EntryTimeDays", names(simDataFullCols)),
        match("BaseTimeDays", names(simDataFullCols)),
        match("RecistChange2", names(simDataFullCols)),
        match("RecistChangeN", names(simDataFullCols))
        )
    )

# Generate the dataFullCols and dataContinuous data frames for the simulated data
tmpDataFullRow = simDataFullCols[1,]
tmpDataFullRow['DateDiagnosis'] = ''
tmpDataFullRow['DateDeath'] = ''
tmpDataFullRow['EntryDate'] = ''
tmpDataFullRow['A_C1D1.date'] = ''
tmpDataFullRow['A_C2D1.date'] = ''
tmpDataFullRow['A_C3D1.date'] = ''
tmpDataFullRow['A_C4D1.date'] = ''
tmpDataFullRow['A_C5D1.date'] = ''
tmpDataFullRow['A_C6D1.date'] = ''
tmpDataFullRow['APetDate'] = ''
tmpDataFullRow['Date1'] = ''
tmpDataFullRow['Date2'] = ''
tmpDataFullRow['Date3'] = ''
tmpDataFullRow['Date4'] = ''
tmpDataFullRow['DateN'] = ''
tmpDataFullRow['Weight'] = 0
tmpDataFullRow['Height'] = 0
tmpDataFullRow['AgeDeath'] = 0
tmpDataFullRow['AgeEntry'] = 0
tmpDataFullRow['AgeDiagnosis'] = 0
tmpDataFullRow['BaseTimeDays'] = 0
tmpDataFullRow['Duration2'] = 0
tmpDataFullRow['Duration3'] = 0
tmpDataFullRow['Duration4'] = 0

tmpDataContinuousRow = simDataFullCols[1, tmpColsOfInterest]
tmpDataContinuousRow$RecistBest = NA
simDataFullCols = tmpDataFullRow
simDataContinuous = tmpDataContinuousRow
for (iPat in 1:nSimulatedPatients) {
    # simulated dataFull
    simDataFullCols[iPat,] = tmpDataFullRow
    simDataFullCols$URN[iPat] = iPat
    simDataFullCols$Status[iPat] = 1
    simDataFullCols$SurvTimeDays[iPat] = simLinRECIST[[as.character(iPat)]]$SurvTimeDays
    simDataFullCols$EntryTimeDays[iPat] = simLinRECIST[[as.character(iPat)]]$EntryTimeDays
    simDataFullCols$BaseTimeDays[iPat] = 30 # EG: Not sure what "BaseTimeDays" entails, putting in 30
    simDataFullCols$RecistChange2[iPat] = simLinRECIST[[as.character(iPat)]]$changeBaseline
    simDataFullCols$RecistChangeN[iPat] = simLinRECIST[[as.character(iPat)]]$changeBaseline
    simDataFullCols$Recist1[iPat] = simLinRECIST[[as.character(iPat)]]$meas[1]
    simDataFullCols$Recist2[iPat] = simLinRECIST[[as.character(iPat)]]$meas[2]

    # simulated dataContinuous
    simDataContinuous[iPat,] = tmpDataContinuousRow
    simDataContinuous$URN[iPat] = iPat
    simDataContinuous$Status[iPat] = 1
    simDataContinuous$SurvTimeDays[iPat] = simLinRECIST[[as.character(iPat)]]$SurvTimeDays
    simDataContinuous$EntryTimeDays[iPat] = simLinRECIST[[as.character(iPat)]]$EntryTimeDays
    simDataContinuous$BaseTimeDays[iPat] = 30 # EG: Not sure what "BaseTimeDays" entails, putting in 30
    simDataContinuous$RecistChange2[iPat] = simLinRECIST[[as.character(iPat)]]$changeBaseline
    simDataContinuous$RecistChangeN[iPat] = simLinRECIST[[as.character(iPat)]]$changeBaseline
    simDataContinuous$RecistBest[iPat] = simLinRECIST[[as.character(iPat)]]$changeMin
}

rm(list = ls(pattern = 'tmp.*'))
rm(list = ls(pattern = 'i[[:upper:]]'))