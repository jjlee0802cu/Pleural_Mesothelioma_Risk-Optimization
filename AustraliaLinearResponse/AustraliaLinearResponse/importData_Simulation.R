# Script to read in anonymized, simulated data
#
# Written by E Gudmundsson, June 2017

rm(list = ls())

# Path to Project directory
tmpDataDir = file.path('C:', 'Users','Justin','Documents','_UC intern','AustraliaLinearResponse')
tmpDataFile = file.path(tmpDataDir, 'SimulatedData_100Patients.rdata')
load(tmpDataFile)


rm(list = ls(pattern = 'tmp.*'))
rm(list = ls(pattern = 'i[[:upper:]]'))