# sometimes I just want to use the data that already exists
# instead of re-calculating it all over again!

tmpBaseDir=file.path('F:','MesoResponse')
tmpProjName='AustraliaLinearResponse'
tmpDataDir=file.path(tmpBaseDir,tmpProjName)
tmpScriptDir=file.path(tmpBaseDir,tmpProjName,tmpProjName)

# import the basic data
source(file.path(tmpScriptDir,'importData.R'))

tmpBaseDir=file.path('F:','MesoResponse')
tmpProjName='AustraliaLinearResponse'
tmpDataDir=file.path(tmpBaseDir,tmpProjName)
# get the regular optimization data
load(file.path(tmpDataDir,'FullDataSetPerformanceEntry.rdata'))

# get the cross-validation data
load(file.path(tmpDataDir,'CVDataSetPerformanceEntry.rdata'))

rm(list=ls(pattern='tmp.*'))