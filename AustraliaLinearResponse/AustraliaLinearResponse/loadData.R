# sometimes I just want to use the data that already exists
# instead of re-calculating it all over again!

tmpBaseDir=file.path('I:','MesoResponse')
tmpProjName='AustraliaLinearResponse'
tmpDataDir=file.path(tmpBaseDir,tmpProjName)
tmpScriptDir=file.path(tmpBaseDir,tmpProjName,tmpProjName)

# import the basic data
source(file.path(tmpScriptDir,'importData.R'))

tmpBaseDir=file.path('I:','MesoResponse')
tmpProjName='AustraliaLinearResponse'
tmpDataDir=file.path(tmpBaseDir,tmpProjName)
# get the regular optimization data
load(file.path(tmpDataDir,'FullDataSetPerformance.rdata'))

# get the cross-validation data
load(file.path(tmpDataDir,'CVDataSetPerformance.rdata'))

rm(list=ls(pattern='tmp.*'))