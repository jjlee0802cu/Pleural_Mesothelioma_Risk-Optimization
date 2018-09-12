rm(list = ls())

# Path to Project directory
tmpDataDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern', 'AustraliaLinearResponse')
tmpDataFile = file.path(tmpDataDir, 'fromR_cohort3D_anonymized.rdata')
load(tmpDataFile)


rm(list = ls(pattern = 'tmp.*'))
rm(list = ls(pattern = 'i[[:upper:]]'))