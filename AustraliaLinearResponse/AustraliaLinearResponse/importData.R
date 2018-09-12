rm(list=ls())
tmpDataDir=file.path('C:','Users','Justin','Documents','_UC intern')
#
# diff. cohorts:
#3a - 2 or more scans with RECIST, baseline not necessary
#3b - 2 or more scans with RECIST, baseline has to be present
#3c - same as 3a, except omit patients with no measureable disease
#3d - same as 3b, except omit patients with no measurebale disease
tmpDataFile=file.path(tmpDataDir,'fromR_cohort3D.rdata')
load(tmpDataFile)



rm(list=ls(pattern='tmp.*'))
rm(list=ls(pattern='i[[:upper:]]'))