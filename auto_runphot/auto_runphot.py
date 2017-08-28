import os, shutil, sys
import fnmatch
import glob
from astropy.io import fits
import pandas as pd

# determine where dolphot lives
dolphotDir = "/".join(shutil.which('dolphot').split('/')[:-1])
if dolphotDir == "":
	sys.stdout.write("DOLPHOT is not in your path. Aborting.\n")
	sys.exit(1)

#set directory for raw data
directory='./raw/'

#create list of flc files in directory
files = [str(x) for x in fnmatch.filter(os.listdir(directory), '*flc.fits')]

#find base name
base_name = os.path.commonprefix(files)
if base_name == "":
	sys.stdout.write("I could not find a common prefix\n")
	sys.exit(1)
else:
	sys.stdout.write("The common prefix is %s\n" % (base_name))

# create a pandas structure
mainDirectory = pd.DataFrame(columns = ['file','filter'])
mainDirectory['file'] = files

# use the file name as index
mainDirectory = mainDirectory.set_index('file')

# parse files to find their filter and add to the pandas structure
for fitsName,row in mainDirectory.iterrows():
    hdulist = fits.open('raw/'+fitsName)
    instrument = hdulist[0].header['INSTRUME']
    if instrument == "WFC3":
        filter = hdulist[0].header['FILTER']
    if instrument == "ACS":
        filter = hdulist[0].header['FILTER1']
        if filter == "" or filter == "CLEAR1L":
            filter = hdulist[0].header['FILTER2']
    mainDirectory.set_value(fitsName,'filter',filter)
    hdulist.close()

# define a filter translation function
def filter_translation(filter):
    if filter in ['F555W','F606W']:
            return "V"
    if filter in ['F814W']:
            return "I"
    return "indef"
# and apply it
mainDirectory['converted'] =mainDirectory.apply(lambda x: filter_translation(x['filter']), axis=1)

# add a counter for the V and I images
vCounter = 1
iCounter = 1
for fitsName,row in mainDirectory.iterrows():
    if mainDirectory.get_value(fitsName,'converted')=='V':
        mainDirectory.set_value(fitsName,'index',vCounter)
        vCounter +=1
    if mainDirectory.get_value(fitsName,'converted')=='I':
        mainDirectory.set_value(fitsName,'index',iCounter)
        iCounter +=1
mainDirectory['index']=mainDirectory['index'].astype(int)


#create text file
f = open('./reduced/runphot5','w')

#write intro 
f.write('#!/bin/tcsh\n')
#f.write('setenv BASE $1\n')
f.write('setenv TARG $1\n')
f.write('setenv DOLPHOT_DIR %s\n' % (dolphotDir))


#first paragraph
for fitsfile,row in mainDirectory.iterrows():
	f.write('\ncp ../raw/%s ${TARG}_%s%d.fits' % (fitsfile,row['converted'],row['index']))
f.write('\n')

# second paragraph 
for fitsfile,row in mainDirectory.iterrows():
 	if instrument == "ACS":
 		f.write('\nnice +19 ${DOLPHOT_DIR}/acsmask ${TARG}_%s%d.fits' % (row['converted'], row['index']))
 	if instrument == "WFC3":
 		f.write('\nnice +19 ${DOLPHOT_DIR}/wfc3mask ${TARG}_%s%d.fits' % (row['converted'], row['index']))

f.write('\n')

# third paragraph
for fitsfile,row in mainDirectory.iterrows():
 	f.write('\nnice +19 ${DOLPHOT_DIR}/calcsky ${TARG}_%s%d 15 35 -128 2.25 2.00'% (row['converted'], row['index']))

f.write('\n')

# write end and close
if instrument == "ACS":
 	f.write('\n\nnice +19 ${DOLPHOT_DIR}/dolphot ${TARG}.phot -pcphot5.param\n')
if instrument == "WFC3":
 	f.write('\n\nnice +19 ${DOLPHOT_DIR}/dolphot ${TARG}.phot -pcphotwfc3.param\n')
#f.write("cat ${TARG}.phot | awk ''$5<=2.5 && $7*$7<=0.09 && $11<=2 && $20>=5 && $24==0 && $33>=5 && $37==0'' > ${TARG}.phot2")
f.close()
