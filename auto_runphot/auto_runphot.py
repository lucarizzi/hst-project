import os
import fnmatch
import glob
from astropy.io import fits

#set directory for raw data
directory='./raw/'

#create list of flc files in directory
files = fnmatch.filter(os.listdir(directory), '*flc.fits')
print files

#get number of flc files dealing with
num = len(files)/2 #change later so it isn't just limited to 2 filters?
print num

#find base name
base_name = os.path.commonprefix(files)
#print base_name


#below is to get filter information  
filter_array=[]
for fitsName in glob.glob('./raw/*flc.fits'):
	hdulist = fits.open(fitsName)
	filters = hdulist[0].header['FILTER1']
	filter_array.append(filters)
	filters = hdulist[0].header['FILTER2']
	filter_array.append(filters)
	if 'CLEAR1L' in filter_array: filter_array.remove('CLEAR1L')
	if 'CLEAR2L' in filter_array: filter_array.remove('CLEAR2L')
	print filter_array
	hdulist.close()

#following manual bit is for purposes of doing it locally
#filter_array = ['F814W', 'F814W', 'F814W', 'F606W', 'F606W', 'F606W']

#create new array to then modify with V, I, etc.
filter_array_short = filter_array

for (i, item) in enumerate(filter_array_short):
    if item == 'F606W':
        filter_array_short[i] = 'V'
    if item == 'F814W':
        filter_array_short[i] = 'I'
print filter_array_short


#create text file
f = open('./reduced/runphot5','w')

#write intro 
f.write('#!/bin/tcsh\n')
f.write('setenv BASE $1\nsetenv TARG $2\n')
f.write('setenv DOLPHOT_DIR /Volumes/External_lr/HST/dolphot2.0/bin\n')

#get rid of base
files = [elem[len(base_name):] for elem in files]
#print file_list

#first paragraph
count=1
count2=1
while (count < 2*num+1):
	f.write('\ncp ../raw/${BASE}%s ${TARG}_%s%s.fits' % (files[count-1], filter_array[count-1], count2))
	count=count+1
	count2=count2+1
	if count2>num:
		count2=1

f.write('\n')

#second paragraph 
count=1
count2=1
while (count < 2*num+1):
	f.write('\nnice +19 ${DOLPHOT_DIR}/acsmask ${TARG}_%s%s.fits' % (filter_array[count-1], count2))
	count=count+1
	count2=count2+1
	if count2>num:
		count2=1

f.write('\n')

#third paragraph
count=1
count2=1
while (count < 2*num+1):
	f.write('\nnice +19 ${DOLPHOT_DIR}/calcsky ${TARG}_%s%s 15 35 -128 2.25 2.00'% (filter_array[count-1], count2))
	count=count+1
	count2=count2+1
	if count2>num:
		count2=1

f.write('\n')

#write end and close
f.write('\n\nnice +19 ${DOLPHOT_DIR}/dolphot ${TARG}.phot -pcphot5.param\n')
f.write('cat ${TARG}.phot | awk ''$5<=2.5 && $7*$7<=0.09 && $11<=2 && $20>=5 && $24==0 && $33>=5 && $37==0'' > ${TARG}.phot2')
f.close()
