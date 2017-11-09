# hst-project
Note: This file goes over the basic reduction procedure. For more details about individual scripts, see their .readme files on Github (currently being written).

First, obtain the data from MAST (https://archive.stsci.edu/hst/). The directory structure should look like:  

ProposalID/GalaxyID/raw/ ->>>raw files go in here  
ProposalID/GalaxyID/reduced/ ->>>cphot/runphot/runfake files go in here  

Then, we need to create input files for DOLHPOT:  

-auto_runphot.py will create the necessary runphot file. auto_runphot.py should itself be placed in ProposalID/GalaxyID/auto_runphot.py and ran from there.  
-The cphot file containing the dithering patterns needs to be created manually.  
-runfake2.0 can be copied over from any other galaxy directory.  

Run the photometry and complete completeness simulations within DOLPHOT from the command line:  
-To run the photometry, do the following from within ../reduced/ (this is an example, items in parenthesis give description of what to input): ./runphot5 ngc5128-s1(TARGET)  
-To run the completeness simulations from within ../reduced/ (this is an example, items in parenthesis give description of what to input): nohup ./runfake CENA-132557(TARGET)  

You should now have many files within ../reduced. The ones we care about (right now) are GalaxyID.phot and GalaxyID.phot.fake. Bring these files out to a separate working directory (currently I am doing the rest of the work on my own computer, but ideally this would be done in ProposalID/GalaxyID/TRGB/ or something similar).  

In /TRGB/, along with the two photometry files, make sure you have both_phot_reader.ipynb (Github), which reads in the photometry and completeness files and returns two files that have certain photometric criterion applied, as well as column headers added (they will be named GalaxyID.phot2 and GalaxyID.phot.fake2).  

You may want to apply a spatial selection to the data- to do this, use the interactiveVI.ipynb (note this is under construction and runs a bit slow through Jupyter- for now, please use the regular .py script). The coordinates will be saved in a text file.  

Next comes the running of the TRGB notebook (TRGB_GA.ipynb) to get the TRGB magnitude, errors, etc. If you want to apply a spatial selection, uncomment the relevant section in the second cell of the notebook and insert the coordinates obtained from the spatial selection tool. You can then run the notebook one cell at a time, or the entire thing at once.  

-One important thing to make sure things go smoothly is to adjust the location of the dolphot installation in the autorunphot script as well as the cphot and runfake files- the automatic detection of the correct location is currently bugged. 

And that should be it! (for now...)
