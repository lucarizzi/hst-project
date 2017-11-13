This python procedure allows you to interactively select a (rectangular) region of an image to perform a fit to. This is useful for situations where the galaxy does not take up the whole field, as well as those where you would like to avoid a selection of blue stars (spiral arms for instance). This thereby reduces contamination. 

To use this tool, simply go to your command line and type:
python interactive.py

The script will ask you for the name of the galaxy- this is just asking for the prefix of the .phot2 file (the photometry file created from both_phot). After reading this file in, a matplotlib window will pop up, which will plot all stars in the image using their x and y positions, and smoothing them with a gaussian kernel to highlight dense regions. You should first click inside the window once, and from then on you can select your rectangular region. Note that you must select this region from the top left to the bottom right, as it outputs the coordinates in a certain order which must be preserved for the TRGB tool. 

You can reselect as many times as you want until you are happy with your selection (which is being outputted to coordinates.txt). Once you are content, close the matplotlib window, and a CMD with the selected stars will appear. You can rerun the script if you still wish to change the selection. The coordinates file The python TRGB tool then takes this text file and uses it to select stars (as can be seen in the third cell of the TRGB jupyter notebook).


Important Notes:

-Make sure you have annotate.py in your directory.

-The coordinates.txt file is rewritten with every selection, so make sure you double-check the contents of that file before proceeding with your analysis.

-Make sure to select from top-left to bottom-right!

-I am currently working on updating this to a free-form selection tool. 
