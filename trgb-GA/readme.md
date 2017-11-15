Now we are ready to run the TRGB tool. This will take in the .phot2 and .phot.fake2 files from both_phot. It will also require some user input and tinkering with certain parameters in order to achieve a stable fit.

In the second cell, change the prefix to match the input files. Stick with the current binning (30), as changing it is unlikely to change the TRGB magnitude or your fit (though you can try if you can't get a good fit). Also, input an initial guess for the TRGB magnitude. If you are not working with F606W and F814W, change band1 and band2 in the next cell. Also in that cell, you will see that it takes in the coordinate.txt file from the interactive_V-I script. Make sure that is in the directory (if you want to select the whole region, you can still do that, but you need the file there still). 

In the 8th cell, you should input the size and location of your box used to fit the TRGB. This is in the line 
PLT ={'cl' : 1.1, 'ch' : 1.8, 'ml' : 23.6, 'mh' : 25.3} (numbers are examples currently), where cl and ch correspond to low and high colors, and ml and mh correspond to low and high magnitudes. 

In cell 13, you must input your initial guesses for a, b, and c, as well as limits for each of these parameters. This may take some time to adjust properly to get a stable fit, depending on your CMD. 

And that is all! There isn't need to change any of the rest of the notebook, though you may adjust which plots and saved to files/how they look.



Notes:

-Use the TRGB_GA_example.ipynb to begin your analysis, and adjust the parameters listed above.


