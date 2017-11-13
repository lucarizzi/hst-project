both_phot_reader is a python program that will take the output files from DOLPHOT photomtery and give out photometry files that are ready for use with the TRGB tool. Along with converting the file to the correct explicit format, it also applies a set of selection criteria based on things such as chi^2, sharpness, etc.

The only parameter that needs to be changed is the prefix, which corresponds to the file name.

Notes:

-both_phot_reader_LR is working and tested for those observations with 2 images in each filter. Use this if this is your working case.

-both_phot_readerAuto is developed to handle instances where there are more than 2 images in each filter (often the case with archival proposals not meant to measure the TRGB). This has been tested lightly and shown to be consistent with others, but could use more scrutiny (upcoming).
