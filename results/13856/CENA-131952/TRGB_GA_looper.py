# coding: utf-8
#import modules
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import ndimage
from scipy.interpolate import BSpline
from scipy.interpolate import interp1d
from scipy.interpolate import spline
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
# get_ipython().magic('matplotlib inline')


f = open('./TRGB_loopTest_writeout','w')

iLoop=21
while iLoop < 51:
    iLoopString = str(iLoop)

    #edit below parameters
    prefix = 'CENA-131952' 
    PLT ={'cl' : 0.9, 'ch' : 1.9, 'ml' : 23, 'mh' : 25 }

    #band names
    band1 ='F606W'
    band2 ='F814W'

    #read in files and create pandas dataframes with header names from first row
    photometryFile = "%s.phot2" % prefix
    completenessFile = "%s.phot.fake2" % prefix
    photDF = pd.read_csv(photometryFile, delim_whitespace=True, header=0) 
    completeDF = pd.read_csv(completenessFile, delim_whitespace=True, header=0) 

    #spatial selection on photoDF - get low & high values from interactive tool
    coordFile = open('coordinates.txt')
    for line in coordFile:
        coords = line.strip().split()
        xLow = int(coords[0])
        xHigh = int(coords[1])
        yHigh = int(coords[2])
        yLow = int(coords[3])
    coordFile.close()

    photDF = photDF[(photDF['x'] > xLow) & (photDF['x'] < xHigh) 
            & (photDF['y'] > yLow) & (photDF['y'] < yHigh)]

    #rename completeDF columns to match other script (for ease)
    completeDF.columns = ['stV', 'stI', 'fakeV', 'fakeI', 'fakeVerr', 'fakeIerr']
    #completeDF[0:3]

    #set 99.999 values to nan values in fakeV and fakeI of completeDF
    completeDF.loc[completeDF.fakeV == 99.999, 'fakeV'] = np.nan
    completeDF.loc[completeDF.fakeI == 99.999, 'fakeI'] = np.nan
    completeDF[13:19]

    #convert to numpy arrays for ease of math
    stV = completeDF.stV.values
    stI = completeDF.stI.values
    fakeV   = completeDF.fakeV.values
    fakeI   = completeDF.fakeI.values
    fakeVerr   = completeDF.fakeVerr.values
    fakeIerr   = completeDF.fakeIerr.values

    #also convert some for real photometry file
    realV = photDF.mag_V
    realVErr = photDF.Verr
    realI = photDF.mag_I
    realIErr = photDF.Ierr

    #additional selection criteria for fake photometry
    #SHOULD LAST TWO BE FAKEIERR INSTEAD?
    fakeSelectionParameter = 0.3
    stV_sel = stV[np.logical_and(fakeVerr < fakeSelectionParameter, fakeV != np.nan)]
    fakeV_sel   = fakeV[np.logical_and(fakeVerr < fakeSelectionParameter, fakeV != np.nan)]
    stI_sel = stI[np.logical_and(fakeVerr < fakeSelectionParameter, fakeI != np.nan)]
    fakeI_sel   = fakeI[np.logical_and(fakeVerr < fakeSelectionParameter, fakeI != np.nan)]


    # In[ ]:

    #define a parameter to be used in GetSmoothedStats (related to bin size, based on # of stars)
    #basedOnStarCounts = int((0.005) * len(photDF)) #coefficient could probably use testing
    #print('basedOnStarCounts:',basedOnStarCounts)


    # In[ ]:

    #define GetSmoothedStats
    def GetSmoothedStats(Xall,X,Y):
        medx = np.median(X)
        w=np.median([abs(x-medx) for x in X])/0.6745
        w=w*(4/(3*len(X)))**(1/5)
        w=max(w/4,0.005)
        #xi = np.linspace(min(X),max(X),basedOnStarCounts)
        xi = np.linspace(min(X),max(X),iLoop)
        xi=np.append(xi,np.max(xi)+0.5*w)
        xi=np.append(xi,np.max(xi)+w)

        dX=Y-X
        nall = np.zeros(len(xi))
        n = nall.copy()
        b = nall.copy()
        s = nall.copy()
        for k in range(len(xi)):
            nall[k]=np.sum(np.e**(-0.5*((xi[k]-Xall)/w)**2))/(w*np.sqrt(2*np.pi))
            p = np.e**(-0.5*((xi[k]-X)/w)**2) / (w*np.sqrt(2*np.pi))
            n[k]   =np.sum(np.e**(-0.5*((xi[k]-X   )/w)**2))/(w*np.sqrt(2*np.pi))
            b[k]=np.sum(dX*p) / n[k]
            s[k]=np.sqrt(np.sum(dX**2*p)/n[k]-b[k]**2)
            if n[k]>=2:
                s[k]=s[k]*np.sqrt(n[k]/(n[k]-1))
        
        c=n.copy()/nall.copy()
        ind = np.logical_and(np.isfinite(c),np.isfinite(b),np.isfinite(s))
        c=c[ind]
        b=b[ind]
        s=s[ind]
        xi=xi[ind]
        return xi,c,b,s


    # In[ ]:

    #calculate 4th order coefficients to be used in error, completeness, and bias functions
    deg=4
    mag_V,complet_V,bias_V,errors_V = GetSmoothedStats(stV,stV_sel,fakeV_sel)
    mag_I,complet_I,bias_I,errors_I = GetSmoothedStats(stI,stI_sel,fakeI_sel)
    centerCoef = np.polyfit(mag_I, bias_I, deg=deg)
    errCoef = np.polyfit(mag_I, errors_I, deg=deg)
    compltCoef = np.polyfit(mag_I,complet_I, deg=deg)


    # In[ ]:

    #create plots of error, completeness, and bias functions
    fig = plt.figure(figsize=(10,7))

    PLTCenters = fig.add_subplot(131)
    PLTCenters.plot(mag_I,bias_I)
    PLTCenters.set_xlabel(band2)
    PLTCenters.set_ylabel("mag")
    PLTCenters.set_title("Centers")
            
    PLTErrors= fig.add_subplot(132)
    PLTErrors.plot(mag_I,errors_I,color='blue')
    PLTErrors.set_xlabel(band2)
    PLTErrors.set_ylabel("mag")
    PLTErrors.set_title("Errors")

    Completeness = fig.add_subplot(133)
    Completeness.plot(mag_I,complet_I)
    Completeness.set_xlabel(band2)
    Completeness.set_ylabel("Fraction")
    Completeness.set_title("Completeness")

    fig.tight_layout()
    #plt.savefig('%sfunctions.ps' %iLoopString, dpi=300, orientation='landscape')


    # In[ ]:

    #CMD with boxed region
    plt.rc('font',family='serif')
    plt.rc('xtick',labelsize='12')
    plt.rc('ytick',labelsize='12')
    cmd = plt.figure(figsize=(6,8))
    ax = cmd.add_subplot(1,1,1)
    ax.plot(realV-realI, realI,'o',markersize=1, color='black')
    # PLT ={'cl' : 1.0, 'ch' : 2.3, 'ml' : 23, 'mh' : 25 }
    ax.plot( [ PLT['cl'], PLT['cl']], [PLT['ml'], PLT['mh'] ] , '--', color='red')
    ax.plot( [ PLT['ch'], PLT['ch']], [PLT['ml'], PLT['mh'] ] , '--', color='red')
    ax.plot( [ PLT['cl'], PLT['ch']], [PLT['ml'], PLT['ml'] ] , '--', color='red')
    ax.plot( [ PLT['cl'], PLT['ch']], [PLT['mh'], PLT['mh'] ] , '--', color='red')
    ax.set_xlabel("%s - %s" %(band1,band2), fontsize=16)
    ax.set_ylabel(band2, fontsize=16)
    VRange=photDF.mag_V.values
    IRange=photDF.mag_I.values
    ax.set_xlim((np.amin(VRange-IRange)+1,np.amax(VRange-IRange)-1))
    ax.set_ylim((np.amax(IRange),np.amin(IRange)+3))
    cmd.tight_layout()
    # plt.savefig('cmd.ps', dpi=300)


    # In[ ]:

    #luminosity function
    binsize = 0.1
    bin = np.arange( PLT['ml']-0.5, PLT['mh']+0.5, binsize)
    magtics = bin[1:] - binsize/2.
    LF, binedeges = np.histogram([i for v,i in zip(realV,realI) 
                                  if (v-i>=PLT['cl'] 
                                    and v-i<=PLT['ch'] 
                                    and i>=PLT['ml'] 
                                    and i<=PLT['mh'])],bins=bin)


    # In[ ]:

    #smoothed luminosity function
    m = np.linspace(PLT['ml']-0.5,PLT['mh']+0.5,50)
    dm = m[1]-m[0]
    smoothLF=np.zeros(len(m))
    for k in range(len(smoothLF)):
        smoothLF[k]=np.sum([(1/(np.sqrt(2*np.pi)*sigma)*np.e**(-1.0*((mag-m[k])**2/(2*sigma**2)))) 
                            for mag,sigma in zip(realI,realIErr)])

    # calculate integral under the curve
    normalization = np.sum(smoothLF)*dm
    smoothLF = smoothLF/normalization
    total = np.sum(smoothLF)

    # reduce normalization
    normalization = normalization / (dm*len(m))

    # SOBEL FILTERING
    LF_sobel = ndimage.sobel(smoothLF)


    # In[ ]:

    #luminosity function and sobel plot
    fig=plt.figure()
    LF_prime = np.diff(LF) / binsize
    magtics_prime = bin[1:-1]

    # smoothed version built by interpolating (spline)
    magtics_interp = np.arange(min(bin),max(bin), 0.01)
    # LF_interp = interp1d(magtics, LF, kind='cubic')
    # LF_prime_interp = interp1d(magtics_prime, LF_prime, kind='cubic')
    LF_interp = spline(magtics, LF, magtics_interp)
    LF_prime_interp = spline(magtics_prime, LF_prime, magtics_interp)

    majorLocator = MultipleLocator(1)
    minorLocator = MultipleLocator(0.1)

    LFplt = fig.add_subplot(211)
    LFplt.plot(magtics_interp,LF_interp,color='blue')
    LFplt.set_title("Interpolated")
    ymaxLF = np.amax(LF)
    ymaxLFprime = np.amax(LF_prime_interp)

    # rescale LFPRIME
    LF_prime_interp = LF_prime_interp / ymaxLFprime * ymaxLF
    LFplt.plot(magtics_interp,LF_prime_interp,':',color='green')
    LFprimemax = np.amax(LF_prime_interp)
    LFplt.xaxis.set_major_locator(majorLocator)
    LFplt.xaxis.set_minor_locator(minorLocator) 
    LFplt.set_xlim((PLT['ml'],PLT['mh']))
    LFplt.set_ylabel("N")
    ymin,ymax = LFplt.get_ylim()
    ymax = np.amax(LF)
    LFplt.set_ylim((0,ymax))
    LFsmoothplt = fig.add_subplot(212)
    fig.subplots_adjust(hspace=.4)
    LFsmoothplt.plot(m,smoothLF)

    # rescale sobel filter
    ymaxLF = np.amax(smoothLF)
    ymaxsobel = np.amax(LF_sobel)
    LF_sobel = LF_sobel/ymaxsobel*ymaxLF
    LFsmoothplt.plot(m,LF_sobel,':',color='green')
    LFsmoothplt.set_title("Smoothed")
    LFsmoothplt.xaxis.set_major_locator(majorLocator)
    LFsmoothplt.xaxis.set_minor_locator(minorLocator)
    LFsmoothplt.set_xlabel(band2)
    LFsmoothplt.set_ylabel("N")
    ymin,ymax = LFsmoothplt.get_ylim()
    LFsmoothplt.set_xlim((PLT['ml'],PLT['mh']))
    LFsmoothplt.set_ylim((0,ymax))
    #plt.savefig('%sLFsmooth.ps' %iLoopString, dpi=300)

    pmagmax = np.amax(realI) #max(x[1] for x in parr)
    pmagmin = np.amin(realI) #min(x[1] for x in parr)
    fmagmax = np.amax(stV) #max(x[1] for x in farr)
    fmagmin = np.amin(stI) #min(x[1] for x in farr)
    deg = 4
    np.random.seed(0)


    # In[ ]:

    #define several needed functions
    def idealLF(m,a,b,c,mTRGB):
            if m >= mTRGB:
                return 10**(a*(m-mTRGB)+b)
            else :
                return 10**(c*(m-mTRGB))

    def Compltfunc(m):
            return compltCoef[0]*m**4 + compltCoef[1]*m**3 + compltCoef[2]*m**2 + compltCoef[3]*m**1 + compltCoef[4]
            #return 1.0

    def mbias(m):
            return centerCoef[0]*m**4 + centerCoef[1]*m**3 + centerCoef[2]*m**2 + centerCoef[3]*m**1 + centerCoef[4]+m
            #return 0.001
    def merr(m):
            return errCoef[0]*m**4 + errCoef[1]*m**3 + errCoef[2]*m**2 + errCoef[3]*m**1 + errCoef[4]    
            #return 0.01
    def ErrFunc(m, m_prime):
            return (1 / (np.sqrt(2*np.pi)* merr(m_prime) )) * np.e**( - ( m - mbias(m_prime) )**2 / ( 2 * merr(m_prime)**2) )    

    def limitFunc(a,b,c,mTRGB):
            " If parameter(s) is out of allowed range, return a big number and add it to LF, preventing fitting procedure choose parameters user doesn't want"
            if a<LIM['al'] or a>LIM['ah'] or b<LIM['bl'] or b>LIM['bh'] or c<LIM['cl'] or c>LIM['ch'] or mTRGB<LIM['mTRGBl'] or mTRGB>LIM['mTRGBh']: return 1.e100
            else: return 0    
        
    def realLF(mlist, a, b, c, mTRGB):  
            binsize = mlist[1]-mlist[0]
            print(a,b,c, mTRGB)
            rLF=[]
            for m in mlist:
                rLF.append(integrate.quad(lambda m_prime: idealLF(m_prime,a,b,c,mTRGB)*Compltfunc(m_prime)*ErrFunc(m,m_prime),pmagmin, pmagmax)[0])
                #rLF.append(integrate.quad(lambda m_prime: idealLF(m_prime,a,b,c,mTRGB)*ErrFunc(m,m_prime),pmagmin, pmagmax)[0])
            return np.asarray(rLF)/(sum(rLF)*binsize) + limitFunc(a,b,c,mTRGB)


    # In[ ]:

    #define arrays for TRGB measurements
    INIT  = {'TRGB' : 99, 'a' : 0.3, 'b' : 0.3, 'c' : 0.3 }
    LIM   = {'TRGBl' : -99, 'TRGBh' : 99, 'al' : -99, 'ah' : 99, 'bl' : -99, 'bh' : 99, 'cl' : -99, 'ch' : 99 }
    RANGE = { 'l' : -1., 'h' : 1. }

    #establish fit guesses
    INIT['TRGB'] = 23.9
    INIT['a'] = 0.25
    INIT['b'] = 0.25
    INIT['c'] = 0.25

    LIM['mTRGBl'] = 22
    LIM['mTRGBh'] = 25
    LIM['al'] = 0
    LIM['ah'] = 5
    LIM['bl'] = 0
    LIM['bh'] = 5
    LIM['cl'] = 0
    LIM['ch'] = 5

    magnitudes_I=[]
    magnitudes_VI=[]
    for x in zip(realV,realI):
        if ( x[0]-x[1] >= PLT['cl'] and x[0]-x[1] <= PLT['ch'] ):
            magnitudes_I.append(x[1])
            magnitudes_VI.append([x[0],x[1]])
    binsize=0.1
    bins = np.arange(INIT['TRGB']+RANGE['l'], INIT['TRGB']+RANGE['h'], binsize)
    obsLF, bin_edges = np.histogram(magnitudes_I, bins=bins, density=True)
    somma = np.sum(obsLF)
    #print(somma)
    bin_centers = (bin_edges[:-1]+bin_edges[1:])/2.
    smoothed_bins = np.arange(min(bin_edges),max(bin_edges)-binsize/2,0.02)
    smoothed_obsLF = spline(bin_centers,obsLF,smoothed_bins)
    smoothed_obsLF = spline(bin_centers,obsLF,smoothed_bins)
    plt.plot(smoothed_bins,smoothed_obsLF)
    plt.plot(bin_centers,obsLF,color='green')
    #print(np.amax(obsLF))


    # In[ ]:

    #TRGB best fit information 
    varini=[ INIT['a'], INIT['b'], INIT['c'], INIT['TRGB'] ]
    TRGBcolor = [0., 0., 0.,]
    BFIT, pcov, info, mesg, ierr = curve_fit(realLF, m, smoothLF, p0=varini,full_output=True)
    # print("\nBest fit", "TRGB: ",BFIT[3],"\na: ",BFIT[0], "\nb: ",BFIT[1],"\nc: ",BFIT[2])
    perr=np.sqrt(np.diag(pcov))


    # In[ ]:

    #luminosity function plot
    initial_guess = realLF(m,INIT['a'],INIT['b'],INIT['c'],INIT['TRGB'])
    calculated_function = realLF(m,a=BFIT[0],b=BFIT[1],c=BFIT[2],mTRGB=BFIT[3])
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(1,1,1)
    ax.plot(m,smoothLF, label='LF')
    ax.plot(m,initial_guess,color="green", label='Initial guess')
    ax.plot(m,calculated_function,color='magenta', label='Best fit')
    legend = ax.legend(loc='upper left', fontsize=14)
    #plt.savefig('%slf_fit.ps' %iLoopString, dpi=300)


    # In[ ]:

    #define functions for best fits
    def bestidealLF( mlist, a, b, c, mTRGB):
            iLF=[]
            for m in mlist:
                if m >= mTRGB: iLF.append( 10**(a*(m-mTRGB)+b ))
                else: iLF.append( 10**(c*(m-mTRGB)) )
            return np.asarray(iLF)
     
    def bestrealLF( mlist, a, b, c, mTRGB):
            rLF=[]
            for m in mlist:
                rLF.append(integrate.quad(lambda m_prime: idealLF(m_prime,a,b,c,mTRGB)*Compltfunc(m_prime)*ErrFunc(m,m_prime),pmagmin, pmagmax)[0])
            return np.asarray(rLF)
        
    bestiLF = bestidealLF(m,BFIT[0],BFIT[1],BFIT[2],BFIT[3])
    bestrLF = bestrealLF(m,BFIT[0],BFIT[1],BFIT[2],BFIT[3]) 


    # In[ ]:

    #TRGB color information
    TRGBcolor = [0., 0., 0.,]
    Tcolor_list = np.asarray([ x[0]-x[1] for x in magnitudes_VI if (x[1]>BFIT[3] and x[1]<BFIT[3]+0.1)])
    nboots = 1000
    sample_size = len(Tcolor_list)
    me = np.zeros(nboots)
    for i in range(nboots):
        sample = Tcolor_list[np.random.randint(0, sample_size, sample_size)]
        me[i] = np.median(sample)
    TRGBcolor[0] = np.median(me)
    TRGBcolor[1] = np.percentile(me,16)
    TRGBcolor[2] = np.percentile(me,84)
    #perr = 3*np.sqrt(np.diag(pcov))

    #print fit information with errors
    color_error = 2*(TRGBcolor[2]-TRGBcolor[0])
    # print("Best fit")
    # print("TRGB: ",BFIT[3], " +/- ", perr[3])
    # print("a: ",BFIT[0], " +/- ", perr[0])
    # print("b: ",BFIT[1], " +/- ", perr[1])
    # print("c: ",BFIT[2], " +/- ", perr[2])
    # print("Color: ", TRGBcolor[0], "+/- ", color_error)


    # In[ ]:

    #plot of luminosity function, sobel filter, and best fit
    arr_roi=[]
    for x in zip(realV,realI):
        if (x[0]-x[1]>=PLT['cl'] and x[0]-x[1]<=PLT['ch'] and x[1]>=INIT['TRGB']+RANGE['l'] and x[1]<=INIT['TRGB']+RANGE['h']):
            arr_roi.append(x)
    arr_roi = np.asarray(arr_roi)
    fig=plt.figure()


    # Scaling fitting result for demonstration #
    for i in range(len(m)):
        if m[i] >= BFIT[3]:
            scale_r = (smoothLF[i]+smoothLF[i-1])/(bestrLF[i]+bestrLF[i-1])
            break

    pltfittics = []
    pltbestrLF = []
    for i in range(len(m)-1):
        if m[i] >= INIT['TRGB']+RANGE['l'] and m[i] <= INIT['TRGB']+RANGE['h']:
            pltfittics.append(m[i])
            pltbestrLF.append(bestrLF[i]*scale_r)
    #                if bestrLF[i]*scale_r<0.001 : print magtics_smooth[i],bestrLF[i]

    lf = plt.figure(figsize=(8,6))
    ax = lf.add_subplot(1,1,1)
    majorLocator = MultipleLocator(0.5)
    minorLocator = MultipleLocator(0.1)
    plt.rc('font',family='serif',size=18)
    ax.set_xlabel(band2)
    ax.set_ylabel("N")
    smoothLF_for_plot = smoothLF * normalization
    ax.plot(m,smoothLF_for_plot, 'k',color='black', label='Luminosity function')
    plt.rc('xtick',labelsize='12')
    plt.rc('ytick',labelsize='12')

    ax.set_xlim((PLT['ml'],PLT['mh']))
    ymax = np.amax(smoothLF_for_plot)
    sobel_max = np.amax(LF_sobel)
    LF_sobel_for_plot = LF_sobel/sobel_max * ymax
    ax.plot(m,LF_sobel_for_plot,':',color='black', label='Sobel filter')
    ax.set_ylim((0,ymax+(ymax/10))) 
    pltbestrLF_for_plot = [x*normalization for x in pltbestrLF]
    ax.plot(pltfittics, pltbestrLF_for_plot, 'k--',color='black', linewidth=2, label='Best fit')
    ax.plot([ INIT['TRGB']+RANGE['l'], INIT['TRGB']+RANGE['l'] ], [0,ymax], '--', color='black')
    ax.plot([ INIT['TRGB']+RANGE['h'], INIT['TRGB']+RANGE['h'] ], [0,ymax], '--', color='black')
    ax.plot([BFIT[3],BFIT[3]],[0,ymax], '-', color='black')
    ax.plot([BFIT[3]-perr[3],BFIT[3]-perr[3]],[0,ymax], '-.', color='black')
    ax.plot([BFIT[3]+perr[3],BFIT[3]+perr[3]],[0,ymax], '-.', color='black')
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator) 
    legend = ax.legend(loc='upper left',fontsize=14)
    # plt.savefig('%slf.ps' %iLoopString, dpi=300) 
 
    TipColorString    = str(TRGBcolor[0])
    ColorErrorString  = str(color_error)

    BFITSTRING0 = str(BFIT[0])
    BFITSTRING1 = str(BFIT[1]) 
    BFITSTRING2 = str(BFIT[2])
    BFITSTRING3 = str(BFIT[3])

    perrString0  = str(perr[0])
    perrString1  = str(perr[1])
    perrString2  = str(perr[2])
    perrString3  = str(perr[3])

    f.write(iLoopString)
    f.write(' TRGB:')
    f.write(BFITSTRING3)
    f.write(" +/- ")
    f.write(perrString3)
    f.write(' a: ')
    f.write(BFITSTRING0)
    f.write(" +/- ")
    f.write(perrString0)
    f.write(' b: ')
    f.write(BFITSTRING1)
    f.write(" +/- ")
    f.write(perrString1)
    f.write(' c: ')
    f.write(BFITSTRING2)
    f.write(" +/- ")
    f.write(perrString2)
    f.write(' Color: ')
    f.write(TipColorString)
    f.write(" +/- ")
    f.write(ColorErrorString)
    f.write('\n')

    plt.close('all')
    
    # print("TRGB: ",BFIT[3], " +/- ", perr[3])
    # print("a: ",BFIT[0], " +/- ", perr[0])
    # print("b: ",BFIT[1], " +/- ", perr[1])
    # print("c: ",BFIT[2], " +/- ", perr[2])
    # print("Color: ", TRGBcolor[0], "+/- ", color_error)

    iLoop = iLoop+1

f.close()


