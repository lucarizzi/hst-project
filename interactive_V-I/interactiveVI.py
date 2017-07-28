
# coding: utf-8

# In[39]:

import csv
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import gaussian_kde
from matplotlib.patches import Rectangle

# In[40]:

#prefix = input("What is the name of the galaxy?")
#for now,:
#prefix='NGC5128-S1'
#prefix='CENA-131952'
prefix='CENA-132302'

#real photometry first
real_file = "%s.phot" % prefix
real_df = pd.read_csv(real_file, delim_whitespace=True, header=None) #create data frame 
real_df = real_df.iloc[:,:37] #select till relevant end of data frame (quality flag of I)

#name columns
columns=['extension','chip','x','y','chi','snr','sharp','round','maj_ax','crowd','type',
           'cts_V','sky_v','nrm_ct_rt_V','nrm_ct_rt_err_V','inst_vega_v','mag_V','Verr','chi_V','snr_v',
           'sharp_V','round_V','crowd_V','flag_V','cts_I','sky_I','nrm_ct_rt_I','nrm_ct_rt_err_I', 
           'inst_vega_I','mag_I','Ierr','chi_I','snr_I','sharp_I','round_I','crowd_I','flag_I']



#assign columns and preview
real_df.columns=columns

#fix y column if not on first chip 
real_df['y'] = real_df.apply(lambda x: x.y+2000 if x.extension>1 else x.y, axis=1)
#real_df


# In[42]:

#define V-I
real_df['V-I'] = real_df['inst_vega_v']-real_df['inst_vega_I']
#real_df


# In[43]:

#establish cuts from below line
#$5<=2.5 && $7*$7<=0.09 && $11<=2 && $20>=5 && $24==0 && $33>=5 && $37==0'' > ${TARG}.phot2")
real_cut = real_df[(real_df['chi'] < 2.5) & (real_df['sharp'] <= 0.09) & (real_df['type'] <= 2) 
        & (real_df['snr_v'] >= 5) & (real_df['flag_V'] == 0) & (real_df['snr_I'] >= 5) 
        & (real_df['flag_I'] == 0)]



# #density plot
x = real_cut['x']
y = real_cut['y']
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=100, edgecolor='')

print('MUST SELECT FROM TOP LEFT TO BOTTOM RIGHT!!!')

# #below adapted from 
# #http://stackmirror.cn/page/qxkzht05cnjt
# should put in its own file to cut down clutter here
class Annotate(object):
  def __init__(self):
      self.ax = plt.gca()
      self.rect = Rectangle((0,0), 1, 1, facecolor='None', edgecolor='green')
      self.x0 = None
      self.y0 = None
      self.x1 = None
      self.y1 = None
      self.ax.add_patch(self.rect)
      self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
      self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
      self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
  def on_press(self, event):
      print('press')
      self.x0 = event.xdata
      self.y0 = event.ydata    
      self.x1 = event.xdata
      self.y1 = event.ydata
      self.rect.set_width(self.x1 - self.x0)
      self.rect.set_height(self.y1 - self.y0)
      self.rect.set_xy((self.x0, self.y0))
      self.rect.set_linestyle('dashed')
      self.ax.figure.canvas.draw()
  def on_motion(self,event):
      if self.on_press is True:
          return
      self.x1 = event.xdata
      self.y1 = event.ydata
      self.rect.set_width(self.x1 - self.x0)
      self.rect.set_height(self.y1 - self.y0)
      self.rect.set_xy((self.x0, self.y0))
      self.rect.set_linestyle('dashed')
      self.ax.figure.canvas.draw()
  def on_release(self, event):
      print('release')
      self.x1 = event.xdata
      self.y1 = event.ydata
      self.rect.set_width(self.x1 - self.x0)
      self.rect.set_height(self.y1 - self.y0)
      self.rect.set_xy((self.x0, self.y0))
      self.rect.set_linestyle('solid')
      self.ax.figure.canvas.draw()
      print(self.x0,self.x1,self.y0,self.y1)
      f = open('coordinates.txt', 'w')
      f.write('%d %d %d %d' % (self.x0, self.x1, self.y0, self.y1))
      f.close()
      return [self.x0,self.x1,self.y0,self.y1]
a = Annotate()
plt.show()

coordinate_file = open('coordinates.txt')
for line in coordinate_file:
    coords = line.strip().split()
    x_low = int(coords[0])
    x_high= int(coords[1])
    y_high = int(coords[2])
    y_low = int(coords[3])
coordinate_file.close()

print(x_low, x_high, y_low, y_high)

real_galaxy = real_cut[(real_cut['x'] > x_low) & (real_cut['x'] < x_high) 
        & (real_cut['y'] > y_low) & (real_cut['y'] < y_high)]

#V-I
mpl.rcParams['agg.path.chunksize'] = 10000
real_galaxy.plot(x='V-I',y='inst_vega_I', marker='o', linestyle='None', markersize='2')
axes = plt.gca()
axes.set_xlim([-2,4])
axes.set_ylim([28,19])
plt.show()






