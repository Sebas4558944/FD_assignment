# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 16:39:50 2019

@author: daanv
"""


from datareader import getFDValues,importExcelData,convertToSec,convertToTimeStr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

keylist,desclist,unitlist,newDict=getFDValues('reference.mat')

time=newDict.get('time')
lat=newDict.get('Gps_lat')
lon=newDict.get('Gps_long')
latsum=0
lonsum=0
length=0
latnew=[]
lonnew=[]
for i in range(len(lat)):
    latsum+=lat[i][0]
    lonsum+=lon[i][0]
    length+=1
    latnew.append(lat[i][0])
    lonnew.append(lon[i][0])
latavg=latsum/length
lonavg=lonsum/length
#
#plt.figure(figsize=(8, 8))
#m = Basemap(projection='lcc', resolution=None, lat_0=latavg, lon_0=lonavg,\
#            width=5E5,height=5E5)
#m.etopo(scale=0.5);
#%%
hsize=5
vsize=5
m = Basemap(projection='mill',llcrnrlat=50.5,urcrnrlat=54,\
            llcrnrlon=3.5,urcrnrlon=7,resolution='h')
m.drawcoastlines()
m.fillcontinents(color='green',lake_color='blue')
         
m.drawcountries(color='black')
m.drawrivers(color='blue')

m.drawmapboundary(fill_color='#FFFFFF')
x,y=m(lonnew,latnew)                  
m.plot(x,y,color='red')
#plt.show()