# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 16:39:50 2019

@author: daanv
"""


from datareader import getFDValues,importExcelData,convertToSec,convertToTimeStr
import matplotlib.pyplot as plt

keylist,desclist,unitlist,newDict=getFDValues('reference.mat')

#%%

x=newDict.get('time')
y=newDict.get('delta_e')
b=32281
e=34279
plt.plot(x[b:e],y[b:e])
plt.show()