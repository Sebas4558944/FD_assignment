# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:12:34 2019

@author: daanv
"""

from datareader import getFDValues,importExcelData,convertToSec,convertToTimeStr

keylist,desclist,unitlist,newDict = getFDValues('reference.mat')

date_of_flight, flight_number, TO_time, LND_time, passengerMass, passengerNames\
, passengerPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2, ACC_Trim,\
 El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, Cg_shift, eigenmotions\
 = importExcelData('Reference_Datasheet.csv')

# eigenmotions - phugoid, short period, dutch roll,\
# dutch roll Yd, aperiodic Roll, Spiral
 
 