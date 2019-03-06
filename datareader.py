# -*- coding: utf-8 -*-
"""
Created on Wed Mar 06 09:10:19 2019

@author: daanv
"""

import numpy as np

def importRawData(f):
    arr=np.genfromtxt(f,delimiter=',',dtype='str')
    #y,x
    date_of_flight=arr[2,3] 
    flight_number=arr[3,3] 
    TO_time=arr[2,5]
    LND_time=arr[3,5]
    #weights = [p1,p2,coord,1L,1R,2L,2R,3L,3R]
    weights=arr[7:16,7]
    weightNames=arr[7:16,3]
    weightPos=arr[7:16,0]
    
    blockfuel=arr[17,3]
    
    #Aircraft config
    ACC_CLCD=arr[22,4]
    
    #time [hrs:min], Elapsed time [sec] #empty, hp [ft] (pressure altitude)\
    #IAS [kts], a [deg], FFL [lbs/hr], FFr [lbs/hr], F. used [lbs], TAT #C
    CL_CD_series1=arr[27:33,1:9]
    CL_CD_series2=arr[40:49,1:9]
    
    #Aircraft config
    ACC_Trim=arr[53,4]
    
    #time [hrs:min], Elapsed time [sec] #empty, hp [ft] (pressure altitude)\
    #IAS [kts], a [deg], de [deg], detr [deg], Fe [N], FFL [lbs/hr], \ 
    #FFr [lbs/hr], F. used [lbs], TAT #C
    El_Trim_Curve=arr[58:64,1:12]
    
    #CG shift: 
    #_shifted all relates to moved person
    name_shifted=arr[69,1]
    pos_shifted=arr[69,1]
    newpos_shifted=[70,4]
    
    #time [hrs:min], Elapsed time [sec] #empty, hp [ft] (pressure altitude)\
    #IAS [kts], a [deg], de [deg], detr [deg], Fe [N], FFL [lbs/hr], \ 
    #FFr [lbs/hr], F. used [lbs], TAT #C
    Cg_shift=arr[74:75,1:12]
    
    #Eigenmotions
    eigenmotions=[]
    phugoid=arr[82,3]
    shortPeriod=arr[83,3]
    dutchRoll=arr[82,6]
    dutchRollYD=arr[83,6]
    aperRoll=arr[82,9]
    spiral=arr[83,9]
    eigenmotions.extend((phugoid,shortPeriod,dutchRoll,dutchRollYD,aperRoll,spiral))

    return date_of_flight, flight_number, TO_time, LND_time, weights, \
    weightNames, weightPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2,\
    ACC_Trim, El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, \
    Cg_shift, eigenmotions

f='Reference_Datasheet.csv'
listed=importRawData(f)