# -*- coding: utf-8 -*-
"""
Created on Wed Mar 06 09:10:19 2019

@author: daanv
"""

import numpy as np
import mat4py
from conversion_helpers import kts_to_ms,ft_to_m,lbs_to_kg
f2m=ft_to_m
l2k=lbs_to_kg
k2m=kts_to_ms

def importExcelData(f):
    #f is filename
    #Process defines wheter or not data is completed and converted to SI \
    # or standard units, default is false.
    
    arr=np.genfromtxt(f,delimiter=',',dtype='str')

    #y,x
    date_of_flight=arr[2,3] 
    flight_number=arr[3,3] 
    TO_time=arr[2,5]
    LND_time=arr[3,5]
    #order = [p1,p2,coord,1L,1R,2L,2R,3L,3R]
    passengerMass=arr[7:16,7] #in kg
    passengerNames=arr[7:16,3]
    passengerPos=arr[7:16,0]
    
    
    blockfuel=arr[17,3] #in lbs
    #Aircraft config
    ACC_CLCD=arr[22,4]
    
    #time [min:sec], Elapsed time [sec] #empty, hp [ft] (pressure altitude)\
    #IAS [kts], a [deg], FFL [lbs/hr], FFr [lbs/hr], F. used [lbs], TAT #C
    CL_CD_series1=arr[27:33,1:9]
    CL_CD_series2=arr[43:49,1:9]
    
    #fill in ET
    for l in (CL_CD_series1,CL_CD_series2):
        for i in range(len(l[:,0])):
            time=CL_CD_series1[i,0]
#            print time
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

    return date_of_flight, flight_number, TO_time, LND_time, passengerMass, \
    passengerNames, passengerPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2,\
    ACC_Trim, El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, \
    Cg_shift, eigenmotions

 
def convertToSec(time):
    splitTime=time.split(':')
    if splitTime[2]=='00':
        h=0
        m=int(splitTime[0])
        s=int(splitTime[1])
    else:
        h=int(splitTime[0])
        m=int(splitTime[1])
        s=int(splitTime[2])        
    t=3600*h+60*m+s
    return t

def convertToTimeStr(h,m,s):
    if int(s)>60:
        unit=s%60
        mchange=int(m)+(int(s)-unit)/(60)
        s=int(s)-int(mchange)*60
        m=int(m)+mchange

    if int(m)>60:
        unit=m%60
        hchange=(int(m)-unit)/(60)
        m=int(m)-int(hchange)*60
        h=int(h)+hchange
    
    hstr=str(h)
    mstr=str(m)
    sstr=str(s)
    if len(hstr)==1:
        hstr='0'+hstr
    if len(mstr)==1:
        mstr='0'+mstr
    if len(sstr)==1:
        sstr='0'+sstr
    
    if int(h)==0:
        timeStr=mstr+':'+sstr+':'+hstr
    else:
        timeStr=hstr+':'+mstr+':'+sstr

    return timeStr


def importFlightData(f):
    data=mat4py.loadmat(f)
    flightdata = data.get('flightdata',{})
    return flightdata

def getFDValues(f):
    data=mat4py.loadmat(f)
    flightdata = data.get('flightdata',{})
    keylist=[]
    desclist=[]
#    datlist=[]
    unitlist=[]
    newdict={}
    for key in flightdata:
        keydict=flightdata.get(key,{})
        keydesc=keydict.get('description',{})
        keydat=keydict.get('data',{})
        keyunits=keydict.get('units',{})
        keylist.append(key)
        desclist.append(keydesc)
#        datlist.append(keydat)
        unitlist.append(keyunits)
        newdict[key] = keydat
    return keylist,desclist,unitlist,newdict



    
f='Reference_Datasheet.csv'
date_of_flight, flight_number, TO_time, LND_time, passengerMass, passengerNames\
, passengerPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2, ACC_Trim,\
 El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, Cg_shift, eigenmotions\
 = importExcelData(f)
 
#print eigenmotions[0]
#print convertToSec(eigenmotions[0])
#print convertToTimeStr(0,0,convertToSec(eigenmotions[0]))

#flightData=importFlightData('reference.mat')
#keydict=flightData.get('Gps_long',{})
#keydata=keydict.get('data',{})

keylist,desclist,unitlist,newDict=getFDValues('reference.mat')