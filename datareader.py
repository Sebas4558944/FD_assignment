# -*- coding: utf-8 -*-
"""
Created on Wed Mar 06 09:10:19 2019

@author: daanv
"""

import numpy as np
import mat4py
import subprocess

def convertToSec(time):
    splitTime=time.split(':')
    try:
        if splitTime[2]=='00':
            if splitTime[0]=='01':
                h=int(splitTime[0])
                m=int(splitTime[1])
                s=int(splitTime[2])
            else:
                h=0
                m=int(splitTime[0])
                s=int(splitTime[1])
        else:
            h=int(splitTime[0])
            m=int(splitTime[1])
            s=int(splitTime[2])
    except IndexError:
            h=0
            m=int(splitTime[0])
            s=int(splitTime[1])
    t=3600*h+60*m+s
    return t

def convertToTimeStr(h,m,s): #hh:mm:ss format
#    s=s-s%0.1 #Get rid of anything smaller than 0.1

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

    #ensure no decimals
    s=int(s)
    m=int(m)
    h=int(h)

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
        timeStr=mstr+':'+sstr+':00'
    else:
        timeStr=hstr+':'+mstr+':'+sstr

    return timeStr


def importExcelData(f):
    # f is filename
    # Process defines wheter or not data is completed and converted to SI \
    # or standard units, default is false.

    arr_str=np.genfromtxt(f,delimiter=',',dtype='str')

    arr=np.genfromtxt(f,delimiter=',',dtype='float')

    #y,x
    date_of_flight=arr_str[2,3]
    flight_number=arr_str[3,3]
    TO_time=arr_str[2,5]
    LND_time=arr_str[3,5]
    #order = [p1,p2,coord,1L,1R,2L,2R,3L,3R]
    passengerMass=arr[7:16,7] #in kg
    passengerNames=arr_str[7:16,3]
    passengerPos=arr_str[7:16,0]


    blockfuel=arr[17,3] #in lbs
    #Aircraft config
    ACC_CLCD=arr_str[22,4]

    #time [min:sec], Elapsed time [sec] #empty, hp [ft] (pressure altitude)\
    #IAS [kts], a [deg], FFL [lbs/hr], FFr [lbs/hr], F. used [lbs], TAT #C
    CL_CD_series1=arr[27:33,1:9]
    CL_CD_series1str=arr_str[27:33,1:9]
    CL_CD_series2=arr[43:49,1:9]
    CL_CD_series2str=arr_str[43:49,1:9]
    
    
    #fill in ET
    for i in range(len(CL_CD_series1[:,0])):
        time=convertToSec(CL_CD_series1str[i,0])
        CL_CD_series1[i,1]=int(time)
        try:
            time=convertToSec(CL_CD_series2str[i,0])
            CL_CD_series2[i,1]=int(time)
        except ValueError:
            pass
            
    #Aircraft config
    ACC_Trim=arr_str[53,4]

    #time [hrs:min], Elapsed time [sec] #empty, hp [ft] (pressure altitude)\
    #IAS [kts], a [deg], de [deg], detr [deg], Fe [N], FFL [lbs/hr], \
    #FFr [lbs/hr], F. used [lbs], TAT #C
    El_Trim_Curve=arr[58:65,1:12]
    El_Trim_Curvestr=arr_str[58:65,1:12]
    
    #fill in ET
    for i in range(len(El_Trim_Curve[:,0])):
        time=convertToSec(El_Trim_Curvestr[i,0])
        El_Trim_Curve[i,1]=int(time)
        
    #CG shift: 
    #_shifted all relates to moved person
    name_shifted=arr_str[69,1]
    pos_shifted=arr_str[69,1]
    newpos_shifted=arr_str[70,4]

    #time [hrs:min], Elapsed time [sec] #empty, hp [ft] (pressure altitude)\
    #IAS [kts], a [deg], de [deg], detr [deg], Fe [N], FFL [lbs/hr], \
    #FFr [lbs/hr], F. used [lbs], TAT #C
    Cg_shift=arr[74:76,1:13]
    Cg_shiftstr=arr_str[74:76,1:13]
    #fill in ET
    for i in range(len(Cg_shift[:,0])):
        time=convertToSec(Cg_shiftstr[i,0])
        Cg_shift[i,1]=int(time)
        

    #Eigenmotions
    eigenmotions=[]
    phugoid=arr_str[82,3]
    shortPeriod=arr_str[83,3]
    dutchRoll=arr_str[82,6]
    dutchRollYD=arr_str[83,6]
    aperRoll=arr_str[82,9]
    spiral=arr_str[83,9]
    eigenmotions.extend((phugoid,shortPeriod,dutchRoll,dutchRollYD,aperRoll,spiral))

    return date_of_flight, flight_number, TO_time, LND_time, passengerMass, \
           passengerNames, passengerPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2, \
           ACC_Trim, El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, \
           Cg_shift, eigenmotions

def getFDValues(f):
    data = mat4py.loadmat(f)
    flightdata = data.get('flightdata', {})
    keylist = []
    desclist = []
    #    datlist=[]
    unitlist = []
    newdict = {}
    for key in flightdata:
        keydict = flightdata.get(key, {})
        keydesc = keydict.get('description', {})
        keydat = keydict.get('data', {})
        keyunits = keydict.get('units', {})
        keylist.append(key)
        desclist.append(keydesc)
        #        datlist.append(keydat) #Costs comp time.
        unitlist.append(keyunits)
        newdict[key] = keydat
    return keylist, desclist, unitlist, newdict


def ThrustingAllDayEveryday(write):
    #write - list of input variables to be written to matlab.dat
    with open('matlab.dat', 'w') as f:
        line=''
        for i in write:
            line+=str(i)+' '
        f.write(line)
        f.close
    #Run exe
    subprocess.call('thrust_new.exe')
    #Read exe output
    out=np.genfromtxt('thrust.dat')
    return out

#%%
#    Test functions
#    
#f='Reference_Datasheet.csv'
#date_of_flight, flight_number, TO_time, LND_time, passengerMass, passengerNames\
#, passengerPos, blockfuel, ACC_CLCD, CL_CD_series1, CL_CD_series2, ACC_Trim,\
# El_Trim_Curve, name_shifted, pos_shifted, newpos_shifted, Cg_shift, eigenmotions \
# = importExcelData(f)


# print eigenmotions[0]
# print convertToSec(eigenmotions[0])
# print convertToTimeStr(0,0,convertToSec(eigenmotions[0]))

# keylist,desclist,unitlist,newDict=getFDValues('reference.mat')
#print convertToSec("01:00:00")