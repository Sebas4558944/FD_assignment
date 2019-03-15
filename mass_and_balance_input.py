# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:57:56 2019

@author: Rick
"""

import numpy as np
from conversion_helpers import inch_to_m, datum_to_LEMAC
from Cit_par import c
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               Inputs
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# The first entry is the mass of the empty aircraft in pounds
# The second entry is the moment due to the empty aircraft
oew = [1.2, 5]

# The first entry is the mass of the fuel in pounds
# The second entry is the moment due to the fuel
fuelload = [[2, 5], 
            [2, 3],
            [1, 4]]


#                               Passengers
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# The first entry is the x-location from the c.g. datum in inches
# The second entry is the mass of the passenger in pounds
seat1 = [1,2]
seat2 = [2,3]
seat3 = [1,2]
seat4 = [2,1]
seat5 = [1,4]
seat6 = [2,1]
seat7 = [1,2]
seat8 = [2,1]
seat9 = [1,2]
seat10 = [2,1]


#                               Baggage
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# The first entry is the x-location from the c.g. datum in inches
# The second entry is the mass of the baggage in pounds
nose = [74, 1]
aft1 = [321, 1]
aft2 = [338, 1]


#                              Total payload
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# Collect all the data from all passengers and baggage
payload = [seat1, seat2, seat3, seat4, seat5, seat6, seat7, seat8, seat9, seat10, nose, aft1, aft2]


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               Computations
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

def cg(oew,payload,fuelload):
    # Add the moment created by the passengers and baggage to the data
    for i in range(len(payload)):
        payload[i].append(payload[i][0]*payload[i][1])
       
    zero_fuel_mass = [oew[0] + sum(payload[:][1]), oew[1] + sum(payload[:][2])] # WE NEED TO ADJUST THIS
    print zero_fuel_mass
    for i in range(len(fuelload)):
        ramp_mass = [zero_fuel_mass[0] + fuelload[i][0], zero_fuel_mass[1] + fuelload[i][1]]
        print ramp_mass
    # Obtain the x-location of the cg w.r.t. the nose of the aircraft, in inches    
        xcg_ramp_mass = ramp_mass[1] / ramp_mass[0]
        print xcg_ramp_mass
    # Obtain the x-location of the cg w.r.t. the LEMAC of the aircraft, in meters 
        cgLEMAC_ramp_mass = xcg_ramp_mass*inch_to_m + datum_to_LEMAC
        
    # Convert the x-location of the cg w.r.t. the LEMAC in percentage of MAC
        xcg_percentage = cgLEMAC_ramp_mass / c * 100.
        print xcg_percentage
    return xcg_percentage

