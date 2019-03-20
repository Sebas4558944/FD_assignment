# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:57:56 2019

@author: Rick
"""

import numpy as np
from conversion_helpers import inch_to_m, datum_to_LEMAC, lbs_to_kg
from Cit_par import c
import datareader
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                               Inputs
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# The first entry is the mass of the empty aircraft in pounds
# The second entry is the moment due to the empty aircraft
oew = [9165, 2677847.5]


#                               Fuel
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# Import the mass of the fuel burnt since take-off
fuel = importExcelData('Post_Flight_Datasheet_13_03_V2.csv')[12]
fuel = fuel[:,-1]

fuel_shifted = importExcelData('Post_Flight_Datasheet_13_03_V2.csv')[16]
fuel_shifted = fuel_shifted[:,10]

# Mass of the fuel at start of engines, in lbs
startfuel = 2800.

# Get a list of the mass of the fuel at start of engines minus the mass of the burnt fuel, in lbs
initfuel = np.full((1,len(fuel)),startfuel)
fuelload = np.array((initfuel - fuel)*lbs_to_kg)
p = np.transpose(fuelload)

initfuel_shifted = np.full((1,len(fuel_shifted)),startfuel)
fuelload_shifted = np.array((initfuel_shifted - fuel_shifted)*lbs_to_kg)
q = np.transpose(fuelload_shifted)

# Moment arms of the fuel given in Table E2, in inch-lbs
mfuel = np.array([[10000],[9000],[8000],[7000],[6000],[5000],[4000]])
mfuel_shifted = np.array([[7500],[8500]])

# Convert to Nm
fuelmoment = np.array(mfuel*inch_to_m*lbs_to_kg)
fuelmoment_shifted = np.array(mfuel_shifted*inch_to_m*lbs_to_kg)

# The first entry is the mass of the fuel in pounds
# The second entry is the moment due to the fuel
fuelload = np.hstack((p,fuelmoment))
fuelload_shifted = np.hstack((q,fuelmoment_shifted))

#                               Passengers
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# The first entry is the x-location from the c.g. datum in inches
# The second entry is the mass of the passenger in pounds
seat1 = [131.,89.]               # Alexander
seat2 = [131.,82.]               # Hans
seat3 = [214.,62.]               # Carmen
seat4 = [214.,70.]               # Rick
seat5 = [251.,74.]               # Daan
seat6 = [251.,65.]               # Sera
seat7 = [288.,80.]               # Martijn
seat8 = [288.,82.]               # Sebastiaan
seat9 = [0.,0]
seat10 = [170.,80.]              # Tigran



#                              Total payload
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# Collect all the data from all passengers and baggage
pay = [seat1, seat2, seat3, seat4, seat5, seat6, seat7, seat8, seat9, seat10]
payload = np.array(pay)

# Prepare for seat shift
seat8_shifted = [134.,82.]
pay_shifted = [seat1, seat2, seat3, seat4, seat5, seat6, seat7, seat8_shifted, seat9, seat10]
payload_shifted = np.array(pay_shifted)

payload[:,0] = payload[:,0]*inch_to_m
payload_shifted[:,0] = payload_shifted[:,0]*inch_to_m

# Prepare list to include moments by passengers and baggage
payl = np.zeros((len(pay),1))
payl_shifted = np.zeros((len(pay_shifted),1))

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                           Computation for regular postions
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

def cg(oew,payload,fuelload):
    # Add the moment created by the passengers and baggage to the data
    for i in range(len(payload)):
        payl[i] = payload[i,0]*payload[i,1]

    Payload = np.hstack((payload,payl))

    zero_fuel_mass = [oew[0] + sum(Payload[:,1]), oew[1] + sum(Payload[:,2])] 
    
    total_mass = []
    
    for i in range(len(fuelload)):
        ramp_mass = [zero_fuel_mass[0] + fuelload[i,0], zero_fuel_mass[1] + fuelload[i,1]]
        total_mass.append(ramp_mass[0])        
    # Obtain the x-location of the cg w.r.t. the nose of the aircraft, in inches    
        xcg_ramp_mass = ramp_mass[1] / ramp_mass[0]
   
    # Obtain the x-location of the cg w.r.t. the LEMAC of the aircraft, in meters 
        cgLEMAC_ramp_mass = xcg_ramp_mass + datum_to_LEMAC

    # Convert the x-location of the cg w.r.t. the LEMAC in percentage of MAC
        xcg_percentage = cgLEMAC_ramp_mass / c * 100.
   
    return xcg_percentage, total_mass

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                           Computation for shifted Sebastiaan
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------

def cg_shifted(oew,payload_shifted,fuelload_shifted):
    # Add the moment created by the passengers and baggage to the data
    for i in range(len(payload_shifted)):
        payl_shifted[i] = payload_shifted[i,0]*payload_shifted[i,1]

    Payload_shifted = np.hstack((payload_shifted,payl_shifted))

    zero_fuel_mass = [oew[0] + sum(Payload_shifted[:,1]), oew[1] + sum(Payload_shifted[:,2])] 
    
    total_mass_shifted = []
    
    for i in range(len(fuelload_shifted)):
        ramp_mass_shifted = [zero_fuel_mass[0] + fuelload_shifted[i,0], zero_fuel_mass[1] + fuelload_shifted[i,1]]
        total_mass_shifted.append(ramp_mass_shifted[0]) 
        
    # Obtain the x-location of the cg w.r.t. the nose of the aircraft, in inches    
        xcg_ramp_mass_shifted = ramp_mass_shifted[1] / ramp_mass_shifted[0]
   
    # Obtain the x-location of the cg w.r.t. the LEMAC of the aircraft, in meters 
        cgLEMAC_ramp_mass_shifted = xcg_ramp_mass_shifted + datum_to_LEMAC

    # Convert the x-location of the cg w.r.t. the LEMAC in percentage of MAC
        xcg_percentage_shifted = cgLEMAC_ramp_mass_shifted / c * 100.
   
    return xcg_percentage_shifted, total_mass_shifted
