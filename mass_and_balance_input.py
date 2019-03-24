# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:57:56 2019

@author: Rick
"""

import numpy as np
from conversion_helpers import inch_to_m, datum_to_LEMAC, lbs_to_kg
from Cit_par import c
from datareader import *

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                               Inputs
# -----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# The first entry is the mass of the empty aircraft in lbs
# The second entry is the moment due to the empty aircraft [lbs-inch]
oew = [9165, 2677847.5]
mac_length = 80.98  # inch
dist_to_lemac = 261.56  # inch

#                               Fuel
# -----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# Import the mass of the fuel burnt since take-off in lbs
fuel = importExcelData('Post_Flight_Datasheet_13_03_V2.csv')[12][:, -2]
fuel_shifted = importExcelData('Post_Flight_Datasheet_13_03_V2.csv')[16][:, 10]

# Mass of the fuel at start of flight, in lbs
start_fuel = 2800.

# Get a list of the mass of the fuel at start of engines minus the mass of the burnt fuel, in kg
init_fuel = np.full((1, len(fuel)), start_fuel)
fuel_load = np.array((init_fuel - fuel) * lbs_to_kg)
p = np.transpose(fuel_load)

init_fuel_shifted = np.full((1, len(fuel_shifted)), start_fuel)
fuel_load_shifted = np.array((init_fuel_shifted - fuel_shifted) * lbs_to_kg)
q = np.transpose(fuel_load_shifted)

# Moment arms of the fuel given in Table E2, in inch-lbs
mfuel = np.array([[5437.01], [5341.72], [5272.02], [5158.23]])
mfuel_shifted = np.array([[4990.32], [4907.79]])

# Convert to Nm
fuelmoment = np.array(mfuel * inch_to_m * lbs_to_kg)
fuelmoment_shifted = np.array(mfuel_shifted * inch_to_m * lbs_to_kg)

# The first entry is the mass of the fuel in kg
# The second entry is the moment due to the fuel in Nm
fuelload = np.hstack((p, fuelmoment))
fuelload_shifted = np.hstack((q, fuelmoment_shifted))

#                               Passengers
# -----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# The first entry is the x-location from the c.g. datum in inches
# The second entry is the mass of the passenger in kg
seat1 = [131., 89.]  # Alexander
seat2 = [131., 82.]  # Hans
seat3 = [214., 62.]  # Carmen
seat4 = [214., 70.]  # Rick
seat5 = [251., 74.]  # Daan
seat6 = [251., 65.]  # Sera
seat7 = [288., 80.]  # Martijn
seat8 = [288., 82.]  # Sebastiaan
seat9 = [0., 0]
seat10 = [170., 80.]  # Tigran

#                              Total payload
# -----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# Collect all the data from all passengers and baggage
pay = [seat1, seat2, seat3, seat4, seat5, seat6, seat7, seat8, seat9, seat10]
payload = np.array(pay)

# Prepare for seat shift
seat8_shifted = [134., 82.]
pay_shifted = [seat1, seat2, seat3, seat4, seat5, seat6, seat7, seat8_shifted, seat9, seat10]
payload_shifted = np.array(pay_shifted)

# Convert locations from inch to m
payload[:, 0] = payload[:, 0] * inch_to_m
payload_shifted[:, 0] = payload_shifted[:, 0] * inch_to_m

# Prepare list to include moments by passengers and baggage
payl = np.zeros((len(pay), 1))
payl_shifted = np.zeros((len(pay_shifted), 1))


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                           Computation for regular postions
# -----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def cg():
    # Add the moment created by the passengers and baggage to the data
    for i in range(len(payload)):
        payl[i] = payload[i, 0] * payload[i, 1]

    Payload = np.hstack((payload, payl))
    zero_fuel_mass = [oew[0]*lbs_to_kg + sum(Payload[:, 1]), oew[1]*lbs_to_kg*inch_to_m + sum(Payload[:, 2])]

    total_mass = []
    cg_location = []

    for i in range(len(fuelload)):
        ramp_mass = [zero_fuel_mass[0] + fuelload[i, 0], zero_fuel_mass[1] + fuelload[i, 1]]
        total_mass.append(ramp_mass[0])
        # Obtain the x-location of the cg w.r.t. the nose of the aircraft, in meters
        xcg_ramp_mass = ramp_mass[1] / ramp_mass[0]

        # Obtain the x-location of the cg w.r.t. the LEMAC of the aircraft, in meters
        cgLEMAC_ramp_mass = xcg_ramp_mass + datum_to_LEMAC

        # Convert the x-location of the cg w.r.t. the LEMAC in percentage of MAC
        xcg_percentage = cgLEMAC_ramp_mass / c * 100.
        cg_location.append(xcg_percentage)
        print "cg location from nose:", xcg_ramp_mass
    print "total mass:", total_mass
        
    return cg_location, total_mass


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                           Computation for shifted Sebastiaan
# -----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def cg_shifted():
    # Add the moment created by the passengers and baggage to the data
    for i in range(len(payload_shifted)):
        payl_shifted[i] = payload_shifted[i, 0] * payload_shifted[i, 1]

    Payload_shifted = np.hstack((payload_shifted, payl_shifted))

    zero_fuel_mass_shifted = [oew[0]*lbs_to_kg + sum(Payload_shifted[:, 1]), oew[1]*lbs_to_kg*inch_to_m + sum(Payload_shifted[:, 2])]

    total_mass_shifted = []
    cg_location = []

    for i in range(len(fuelload_shifted)):
        ramp_mass_shifted = [zero_fuel_mass_shifted[0] + fuelload_shifted[i, 0], zero_fuel_mass_shifted[1] + fuelload_shifted[i, 1]]
        total_mass_shifted.append(ramp_mass_shifted[0])

        # Obtain the x-location of the cg w.r.t. the nose of the aircraft, in meters
        xcg_ramp_mass_shifted = ramp_mass_shifted[1] / ramp_mass_shifted[0]

        # Obtain the x-location of the cg w.r.t. the LEMAC of the aircraft, in meters
        cgLEMAC_ramp_mass_shifted = xcg_ramp_mass_shifted + datum_to_LEMAC

        # Convert the x-location of the cg w.r.t. the LEMAC in percentage of MAC
        xcg_percentage_shifted = cgLEMAC_ramp_mass_shifted / c * 100.
        cg_location.append(xcg_percentage_shifted)

    return cg_location, total_mass_shifted

q = cg()
print q
