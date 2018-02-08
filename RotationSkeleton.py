'''
# RotationSkeleton.py
# M. Grant & J. Nevin MPhys 2018
# Monte Carlo simulation of vortex ring transport in rotational system.
'''

import matplotlib.pyplot as plt
import numpy as np
import random as rnd


# FUNCTIONS

def RadiusChange(E):
    # Calculate radius in electric field
    R = (2*e*E)/(rho*(kappa**2)*(Lam - 1.5))
    return R

def RingVelocity(R):
    # Calculate ring velocity
    v = (kappa/(4*np.pi*R))*(Lam - 0.5)
    return v

#  TO DO
# Functions: Lambda, RNG, Euler Method, charge escape


# MAIN

# Constants

rho = 145 # Superfluid helium-4 density
kappa = 9.98e-8 # Quantum of circulation
e = 1.602e-19 # Electron charge

# Parameters

E = 1 # Electric Field
X = 0.045 # Box Size
a = 1e-10 # Vortex Core Size
R_0 = 0.1e-6 # Initial ring radius
dt = 0.0001 # Time step

# Empty arrays for stored values


# While loop (boundary conditions)

# Calculate P of interaction and RNG, compare to determine event

# IF statement: collision
    # RNG impact param ---> determine final Radius using fit
    # Check if charge escaped
    # if charge on ring take step and repeat
    # Reset interaction time


# ELSE : No collision
    # Calulate radius due to E Field
    # add step to time and interaction time

# Store time of absorption (time of flight)
=======

 
