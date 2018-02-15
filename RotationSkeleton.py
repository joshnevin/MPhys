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

def FinalState(b):
    # Final state from impact parameter
    if b > 0.95:
        # Miss
        return -1
    else:
        # Reconnection
        R_e = (0.7731*(b**9) + 0.5533*(b**8) - 1.523*(b**7) - 1.245*(b**6) +
            0.7401*(b**5) + 0.5716*(b**4) - 0.3095*(b**3) - 0.2636*(b**2) -
                0.2434*(b) + 0.821)

        return R_e

#  TO DO
# Functions: Lambda, RNG, Euler Method, charge escape

# MAIN

# Constants

rho = 145 # Superfluid helium-4 density
kappa = 9.98e-8 # Quantum of circulation
e = 1.602e-19 # Electron charge

# Parameters

E = 2000 # Electric Field
X = 0.045 # Box Size
a = 1e-10 # Vortex Core Size
R_0 = 0.1e-6 # Initial ring radius
dt = 0.01 # Time step
omega = 1 # rad per second, rotational velocity
L = 2*omega/kappa
sigma = 2 # Interaction size in units of radius

t_f = []
R_max = []

# Loop over n rings and record current
for q in range(100):
    print(q)

    # Arrays

    R = []
    t = []
    x = []

    R.append(R_0)
    t.append(0)
    x.append(0)

    # Counter

    n = 0
    z = 0

    # Interaction timer
    t_r = 0

    # Simulation

    while x[n] < X:

        Lam = np.log(8*R[n]/a)

        vel = RingVelocity(R[n])

        # Calculate P of interaction and RNG, compare to determine event
        # Using 2R as the interaction size
        P_int = sigma*R[n]*L*vel*dt
        P_test = np.random.uniform(low=0.000, high=1.0)

        print(P_int)
        print(P_test)


        # Did an interaction happen?
        if P_test < P_int:
            print('Interaction')
            # IF statement: collision
            # RNG impact param ---> determine final Radius using fit
            # Check if charge escaped
            # if charge on ring take step and repeat
            # Reset interaction time
            # Collison happened so determine final state
            b = np.random.uniform(low=-1.0, high=1.0)
            # Determine final state
            R_e = FinalState(b)
            if R_e != -1:
                R[n] = R_e*R[n]
                # Check if charge escaped
                P_escape = R_e # ratio of emission to incident ring radii
                P_test2 = np.random.uniform(low=0.0, high=1.0)
                if P_test2 < P_escape:
                    print('Successful step')
                    # take step
                    v = RingVelocity(R[n])
                    x.append(x[n] + dt*v)
                    t.append(t[n] + dt)
                    t_r = 0 # reset interaction timer
                    dx = dt*v
                    g = RadiusChange(E)
                    R.append(R[n] + dx*g)
                else:
                    continue
                    print('Lost')


            else:
                # Miss so take normal steps
                #print('Onward')
                v = RingVelocity(R[n])
                x.append(x[n] + dt*v)
                t.append(t[n] + dt)
                t_r = 0 # reset interaction timer
                dx = dt*v
                g = RadiusChange(E)
                R.append(R[n] + dx*g)


        else:
            #print('Onward')
            # ELSE : No collision
            # Calulate radius due to E Field
            # add step to time and interaction time
            # No interaction so take normal step
            v = RingVelocity(R[n])
            x.append(x[n] + dt*v)
            t.append(t[n] + dt)
            t_r =+ dt
            dx = dt*v
            g = RadiusChange(E)
            R.append(R[n] + dx*g)


        n = n + 1

    # End

    # Store time of absorption (time of flight)
    t_f.append(t[n])
    R_max.append(R[n])

# End


fig, ax = plt.subplots()
ax.plot(t_f, '.', label='Flight Time')
#plt.ylabel('Radius (m)')
#plt.xlabel('Distance (m)')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#legend = ax.legend(loc='lower right', shadow=False, fontsize='x-large')
#plt.savefig('RvD.png', format='png', dpi=900)
plt.show()


# Calculate current by a Histogram
fig, ax = plt.subplots()
plt.hist(t_f)
plt.show()




    # ISSUES:
    # Still records t(n) etc even if ring lost, just stops earlier
