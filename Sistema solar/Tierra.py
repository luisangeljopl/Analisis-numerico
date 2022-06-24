import numpy as np
import matplotlib.pyplot as plt
#Set parameters:
N = 365                     # Earth days in a year
dt = 1.00/N                 # Time Step: Fractions of a year - 1 Earth day (i.e. 1/365)
mu = 4 * np.pi**2           # mu=4pi^2 is the Gravitational Parameter: mu = GM where G=6.67e-11 is the Universal Gravitational Constant and M is the mass of the body
rEar = 1

#Create an array, for all variables, of size N with all entries equal to zero:
xEar = np.zeros((N,))
yEar = np.zeros((N,))
vxEar = np.zeros((N,))
vyEar = np.zeros((N,))

# Initial Conditions:
xEar[0] = rEar                   # (x0 = r, y0 = 0) in AU
vyEar[0] = np.sqrt(mu/rEar)      # (vx0 = 0, v_y0 = sqrt(mu/r)) AU/yr

#Implement Verlet Algorithm:
for k in range(0,N-1):
    vxEar[k+1] = vxEar[k] - (mu*xEar[k]) / (rEar**3)*dt
    xEar [k+1] = xEar[k] + vxEar[k+1]*dt
    vyEar[k+1] = vyEar[k] - (mu*yEar[k]) / (rEar**3)*dt
    yEar [k+1] = yEar[k] + vyEar[k+1]*dt

#Plot:
plt.plot(xEar, yEar, 'go')
plt.title ('Circular Orbit of Earth')
plt.xlabel ('x')
plt.ylabel ('y')
plt.axis('equal')
plt.show()

