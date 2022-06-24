import numpy as np
import matplotlib.pylab as plt
import math
from pylab import *

""" example using eulers method for solving the ODE
y’(x) = f(x, y) = y
y(0) = 1
Eulers method:
y^(n + 1) = y^(n) + h*f(x, y^(n)), h = dx
"""

N = 19
x = np.linspace(0,2.5, N)
h = 0.1 # steplength
y_0 = 0 # initial condition
Y = np.zeros_like(x) # vector for storing y values
q= arange(19.)
p=y=(3/125)*q**2+(1/5)*q
Y[0] = y_0 # first element of y = y(0)
for n in range(N):
    f = -9.8+(2000/(200-x))+(2/(x-200))*Y
f[n] =  f[n] + h*f[n]
f_analytic = np.exp(x)-1
print(x)
print(f_analytic)


# construcci´on de gr´afica
plt.title("VELOCIDAD DE UN COHETE")
plt.plot( x ,f_analytic , linewidth=1.0)
plt.plot(x,p,color="red")
plt.xlabel("t")
plt.ylabel("v(t)")
plt.plot(x, f_analytic , '-r', label='solución númerica',color="blue")
plt.plot(x, p, '-r', label='solución analítica')
plt.legend(loc='upper left')
plt.ylim(0, 15)
#plt.savefig(’../fig-ch1/euler_simple.pdf’, transparent=True)       
plt.show()


