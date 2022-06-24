f = lambda x : (x**2)
a = 1; b = 3; N = 50
n = 10 # Use n*N+1 points to plot the function smoothly
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(a,b,N+1)
y = f(x)

X = np.linspace(a,b,n*N+1)
Y = f(X)

plt.figure(figsize=(15,5))

plt.subplot(1,3,1)
plt.plot(X,Y,'b')
x_left = x[:-1] # Left endpoints
y_left = y[:-1]
plt.plot(x_left,y_left,'b.',markersize=10)
plt.bar(x_left,y_left,width=(b-a)/N,alpha=0.2,align='edge',edgecolor='b')
plt.title('Suma de Riemann por la izquierda, N = {}'.format(N))
plt.xlabel('x')
plt.ylabel('f(x)')

plt.subplot(1,3,2)
plt.plot(X,Y,'b')
x_mid = (x[:-1] + x[1:])/2 # Midpoints
y_mid = f(x_mid)
plt.plot(x_mid,y_mid,'b.',markersize=10)
plt.bar(x_mid,y_mid,width=(b-a)/N,alpha=0.2,edgecolor='b')
plt.title('Suma de Riemann punto medio, N = {}'.format(N))
plt.xlabel('x')
plt.ylabel('f(x)')

plt.subplot(1,3,3)
plt.plot(X,Y,'b')
x_right = x[1:] # Left endpoints
y_right = y[1:]
plt.plot(x_right,y_right,'b.',markersize=10)
plt.bar(x_right,y_right,width=-(b-a)/N,alpha=0.2,align='edge',edgecolor='b')
plt.title('Suma de Riemann por la derecha, N = {}'.format(N))
plt.xlabel('x')
plt.ylabel('f(x)')

plt.show()
