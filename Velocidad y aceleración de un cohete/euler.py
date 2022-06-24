import numpy as np
import matplotlib.pylab as plt
import math
from pylab import *

def f(x):
    return (3/125)*x**2+(1/5)*x

def g(x,y):
    return -9.8+(2000/(200-x))+(2/(x-200))*y

x = range(0, 10)
h=0.1
s=2.5
N = int((s/h)+1)
x = np.zeros(N)
y=np.zeros(N)
x[0]=0
y[0]=0
print(x[0],y[0],)
for i in np.arange(1,N):
    y[i]=y[i-1]+(g(x[i-1],y[i-1]))*h
    x[i]=x[i-1]+h
    print(x[i],y[i])

# construcci´on de gr´afica
plt.title("Comparación entre gráficas")
plt.plot( x ,y , linewidth=1.0)
plt.xlabel("tiempo")
plt.ylabel("velocidad")
plt.plot(x, [f(i) for i in x],color="red")
plt.plot(x,y, '-r', label='Gráfica por método númerico',color="blue")
plt.plot(x, [f(i) for i in x], '-r', label='Función analítica')
plt.legend(loc='upper left')

# Mostramos el grafico.
plt.show()
