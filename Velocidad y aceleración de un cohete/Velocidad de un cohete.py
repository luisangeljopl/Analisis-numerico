import numpy as np
import matplotlib.pylab as plt
import math

def g(x,y):
    return -9.8+(2000/(200-x))+(2/(x-200))*y


h=0.1
s=10
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
plt.title("VELOCIDAD DE UN COHETE")
plt.plot( x ,y , linewidth=1.0)
plt.xlabel("t")
plt.ylabel("v(t)")
plt.ylim(0,15)
plt.show()
