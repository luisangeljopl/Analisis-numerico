from matplotlib.pylab import*
# ahora con velocidad variable

x0=0
delta_t=0.01
puntos=200

# se define la funci´on velocidad
def v(t):
    return 3.4*t**2

# un par de arreglos
t=[n*delta_t for n in range (puntos)]
x=[0 for n in range(puntos)]
x[0]= x0
for n in range(puntos-1):
     x[n+1]=x[n]+delta_t*v(t[n]) # la f´ormula ()
print(t,x)
# construcci´on de gr´afica
title("MOVIMIENTO ACELERADO")
xlabel("t")
ylabel("x(t)")
plot(t,x,linewidth=2.0)
#plt.savefig(’../fig-ch1/euler_simple.pdf’, transparent=True)
show()
