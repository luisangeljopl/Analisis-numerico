# Integración: Regla de los trapecios
import numpy as np
import matplotlib.pyplot as plt

#Función a evaluar
def funcionx(x):
    fxi = x**2
    return(fxi)

# PROGRAMA
a =float(input('a: '))
b = float(input('b: '))
tramos = int(input('tramos: '))

# PROCEDIMIENTO
h = (b-a)/tramos
x = a
suma = funcionx(x)
for i in range(0,tramos-1,1):
    x = x+h
    suma = suma + 2*funcionx(x)
suma = suma + funcionx(b)
area = h*(suma/2)

# SALIDA
print('tramos: ', tramos)
print('Integral: ', area)
# Puntos para la gráfica
muestras = tramos+1
xi = np.linspace(a,b,muestras)
fi = funcionx(xi)

# Gráfica
# Referencia función contínua
xig = np.linspace(a,b,muestras*10)
fig = funcionx(xig)
plt.plot(xig,fig)
# Trapecios
plt.fill_between(xi,0,fi, color='g')
plt.title('Integral: Regla de Trapecios')
for i in range(0,muestras,1):
    plt.axvline(xi[i], color='w')
plt.plot(xi,fi, 'o',) # puntos muestra
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend(loc="Upper right")
plt.show()
