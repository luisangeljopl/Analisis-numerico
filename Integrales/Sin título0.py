# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 21:39:27 2019

@author: Luis Angel
"""

import numpy as np
import matplotlib.pyplot as plt

#Función a evaluar
def funcionx(x):
    fxi = x**2
    return(fxi)

# PROGRAMA
a = 1 # float(input('a: '))
b = 3 # float(input('b: '))
n = 1000 # int(input('tramos: '))
ra=1
rb=3

# PROCEDIMIENTO
h = (b-a)/n
x = a
suma = funcionx(x)
for i in range(0,n-1,1):
    x = x+h
    suma = suma + 2*funcionx(x)
suma = suma + funcionx(b)
area = h*(suma/2)

# SALIDA
print('tramos: ', n)
print('Integral: ', area)

def riemannplot(f, a, b, ra, rb, n):
    # f es la función 
    # a y b son los limites del eje x para graficar la funcion f
    # ra y rb son los limites del intervalo en el eje x del que queremos calcular la suma
    # n es el numero de rectangulos que calcularemos
    atenuacion = (b-a)/n
    x = np.arange(a, b+atenuacion, atenuacion)
    plt.plot(x, f(x), color='red')
    n = 1000 # int(input('tramos: '))
    delta_x = (rb-ra)/n
   
    riemannx = np.arange(ra, rb, delta_x)
    riemanny = f(riemannx)
    plt.bar(riemannx,riemanny,width=delta_x,alpha=0.5,facecolor='red')
    plt.plot(x,f(x), 'o')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('n=1000')
    plt.figtext(0.3,0.7, "A= " + str(area) , color='blue')
    plt.show()



def f(x):
    return x**2

riemannplot(f, 0, 1.1, 0, 1, 10)