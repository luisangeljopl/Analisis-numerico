# coding: utf-8
# Your code here!

def funcion(x):
 f=x*2
 return(f)
s=15 #límite superior
i=0 #límite inferior
n=75 #particiones
l = (s-i)/n  #largo de partición
a=0
while i < s :
    r = funcion(i) #evaluar la funcion
    p = r * l 
    a += p
    i += l 

print("El area es",a)
