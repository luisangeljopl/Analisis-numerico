import matplotlib.pyplot as plt
import numpy as np

n = int(input("Ingrese el numero de fila: \n"))
m = int(input("Ingrese el numero de columna: \n"))

matriz = np.random.randint(0, 100, size=(n, m))
print(matriz)
eigenvalues=np.linalg.eigvals(matriz)
print(eigenvalues)
x=(eigenvalues).real
y=(eigenvalues).imag
print(x)
print(y)
plt.scatter(x,y)
plt.xlabel('parte real')
plt.ylabel('parte imaginaria')
plt.show()
