import numpy as np
A = np.random.rand(6,7)
print(A)
print("Mostramos la primera y segunda columna", A[:, 0:2])
A = np.arange(15).reshape(3,5)
print("Mostramos la ultima fila", A[-1])
print("mostramos las ultimas 3 columnas", A[:, -3:])
print("tomamos una matriz menor dentro de A", A[0:3, 3:6])
