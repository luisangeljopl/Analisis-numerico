n=4999
if n > 0:

    lista=[]

    for i in range(2,n):
        

        creciente = 2

        esPrimo = True

        while esPrimo and creciente < i:

            if i % creciente == 0:

                esPrimo = False

            else:

                creciente += 1

        if esPrimo:

            

            
            print(i)
            



