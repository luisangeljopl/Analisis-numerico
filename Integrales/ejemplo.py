def funcion ( x ):
    f =2* x **3
    return ( f )
s =9
i =0
n =5
l = (s - i )/ n 
a =0
while i < s :
   r = funcion ( i )
   p = r * l
   a += p
   i += l
print ( " El area es " ,a )
