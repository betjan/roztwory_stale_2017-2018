import pylab
import math
import numpy as np

lista_x=[0.01,0.02,0.03,0.04,0.05,0.06]
lista_y=[377,430,481,529,573,614]

a=np.polyfit(lista_x,lista_y,2)
print a

x2=a[0]
x1=a[1]
x0=a[2]

def func(x,a,b,c):
    return a*x**2+b*x+c


'''
def func(x):
    return math.ln(x)
'''   
y_list=[]
a=np.arange(0.01,0.061,0.002)
for i in a:
    y_list.append(func(i,x2,x1,x0))
print a

pylab.xlabel('Delta E')
pylab.ylabel('Liczba krokow')
pylab.xlim(0,0.07)
pylab.title('Dopasowana funkcja kwadratowa')
pylab.plot(lista_x,lista_y,'.')
pylab.plot(a,y_list,'-')
pylab.show()
