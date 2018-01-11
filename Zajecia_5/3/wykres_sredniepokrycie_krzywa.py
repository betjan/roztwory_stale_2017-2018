import pylab
import math
import numpy as np
kB = 8.617*10**(-5)
T=298
lista_x=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]
lista_y=[59.51661,68.60783,76.50208,82.69358,87.68669,91.21901,93.93834,95.74765,97.12576,98.00634]

def func(x):
    return 100*(1/(math.exp(-x/(kB*T))+1))
    #return 100*(1/(math.exp(-x/(8.617*10**(-5)*298))+1))
    
y_list=[]
a=np.arange(0,0.5,0.001)
for i in a:
    y_list.append(func(i))
print a
    

pylab.xlim(0.00,0.11)
#pylab.ylim(50.00,100.00)
pylab.title('zaleznosc sredniego pokrycia od energii')
pylab.xlabel('delta E(eV)')
pylab.ylabel('srednie pokrycie')

pylab.plot(lista_x,lista_y,'.')
pylab.plot(a,y_list,'-')
pylab.show()
