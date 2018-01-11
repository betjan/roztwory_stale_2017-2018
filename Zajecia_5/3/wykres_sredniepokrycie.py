import pylab

lista_x=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]
lista_y=[59.51661,68.60783,76.50208,82.69358,87.68669,91.21901,93.93834,95.74765,97.12576,98.00634]

pylab.xlim(0.00,0.11)
pylab.ylim(50.00,100.00)
pylab.title('zaleznosc sredniego pokrycia od energii')
pylab.xlabel('delta E(eV)')
pylab.ylabel('srednie pokrycie')

pylab.plot(lista_x,lista_y,'.')
pylab.show()