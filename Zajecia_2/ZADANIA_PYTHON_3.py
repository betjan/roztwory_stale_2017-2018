# -*- coding: utf-8 -*-
# Zestaw III: Przybliżenie Bragga-Williamsa

# Ten moduł elminuje problem z dzieleniem int/int
from __future__ import division
import math
# Biblioteka do obliczeń nuemrycznych
import numpy as np
import os
import pylab

# Wartość referencyjna, slużąca do sprowadzania floatów do zera
# oraz eliminowaniu problemów z math.log(0)
eps=10**(-23)
# Stała Boltzmanna wyrażona w eV
kB = 8.617*10**(-5)

# Funkcja do obliczania entropii dla zadanej listy stężeń clist
# INPUT: clist - lista stężeń 
# OUTPUT: entropia(float) 
def getEntropy(clist):
    suma_S=0
    for c in clist:
        if c==0:
            suma_S+=eps*math.log(eps)
        else:
            suma_S+=c*math.log(c)
    return -kB*suma_S

# Funkcja do obliczania energii konfiguracyjnej w modelu Isinga
# INPUT: clist - lista stężeń, vmatrix - macierz potencjałów dwuciałowych, z - liczba koordynacyjna   
# OUTPUT: energia(float)
def getEnergy(clist,vmatrix,z):
    suma_E=0
    for i in range(0,len(clist)):
        suma_E+=clist[i]**2*vmatrix[i][i]/2 
        for j in range(i+1,len(clist)):
            suma_E+=clist[i]*clist[j]*vmatrix[i][j] 
    return suma_E*z

# Funkcja do obliczania energii swobodnej (argumenty jak wyżej)
def getHelmholtz(clist,vmatrix,z,T):
    return getEnergy(clist,vmatrix,z)-T*getEntropy(clist)

# Funkcja do generacji listy stężeń - jest to tzw. generator, który wymaga specjalnej inicjalizacji
# przy użyciu pętli for, np. "for clist in Clist_gen(4,0,0.02):", aby zwracał kolejne kombinacje stężeń
# Co więcej, jest to generator rekurencyjny, tj. funkcja zleca wykonywanie samej siebie.
# Ciężko się połapać jak to działa i nie jest to konieczne. Jest to również jedna z możliwych odpowiedzi
# na zadanie bardziej specjalne z poprzedniego zestawu        
def Clist_gen(cnum,history,mesh):
    if cnum>1:
        for i in np.arange(0,1-history+0.00001,mesh):
            history+=i
            cnum-=1
            for value in Clist_gen(cnum,history,mesh):
                yield [i]+value
            history-=i
            cnum+=1
    else:
        if abs(1-history)<(mesh/2):
            yield [0]
        else:
            yield [1-history]

for clist in Clist_gen(4,0,0.02):
    print clist

#Zaczniemy od najprostszej formy przybliżenia BW, gdzie nie dokonuje się żadnego rozróżnienia między węzłami sieci
#W oczywisty sposób, jest to podejście, w którym ilość informacji otrzymywanej o układzie jest minimalna. 
#Klasycznym i elementarnym zastosowaniem przybliżenia BW w tej formie jest układ spinów (+1/-1) w sztywnej sieci.
#Przyjmujemy oddziałwyanie wyłacznie między najbliższymi sąsiadami (NN - nearest neighbours), -J dla zgodnej orientacji spinów, +J dla przeciwnej 
#(przy takiej konwencji J>0 dla ferromagnetyka). Wiadomo, że układ taki jest dwuskładnikowy i mamy stężenia (c+,c-). Obserwablę stanowi 
#magnetyzacja próbki M, proporcjonalna do (c+)-(c-). W ramach tutorialu, zbadamy zachowanie takiego układu 
#bez przyłożonego zewnętrznego pola magnetycznego B, w funkcji temperatury T. 

#W pierwszym kroku, zobaczmy, jak wygląda zależność F(c+) dla pojedynczej temperatury
#Przy ustalanej stałej sprzężenia J; definiujemy funkcję:

def show_F_diag(T,J,z):
    lista_C=[]
    lista_F=[]
    #Macierz oddziaływań wprowadzana przy użyciu struktury numpy.array()
    vmatrix=np.array([[-J,J],[J,-J]])
    for c in Clist_gen(2,0,0.001):
        lista_C.append(c[0])
        lista_F.append(getHelmholtz(c,vmatrix,z,T))
    #Rysowanie wykresu lista_F(lista_C) przy użyciu punktów (".")
    pylab.plot(lista_C,lista_F)
    pylab.xlabel('c+')
    pylab.ylabel('F')
    pylab.title("Energia Helmholtza, "+"T="+str(T) + ", J=" + str(J) + ", z=" + str(z))
    pylab.show()
                     
#show_F_diag(1000,0.01,6)   

#Ponieważ będzie to potrzebne w późniejszych zadaniach, definiujemy prostą funkcję
# do znajdowania najmniejszej wartości, która zwraca odpowiadający jej indeks z listy                                                                                                                                                                                                                                                                                                                                   
def get_min_ind(Helm_list):
    mem_helm=Helm_list[0]
    mem_ind=0
    for i in range(0,len(Helm_list)):
        test_helm=Helm_list[i]
        if test_helm<(1.00000001*mem_helm):
            mem_helm=test_helm
            mem_ind=i
    return mem_ind
    
#------PPZYKLADOWE ZADANIE: PRZYBLIŻENIE BW DLA SPINÓW-------     

#_1)Przepisz funkcję show_F_diag do nowej funkcji get_F_diag, tak, aby zwracała (lista_C,lista_F) i nic nie rysowała
#_2)Napisz nową funkcję get_T_magn(Tstart,Tend,step,J,z), która, iterując po range(Tstart,Tend,step), będzie wykonywać get_F_diag(),
#   przy pomocy get_min_ind() wyznaczać położenie minimum "na wykresie" (czyli c+ odpowiadające najniższej F), obliczać "magnetyzację" M
#   dla tego punktu (właściwie różnicę populacji spinów) i zwracać listę temperatur (lista_T) oraz listę M(T) (lista_mag): return (lista_T,lista_mag)
#_3)Dobierz parametry tak, aby widoczne było przejście ferromagnetyk-paramagnetyk i narysuj wykres przy pomocy plot_mag_T(lista_T,lista_mag).                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                            
#-----------------------                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

#------ROZWIĄZANIE ZADANIA------
def get_F_diag(T,J,z):
    lista_C=[]
    lista_F=[]
    vmatrix=np.array([[-J,J],[J,-J]])
    for c in Clist_gen(2,0,0.001):
        lista_C.append(c[0])
        lista_F.append(getHelmholtz(c,vmatrix,z,T))
    return (lista_C,lista_F)
    
def get_T_magn(Tstart,Tend,step,J,z):
    lista_T=[]
    lista_mag=[]
    for t in range(Tstart,Tend,step):
        lista_T.append(t)
        tup=get_F_diag(t,J,z)
        ind=get_min_ind(tup[1])
        c=tup[0][ind]
        lista_mag.append(abs(2*c-1))   
    return (lista_T,lista_mag)   
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  