# -*- coding: utf-8 -*-
from __future__ import division
import os

os.chdir('C:/Users/Jan/Desktop/STUDIA_DOK/ROZTWORY') 

from ZADANIA_PYTHON_3 import *

def get_Fline_B2(T,vmatrix,s,mesh=0.001,z=8):
    eta_list=[]
    lista_CA1=[]
    lista_CA2=[]
    lista_F=[]
    #(3/4)*ca1+(1/4)*ca2=s
    #ca2 = 4*s-3*ca1<1
    #ca1 > (4s-1)/3
    for ca1 in np.arange(max(2*s-1,0),min(1,2*s),mesh):
        ca2 = 2*s-ca1
        cb1 = 1-ca1
        cb2 = 1-ca2 
        c1=[ca1,cb1]
        c2=[ca2,cb2]
        c_list=c1+c2
        lista_CA1.append(ca1)
        lista_CA2.append(ca2)
        H=getHelmholtz(c_list,vmatrix,z,T)/2
        lista_F.append(H)
        eta_list.append(get_eta_B2(s,ca1))
    #pylab.figure()
    #pylab.plot(eta_list,lista_F,'.')
    #pylab.plot(lista_CA2,lista_F,'.')
    #pylab.show()
    return (eta_list,lista_F)

def get_T_eta_B2(Tstart,Tend,vmatrix,Tmesh,s):
    eta_min=[]
    T_list=[]
    for T in np.arange(Tstart,Tend,Tmesh):
        T_list.append(T)
        eta,F=get_Fline_B2(T,vmatrix,s)
        ind=get_min_ind(F)
        eta_min.append(eta[ind])
    pylab.figure()
    pylab.plot(T_list,eta_min,'.')
    #pylab.plot(lista_CA2,lista_F,'.')
    pylab.show()    
        
def get_eta_B2(s,ca1):
    eta=abs((ca1-s)/s)
    return eta

vmatrix=np.array([[0,0,-0.12,-0.125],[0,0,-0.125,-0.06],[-0.12,-0.125,0,0],[-0.125,-0.06,0,0]])
#get_T_eta_B2(200,2000,vmatrix,100,0.5)
#get_T_eta_B2(1600,1650,vmatrix,2,0.5)
#get_Fline_B2(1626,vmatrix,0.5)
