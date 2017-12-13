# -*- coding: utf-8 -*-

# ---SPECIFING LIBRARIES USED IN SCRIPT---
from __future__ import division
import os
# ---END---

# CHANGING TO PROGRAM DIRECTORY
os.chdir('C:/Users/Jan/Desktop/STUDIA_DOK/ROZTWORY') 

# LOADING ANOTHER SCRIPT TO INHERIT SOME FUNCTIONS
from ZADANIA_PYTHON_3 import *

# FUNCTION FOR CONFIGURATIONAL ENERGY CALCULATION
# INPUT: c_list= list of concentrations(len=4), vmatrix=matrix of Ising potentials, z=configurational number
# OUTPUT: E=configuration energy
def get_subl_ia(c_list,vmatrix,z=4):
    E=z*(c_list[0]*c_list[2]*vmatrix[0][2]+c_list[1]*c_list[3]*vmatrix[1][3]+c_list[0]*c_list[3]*vmatrix[0][3]+c_list[1]*c_list[2]*vmatrix[1][2]) 
    return E

# FUNCTION FOR FREE ENERGY PLOTTING VS. LRO PARAMETER (stoichiometry fixed)
# INPUT: T=temperature, vmatrix=matrix of Ising potentials, s=stoichiometry
# OUTPUT: FE_GRAPH
def get_Fline_L12(T,vmatrix,s,mesh=0.001,z=8):
    eta_list=[]
    lista_CA1=[]
    lista_CA2=[]
    lista_F=[]
    #(3/4)*ca1+(1/4)*ca2=s
    #ca2 = 4*s-3*ca1<1
    #ca1 > (4s-1)/3
    #ca1 < 1
    ca1_max= min((4/3)*s,1)
    ca1_min= max((4*s-1)/3,0)
    for ca1 in np.arange(ca1_min,ca1_max+mesh/2,mesh):
        ca2 = 4*s-3*ca1
        cb1 = 1-ca1
        cb2 = 1-ca2 
        c1=[ca1,cb1]
        c2=[ca2,cb2]
        c_list=c1+c2
        lista_CA1.append(ca1)
        lista_CA2.append(ca2)
        #S1=3*getEntropy(c1)/4
        H1=3*getHelmholtz(c1,vmatrix,z,T)/4
        S2=-T*getEntropy(c2)/4
        E12=3*get_subl_ia(c_list,vmatrix)/4
        TOTAL_F=H1+S2+E12
        #TOTAL_F=S1+S2
        lista_F.append(TOTAL_F)
        eta_list.append(get_eta(s,ca1))
        print ca1,cb1,ca2,cb2,TOTAL_F,H1,S2,E12
    #pylab.figure()
    #pylab.plot(eta_list,lista_F,'.')
    #pylab.plot(lista_CA2,lista_F,'.')
    #pylab.show()
    return (eta_list,lista_F)

# LONG-RANGE-ORDER PARAMETER ETA 
# INPUT: Tstart=lowest temperature, Tend=highest temperature, vmatrix=matrix of Ising potentials, T_mesh=temperature increment, s=stoichiometry
# OUTPUT: E=configuration energy
def get_T_eta_L12(Tstart,Tend,vmatrix,Tmesh,s):
    eta_min=[]
    T_list=[]
    for T in np.arange(Tstart,Tend,Tmesh):
        T_list.append(T)
        eta,F=get_Fline_L12(T,vmatrix,s)
        ind=get_min_ind(F)
        eta_min.append(eta[ind])
    pylab.figure()
    pylab.plot(T_list,eta_min,'.')
    #pylab.plot(lista_CA2,lista_F,'.')
    pylab.show()    
        

def get_eta(s,ca1):
    eta=2*(ca1-s)/(-ca1+2*s)
    return eta


vmatrix=np.array([[-0.12,-0.125,-0.12,-0.125],[-0.125,-0.06,-0.125,-0.06],[-0.12,-0.125,0,0],[-0.125,-0.06,0,0]])
#get_T_eta_L12(400,950,vmatrix,20,0.75)
#get_T_eta_L12(620,720,vmatrix,4,0.75)
#get_Fline_L12(570,vmatrix,0.75)


