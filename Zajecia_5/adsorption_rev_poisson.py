# -*- coding: utf-8 -*-
from __future__ import division
import random
import math
import pylab
import numpy as np
import matplotlib.pyplot as plt

# Prosty program do testowania lańcuchów Markowa: równowagowa adsorpcja

#=========EXTERNAL_PARAMETERS==============
#Parametry zewnetrzne: stala Boltzmanna, temperatura, roznica energii miedzy poziomami: wezel pusty/wezel obsadzony; czastkowe czestosci przejsc (analogia procesu Poissona)
kB = 8.617*10**(-5)
T=298
delta_E=0.01 # Konwencja: delta_E>0
lamb_a=1
lamb_d=math.exp(-delta_E/(kB*T))
#==========================================

#=========SYSTEM_PARAMETERS================
#Poczatkowe pokrycie (ilosc wezlow)
start_cov=0
#Calkowita liczba wezlow
num_of_sites=100
#Calkowita liczba krokow
num_of_steps=100000
#==========================================

#==========FUNCTIONS=======================

#Funkcja do symulacji lancucha Markowa przez n-krokow
#OUTPUT: lista_numer_kroku, lista_czas_MC, lista_liczba_obsadzonych_wezlow
def simul_chain(start_cov,num_of_sites,num_of_steps):
    n=0
    t=0
    n_list=[]
    cov_list=[]
    t_list=[]
    n_list.append(n)
    t_list.append(t)
    cov_list.append(start_cov)
    while n<num_of_steps:
        uncov=num_of_sites-start_cov
        delta_t=1/(start_cov*lamb_d+uncov*lamb_a)
        toss=random.random()
        p_ads=uncov*lamb_a*delta_t
        if toss<p_ads:
            start_cov+=1
        else:
            start_cov-=1
        t+=delta_t
        n_list.append(n)
        t_list.append(t)
        cov_list.append(start_cov)
        #print t,start_cov
        n+=1
    return n_list,t_list,cov_list

#Funkcja do symulacji lancucha Markowa przez n-krokow
#OUTPUT: lista_prawdopodobienstw_stanow, lista_prawdopodobienstw_stanow_wazonych_czasem_MC, lista_czasow_przejsc_ze_stanow
def longterm_anal_chain(num_of_sites):
    # W tak malym ukladzie mozna dokonac obliczen analitycznych dla prawdopodobienstw dlugoterminowych
    delta_T_list=[]
    num_all_states=num_of_sites+1
    a_array=np.zeros((num_all_states,num_all_states))
    #print a_array
    a_array[0].fill(1)
    deltaT0=1/(num_of_sites*lamb_a)
    delta_T_list.append(deltaT0)
    #a_array[0][1]=lamb_d/(lamb_d+(num_of_sites-1)*lamb_a)
    a_array[num_of_sites][num_of_sites]=-1
    a_array[num_of_sites][num_of_sites-1]=lamb_a/(lamb_a+lamb_d*num_of_sites)
    #a_array[num_of_sites][num_of_sites-1]=lamb_a
    for i in range(1,num_of_sites):
        a_array[i][i-1]=lamb_a*(num_of_sites-(i-1))/(lamb_a*(num_of_sites-(i-1))+lamb_d*(i-1))
        #a_array[i][i-1]=lamb_a*(num_of_sites-(i-1))
        a_array[i][i]=-1
        #a_array[i][i+1]=lamb_d*(i+1)
        a_array[i][i+1]=lamb_d*(i+1)/(lamb_a*(num_of_sites-(i+1))+lamb_d*(i+1))
        deltaT=1/(lamb_a*(num_of_sites-i)+lamb_d*i)
        delta_T_list.append(deltaT)
    deltaTN=1/(num_of_sites*lamb_d)
    delta_T_list.append(deltaTN)
    b=np.zeros(num_all_states)
    b[0]=1
    print a_array
    probs=np.linalg.solve(a_array,b)
    unnorm_list=[]
    for i in range(0,len(delta_T_list)):
        z=delta_T_list[i]*probs[i]
        unnorm_list.append(z)
    #print probs
    norm_list = [i/sum(unnorm_list) for i in unnorm_list]
    #print norm_list
    return probs,norm_list,delta_T_list

#Funkcja do rysowania wyniku simul_chain()
def all_cov_plot(t_list,cov_list):
    pylab.figure()
    pylab.title('Simulation, T='+str(T)+', E='+str(delta_E))
    pylab.xlabel('Time[a.u.]')
    pylab.ylabel('Coverage')
    pylab.plot(t_list,cov_list)
    pylab.show()

#Funkcja do porownania prawdopodbienstw(kroku), bez wazenia czasem MC, w symulacji simul_chain()
# ,oraz teoretycznie
def long_term_dist_hist(cov_list,probs):
    satur_list=cov_list[-int(len(cov_list)/5):]
    up_prob=max(satur_list)
    low_prob=min(satur_list)
    prob_shortlist=probs[low_prob:up_prob+1]
    pylab.figure()
    pylab.hist(satur_list, bins=np.arange(low_prob-1/2, up_prob+1, 1), normed=True)
    pylab.figure()
    pylab.plot(range(low_prob,up_prob+1),prob_shortlist,'.')
    pylab.show()

#Funkcja do porownania prawdopodbienstw(kroku), bez wazenia czasem MC,
# w symulacji simul_chain(),oraz teoretycznie.
def long_term_t_weighted(cov_list,norm_list,delta_T_list):
    satur_list=cov_list[-int(len(cov_list)/5):]
    print satur_list
    up_prob=max(satur_list)
    low_prob=min(satur_list)
    print low_prob,up_prob
    #weigh_list=t_list[low_prob:up_prob+1]
    prob_shortlist=norm_list[low_prob:up_prob+1]
    pylab.figure()
    n_small_list=range(low_prob,up_prob+1)
    #hist_dict=collections.Counter(satur_list)
    #hist_list=hist_dict.values()
    hist_list=[]
    for n in n_small_list:
        prob=satur_list.count(n)
        freq=delta_T_list[n]*prob
        hist_list.append(freq)
    norm_hist_list = [i/sum(hist_list) for i in hist_list]
    pylab.bar(n_small_list,norm_hist_list,align='center')
    pylab.title('Simulation: frequency of states, T='+str(T)+'K, E='+str(delta_E))
    pylab.xlabel('Coverage')
    pylab.ylabel('Relative frequency')
    pylab.figure()
    pylab.plot(range(low_prob,up_prob+1),prob_shortlist,'.')
    pylab.title('Calculation: frequency of states, T='+str(T)+', E='+str(delta_E))
    pylab.xlabel('Coverage')
    pylab.ylabel('Relative frequency')
    pylab.show()

#Zmiany rozkladu prawdopodobienstw wystepowania stanow(teoretyczne), rysowane na zywo
#UWAGA: Nie uzywac dla ukladow wiekszych niz N=100 ze wzgledu na dlugi czas procesu
def live_plot_ME(num_of_sites,start_cov):
    delta_T_list=[]
    num_all_states=num_of_sites+1
    trans_matrix=np.zeros((num_all_states,num_all_states))
    deltaT0=1/(num_of_sites*lamb_a)
    delta_T_list.append(deltaT0)
    trans_matrix[0][1]=lamb_d/(lamb_d+(num_of_sites-1)*lamb_a)
    trans_matrix[num_of_sites][num_of_sites-1]=lamb_a/(lamb_a+lamb_d*num_of_sites)
    for i in range(1,num_of_sites):
        trans_matrix[i][i-1]=lamb_a*(num_of_sites-(i-1))/(lamb_a*(num_of_sites-(i-1))+lamb_d*(i-1))
        #a_array[i][i-1]=lamb_a*(num_of_sites-(i-1))
        #a_array[i][i+1]=lamb_d*(i+1)
        trans_matrix[i][i+1]=lamb_d*(i+1)/(lamb_a*(num_of_sites-(i+1))+lamb_d*(i+1))
        deltaT=1/(lamb_a*(num_of_sites-i)+lamb_d*i)
        delta_T_list.append(deltaT)
    deltaTN=1/(num_of_sites*lamb_d)
    delta_T_list.append(deltaTN)
    prob_vector=np.zeros(num_all_states)
    prob_vector[start_cov]=1
    prob_ref_vector=np.zeros(num_all_states)
    n=0
    #pylab.figure()
    plt.axis([0, num_of_sites, 0, 1])
    plt.ion()
    plt.scatter(range(0,num_all_states),prob_vector)
    plt.ylim([0,1])
    plt.xlim([0,num_of_sites])
    plt.pause(0.2)
    plt.clf()
    while True:
        prob_vector=np.dot(trans_matrix,prob_vector)
        #print prob_vector
        plt.scatter(range(0,num_all_states),prob_vector)
        plt.ylim([0,1])
        plt.xlim([0,num_of_sites])
        plt.xlabel('Coverage')
        plt.ylabel('Probability')
        plt.title('Probability distribution, n='+str(n))
        plt.pause(0.001)
        if np.linalg.norm(np.subtract(prob_vector,prob_ref_vector))<5e-5:
            break
        plt.clf()
        if n%2==0:
            prob_ref_vector=prob_vector
        n+=1
    prob_ref_vector=np.dot(trans_matrix,prob_vector)
    unperiod_vector=np.add(prob_vector,prob_ref_vector)
    for i in range(0,len(unperiod_vector)):
        unperiod_vector[i]/=2
    plt.clf()
    pylab.title('Calculation: averaged probabilty of states, T='+str(T)+', E='+str(delta_E))
    plt.xlabel('Coverage')
    plt.ylabel('Probability')
    plt.plot(range(0,num_all_states),unperiod_vector,'.')
    plt.ylim([0,1])
    plt.xlim([0,num_of_sites])
    pylab.show()

#==========================================


#============EXECUTIVES====================
n_list,t_list,cov_list=simul_chain(start_cov,num_of_sites,num_of_steps)
#all_cov_plot(t_list,cov_list)
#all_cov_plot(n_list,cov_list)
probs,norm_list,delta_T_list=longterm_anal_chain(num_of_sites)
#long_term_dist_hist(cov_list,probs)
long_term_t_weighted(cov_list,norm_list,delta_T_list)
#live_plot_ME(num_of_sites,start_cov)
#==========================================
