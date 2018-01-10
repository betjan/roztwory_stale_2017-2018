# -*- coding: utf-8 -*-
from __future__ import division
import random
import math
import pylab
import numpy as np
import matplotlib.pyplot as plt

# Prosty program do testowania lańcuchów Markowa: równowagowa adsorpcja z wykorzystaniem algorytmu Metropolisa

#=========EXTERNAL_PARAMETERS==============
#Parametry zewnetrzne: stala Boltzmanna, temperatura, roznica energii miedzy poziomami: wezel pusty/wezel obsadzony
kB = 8.617*10**(-5)
T=298
delta_E=0.02 # Konwencja: delta_E>0 oznacza preferencję stanu obsadzonego
#==========================================

#=========SYSTEM_PARAMETERS================
#Poczatkowe pokrycie (w jednostkach liczby węzłów)
start_cov=18
#Calkowita liczba wezlow
num_of_sites=100
#Calkowita liczba krokow
num_of_steps=500000
#==========================================

#back_prob = elementarne prawdopodobieństwo odczepienia zaadsorbowanej cząsteczki
back_prob=math.exp(-delta_E/(kB*T))
#stay_prob = elementarne prawdopodobieństwo pozostania w obecnym stanie pokrycia
stay_prob=1-back_prob
#Pytanie: Ile wynosi elementarne prawdopodobieństwo adsorpcji (przyłączenia) cząsteczki?


#==========FUNCTIONS=====================
def simul_chain(start_cov,num_of_sites,num_of_steps):
    n=0
    n_list=[]
    cov_list=[]
    n_list.append(n)
    cov_list.append(start_cov)
    while n<num_of_steps:
        coverred = start_cov/num_of_sites
        if start_cov==num_of_sites:
            arrow='back'
            toss=random.random()
        elif start_cov==0:
            arrow='up'
        else:
            can_down = random.random()
            if can_down<coverred:
                arrow='back'
                toss=random.random()
            else:
                arrow='up'
            #arrow=random.choice(('up','back'))
        if arrow=='up':
            start_cov+=1
        else:
            if toss<back_prob:
                start_cov-=1
        cov_list.append(start_cov)
        n_list.append(n)
        n+=1
    return n_list, cov_list

#Funkcja do analitycznego obliczania prawdopodobieństw długoterminowych
def longterm_anal_chain(num_of_sites):
    num_all_states=num_of_sites+1
    a_array=np.zeros((num_all_states,num_all_states))
    a_array[0].fill(1)
    a_array[num_of_sites][num_of_sites]=stay_prob
    a_array[num_of_sites][num_of_sites-1]=1/num_of_sites
    a_array[1][0]=1
    a_array[1][1]=stay_prob*(1/num_of_sites)-1
    a_array[1][2]=back_prob*(2/num_of_sites)
    a_array[num_of_sites-1][num_of_sites-2]=(2/num_of_sites)
    a_array[num_of_sites-1][num_of_sites-1]=stay_prob*((num_of_sites-1)/num_of_sites)-1
    a_array[num_of_sites-1][num_of_sites]=back_prob
    for i in range(2,num_of_sites-1):
        a_array[i][i-1]=((num_of_sites-(i-1))/num_of_sites)
        a_array[i][i]=stay_prob*(i/num_of_sites)-1
        a_array[i][i+1]=back_prob*((i+1)/num_of_sites)
    b=np.zeros(num_all_states)
    b[0]=1
    #print a_array
    probs=np.linalg.solve(a_array,b)
    return probs

def all_cov_plot(n_list,cov_list):
    pylab.figure()
    pylab.title('Simulation, T='+str(T)+', E='+str(delta_E))
    pylab.xlabel('Step')
    pylab.ylabel('Coverage')
    pylab.plot(n_list,cov_list)
    pylab.show()

def get_satur_params(cov_list):
    satur_list=cov_list[-int(len(cov_list)/5):]
    cov_mean=np.mean(satur_list)
    std_dev=np.std(satur_list)
    return cov_mean,std_dev

def long_term_dist_hist(cov_list,probs):
    satur_list=cov_list[-int(len(cov_list)/5):]
    up_prob=max(satur_list)
    low_prob=min(satur_list)
    prob_shortlist=probs[low_prob:up_prob+1]
    pylab.figure()
    pylab.hist(satur_list, bins=np.arange(low_prob-1/2, up_prob+1, 1), normed=True)
    pylab.title('Simulation: frequency of states, T='+str(T)+', E='+str(delta_E))
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
    num_all_states=num_of_sites+1
    trans_matrix=np.zeros((num_all_states,num_all_states))
    trans_matrix[num_of_sites][num_of_sites]=stay_prob
    trans_matrix[num_of_sites][num_of_sites-1]=1/num_of_sites
    trans_matrix[1][0]=1
    trans_matrix[1][1]=stay_prob*(1/num_of_sites)
    trans_matrix[1][2]=back_prob*(2/num_of_sites)
    trans_matrix[num_of_sites-1][num_of_sites-2]=(2/num_of_sites)
    trans_matrix[num_of_sites-1][num_of_sites-1]=stay_prob*((num_of_sites-1)/num_of_sites)
    trans_matrix[num_of_sites-1][num_of_sites]=back_prob
    for i in range(2,num_of_sites-1):
        trans_matrix[i][i-1]=((num_of_sites-(i-1))/num_of_sites)
        trans_matrix[i][i]=stay_prob*(i/num_of_sites)
        trans_matrix[i][i+1]=back_prob*((i+1)/num_of_sites)
    prob_vector=np.zeros(num_all_states)
    prob_vector[start_cov]=1
    prob_ref_vector=np.zeros(num_all_states)
    n=0
    plt.axis([0, num_of_sites, 0, 1])
    plt.ion()
    plt.scatter(range(0,num_all_states),prob_vector)
    plt.ylim([0,1])
    plt.xlim([0,num_of_sites])
    plt.xlabel('Coverage')
    plt.ylabel('Probability')
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
        prob_ref_vector=prob_vector
        n+=1
#====================================================

#===============EXECUTIVES======================
#print longterm_anal_chain(num_of_sites)
#n_list,cov_list=simul_chain(start_cov,num_of_sites,num_of_steps)
#all_cov_plot(n_list,cov_list)
#probs=longterm_anal_chain(num_of_sites)
#long_term_dist_hist(cov_list,probs)
#long_term_t_weighted(cov_list,norm_list,delta_T_list) tego nie lubimy
#mean_cov,std_dev=get_satur_params(cov_list)
#print mean_cov,std_dev
#live_plot_ME(num_of_sites,start_cov)
#===============EXECUTIVES======================
