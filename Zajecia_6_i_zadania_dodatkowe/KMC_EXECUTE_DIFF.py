# -*- coding: utf-8 -*-
from __future__ import division
import os
import numpy as np
#os.chdir('C:/Users/Jan/Desktop/STUDIA_DOK/ROZTWORY/ZAJECIA_8')

#Program do symulacji KMC w warunkach 2D ukladu a`la NiAl
kB = 8.617*10**(-5)
T=1000
dim=30

#Zbiór potencjałów Isinga, V=wakancje
#VAB,VAA,VBB,VAV,VBV,VVV=(0,0,0,0,0,0)
VAB,VAA,VBB,VAV,VBV,VVV=(-0.125,-0.12,-0.06,0.04,-0.04,0)

#Potencjały oddziaływania z danym składnikiem zapisywane w formie wektorów
vec_int_A=np.array([VAA,VAB,VAV])
vec_int_B=np.array([VAB,VBB,VBV])
vec_int_V=np.array([VAV,VBV,VVV])

vec_int_tup=(vec_int_A,vec_int_B,vec_int_V)

from KMC_RTA_2D import *

#=============FUNCTIONS=======================
def get_period_hop(tup1,tup2,vec_period,dim):
    lim_crit=2*dim-1
    if tup1[0]%lim_crit!=0 and tup1[1]%lim_crit!=0:
        return vec_period
    x_det=tup2[0]-tup1[0]
    y_det=tup2[1]-tup1[1]
    if x_det==-lim_crit:
        vec_period[0]+=1
        if y_det==-lim_crit:
            vec_period[1]+=1
        elif y_det==lim_crit:
            vec_period[1]-=1
    elif x_det==lim_crit:
        vec_period[0]-=1  
        if y_det==-lim_crit:
            vec_period[1]+=1
        elif y_det==lim_crit:
            vec_period[1]-=1  
    else:
        if y_det==-lim_crit:
            vec_period[1]+=1
        elif y_det==lim_crit:
            vec_period[1]-=1  
    return vec_period

def get_distance_R_R2(tup0,tup_fin,dim,vec_period):
    suma_dis_R=[]
    suma_dis_R2=0
    for i in range(0,len(tup0)):
        suma_dis_R.append(tup_fin[i]+2*dim*vec_period[i]-tup0[i])
        suma_dis_R2+=(tup_fin[i]+2*dim*vec_period[i]-tup0[i])**2
    return suma_dis_R[0],suma_dis_R2
        
def diff_monitor(num_of_steps,dim,T,v_int_tup,num_of_trials):
    time_list=[]
    final_time_list=[]
    all_R_list=np.zeros(num_of_steps)
    all_R2_list=np.zeros(num_of_steps)
    #n_list=range(0,num_of_steps)
    d=gen_ordered_lattice(dim,1,0)
    d2=d.copy()
    sample=0
    while sample<num_of_trials:
        d=d2.copy()
        print sample
        time_MC=0
        R_list=[]
        R2_list=[]
        match_dict,group_prob_list,sum_prob=get_trans_probs(d,dim,T,v_int_tup)
        hop_vacancy,i=choose_trans(match_dict,group_prob_list)
        R_list.append(0)
        R2_list.append(0)
        tup0=hop_vacancy
        tup_last=hop_vacancy
        n=1
        vec_period=[0,0]
        time_MC+=(1/sum_prob)
        while n<num_of_steps:
            make_trans(d,hop_vacancy,i,dim)
            match_dict,group_prob_list,sum_prob=get_trans_probs(d,dim,T,v_int_tup)
            hop_vacancy,i=choose_trans(match_dict,group_prob_list)
            tup_fin=hop_vacancy
            vec_period=get_period_hop(tup_last,tup_fin,vec_period,dim)
            tup_last=hop_vacancy
            R_curr,R2_curr=get_distance_R_R2(tup0,tup_fin,dim,vec_period)
            R_list.append(R_curr)
            R2_list.append(R2_curr)
            n+=1
            time_MC+=(1/sum_prob)
        all_R_list=np.add(R_list,all_R_list)
        all_R2_list=np.add(R2_list,all_R2_list)
        sample+=1
        time_list.append(time_MC)
        print time_MC    
    mean_time=np.mean(time_list)
    step=mean_time/num_of_steps
    for i in np.arange(0,mean_time,step):
        final_time_list.append(i)
    two_mt_list=[]
    two_R2_list=[]
    #norm_R_list = [i/num_of_trials for i in all_R_list]
    norm_R2_list = [i/num_of_trials for i in all_R2_list]
    for i in range(0,len(final_time_list)-2):
        two_R2_list.append(norm_R2_list[i]+norm_R2_list[i+1])
        two_mt_list.append(final_time_list[i]+final_time_list[i+1])
    #plt.figure()
    #plt.title('Self diffusion, X displacement, num_of_trials='+str(num_of_trials))
    #plt.plot(n_list,norm_R_list,'.')
    #plt.figure()
    #plt.title('Self diffusion, averaged R2, num_of_trials='+str(num_of_trials))
    #plt.plot(n_list,norm_R2_list,'.')
    plt.figure()
    plt.title('Self diffusion, averaged R2, num_of_trials='+str(num_of_trials)+', T='+str(T)+'K')
    plt.plot(two_mt_list,two_R2_list,'.')
    plt.show()
    #satur_list_time=final_time_list[-int(len(norm_R2_list)/5):]
    #satur_list_R2=norm_R2_list[-int(len(norm_R2_list)/5):]
    satur_list_time=two_mt_list[:]
    satur_list_R2=two_R2_list[:]
    return satur_list_time, satur_list_R2

def estimate_D(time_list,r2_list):
    a,b=np.polyfit(time_list,r2_list,1)
    D=a/4
    return D
    
def diff_temp_analysis(dim,vec_int_tup):
    T_range=range(500,1100,100)
    D_list=[]
    for T in T_range:
        time_list,r2_list=diff_monitor(50,dim,T,vec_int_tup,1000)
        D=estimate_D(time_list,r2_list)
        D_list.append(D)
    #plt.title('Self diffusion, averaged R2, num_of_trials='+str(num_of_trials))
    print T_range
    print D_list
    plt.figure()
    plt.title('Diffusion constant temperature dependence')
    plt.xlabel('T')
    plt.ylabel('D')
    plt.plot(T_range,D_list,'.')
    plt.show()
    log_list=[]
    rev_T_list=[]
    for i in range(0,len(D_list)):
        logD=math.log(D_list[i])
        log_list.append(logD)
        rev_T=1/T_range[i]
        rev_T_list.append(rev_T)
    slope,lnD0=np.polyfit(rev_T_list,log_list,1)
    delta_E=-kB*slope
    plt.figure()
    plt.title('Arrhenius plot')
    plt.xlabel('1/T')
    plt.ylabel('log(D)')
    plt.plot(rev_T_list,log_list,'.')
    plt.show() 
    return delta_E, lnD0

def find_vac(d):
    for i in d.keys():
        if d[i][0]==vec_V[0]:
            if d[i][1]==vec_V[1]:
                print i
#print time_list, r2_list

def gauss_vac_anal(num_of_steps,dim,T,v_int_tup,num_of_trials):
    d=gen_vac_mid_str(dim)
    #print find_vac(d)
    d2=d.copy()
    sample=0
    x_pos_list=[]
    while sample<num_of_trials:
        time_MC=0
        d=d2.copy()
        #print find_vac(d), find_vac(d2)
        match_dict,group_prob_list,sum_prob=get_trans_probs(d,dim,T,v_int_tup)
        hop_vacancy,i=choose_trans(match_dict,group_prob_list)
        tup0=hop_vacancy
        tup_last=hop_vacancy
        time_MC+=(1/sum_prob)
        n=1
        vec_period=[0,0]
        while n<num_of_steps:
            make_trans(d,hop_vacancy,i,dim)
            match_dict,group_prob_list,sum_prob=get_trans_probs(d,dim,T,v_int_tup)
            hop_vacancy,i=choose_trans(match_dict,group_prob_list)
            tup_fin=hop_vacancy
            vec_period=get_period_hop(tup_last,tup_fin,vec_period,dim)
            tup_last=hop_vacancy
            time_MC+=(1/sum_prob)
            n+=1
        sample+=1
        x_pos_list.append(tup_last[0])
        print tup_last
    plt.figure()
    plt.title('Vac position after '+str(num_of_steps)+' steps')
    plt.xlabel('x_position')
    plt.ylabel('probability')
    plt.hist(x_pos_list, bins=range(min(x_pos_list), max(x_pos_list) + 1, 1))
    plt.show()
    return x_pos_list
    
def gauss_vac_anal_full(num_of_steps,dim,T,v_int_tup,num_of_trials):
    d=gen_vac_mid_str(dim)
    #print find_vac(d)
    d2=d.copy()
    sample=0
    x_pos_list=[]
    while sample<num_of_trials:
        d=d2.copy()
        #print find_vac(d), find_vac(d2)
        match_dict,group_prob_list,sum_prob=get_trans_probs(d,dim,T,v_int_tup)
        hop_vacancy,i=choose_trans(match_dict,group_prob_list)
        tup0=hop_vacancy
        tup_last=hop_vacancy
        n=1
        vec_period=[0,0]
        while n<num_of_steps:
            make_trans(d,hop_vacancy,i,dim)
            match_dict,group_prob_list,sum_prob=get_trans_probs(d,dim,T,v_int_tup)
            hop_vacancy,i=choose_trans(match_dict,group_prob_list)
            tup_fin=hop_vacancy
            vec_period=get_period_hop(tup_last,tup_fin,vec_period,dim)
            tup_last=hop_vacancy
            n+=1
        sample+=1
        x_pos_list.append(tup_last[0])
        print tup_last
    sample=0
    while sample<num_of_trials:
        d=d2.copy()
        #print find_vac(d), find_vac(d2)
        match_dict,group_prob_list,sum_prob=get_trans_probs(d,dim,T,v_int_tup)
        hop_vacancy,i=choose_trans(match_dict,group_prob_list)
        tup0=hop_vacancy
        tup_last=hop_vacancy
        n=1
        vec_period=[0,0]
        while n<num_of_steps+1:
            make_trans(d,hop_vacancy,i,dim)
            match_dict,group_prob_list,sum_prob=get_trans_probs(d,dim,T,v_int_tup)
            hop_vacancy,i=choose_trans(match_dict,group_prob_list)
            tup_fin=hop_vacancy
            vec_period=get_period_hop(tup_last,tup_fin,vec_period,dim)
            tup_last=hop_vacancy
            n+=1
        sample+=1
        x_pos_list.append(tup_last[0])
        print tup_last
    plt.figure()
    plt.title('Vac position after '+str(num_of_steps)+' steps')
    plt.xlabel('x_position')
    plt.ylabel('probability')
    plt.hist(x_pos_list, bins=range(0, 2*dim + 1, 1))
    plt.show()
    return x_pos_list

#=================================================

#ZADANIE 1

# Przy potencjalach dla układu modelującego NiAl ustawić T=500K, dim=8 i wykonać live_ordering_plot(60,d,dim,T,vec_int_tup). 

#ZADANIE 2

# Ustawic potencjaly Isinga dla ukladu uporzadkowanego, dim=30. Wykonac symulacje sledzenia pozycji Vac po num_of_steps krokach wywolujac funkcje:
# gauss_vac_anal(50,dim,1000,vec_int_tup,1000) (zajmie to chwilę)
# Zapisac obrazek. Nastepnie wykonac taka sama symulacje dla zerowych potencjalow Isinga. Zapisac i porownac obrazki. Jakie różnice można zaobserwować i z czego one wynikają? 
# W którym przypadku współczynnik dyfuzji przyjmuje wyższą wartość? Czy wszystkie spodziewane paski są obecne, z czego to wynika i jak to naprawić?

#ZADANIE 3

# Ustawic potencjaly Isinga dla ukladu uporzadkowanego, dim=30. Wykonac symulację pomiaru średniego kwadratu przemieszczenia w czasie w zależności od temperatury: diff_temp_analysis(dim,vec_int_tup). Jaki kształt ma obserwowany wykres? Ile wynosi energia aktywacji?
# Napisz osobny skrypt, w którym porównasz przebieg punktów otrzymanych w symulacji, z funkcją eksponencjalną Arrheniusa o delta_E i logD0 podanych na końcu outputu diff_temp_analysis.

#=================EXECUTIVES======================
#d=gen_ordered_lattice(dim,1,0)
#live_ordering_plot(60,d,dim,T,vec_int_tup)
#print diff_temp_analysis(dim,vec_int_tup)
#gauss_vac_anal(50,dim,1000,vec_int_tup,1000)
#=================================================