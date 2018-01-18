# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import math
import random
import matplotlib.pyplot as plt

#Program do symulacji KMC w warunkach 2D ukladu a`la NiAl
#T=500
#dim=30

#Zbiór potencjałów Isinga, V=wakancje
#VAB,VAA,VBB,VAV,VBV,VVV=(0,0,0,0,0,0)
#VAB,VAA,VBB,VAV,VBV,VVV=(-0.125,-0.12,-0.06,0.04,-0.04,0)
#VAB,VAA,VBB,VAV,VBV,VVV=(-0.12,-0.06,-0.06,0,0,0)

#Potencjały oddziaływania z danym składnikiem zapisywane w formie wektorów
#vec_int_A=np.array([VAA,VAB,VAV])
#vec_int_B=np.array([VAB,VBB,VBV])
#vec_int_V=np.array([VAV,VBV,VVV])

kB = 8.617*10**(-5)
#Wektory jednostkowe reprezentujące obsadzenie pojedynczego wezla
vec_A=np.array([1,0,0])
vec_B=np.array([0,1,0])
vec_V=np.array([0,0,1])
vec_empty=np.array([0,0,0])

#Funkcja generująca "powierzchnię" BCC ("spłaszczone" dwie warstwy)
#z zadaną liczbą wakancji i stechiometrią N_A/(N_A+N_B)
def gen_random_lattice(dim,stech,cv_perc):
    d={}
    num_of_sites=2*dim**2
    #print num_of_sites
    n=0
    for x in range(0,dim*2,2):
        for y in range(0,dim*2,2):
            d[(x,y)]=vec_V
            d[(x+1,y+1)]=vec_V
            n+=2
    site_list=d.keys()
    num_A=int((1-cv_perc)*stech*num_of_sites)
    num_B=int((1-cv_perc)*(1-stech)*num_of_sites)
    n=0
    while n<num_A:
        p=random.choice(site_list)
        d[p]=vec_A
        site_list.remove(p)
        n+=1
    n=0
    while n<num_B:
        p=random.choice(site_list)
        d[p]=vec_B
        site_list.remove(p)
        n+=1
    return d

def gen_ordered_lattice(dim,nv_latt_a,nv_latt_b):
    d={}
    #print num_of_sites
    n=0
    a_list=[]
    b_list=[]
    for x in range(0,dim*2,2):
        for y in range(0,dim*2,2):
            a_list.append((x,y))
            b_list.append((x+1,y+1))
            d[(x,y)]=vec_A
            d[(x+1,y+1)]=vec_B
            n+=2
    vac_a=random.sample(a_list, nv_latt_a)
    vac_b=random.sample(b_list, nv_latt_b)
    for i in vac_a:
        d[i]=vec_V
    for j in vac_b:
        d[j]=vec_V
    return d
    
def gen_vac_mid_str(dim):
    d={}
    #print num_of_sites
    n=0
    a_list=[]
    b_list=[]
    for x in range(0,dim*2,2):
        for y in range(0,dim*2,2):
            a_list.append((x,y))
            b_list.append((x+1,y+1))
            d[(x,y)]=vec_A
            d[(x+1,y+1)]=vec_B
            n+=2
    vac_a=[(dim,dim)]
    for i in vac_a:
        d[i]=vec_V
    return d

def index_move_translation(index):
    if index==0:
        return [-1,-1]
    elif index==1:
        return [-1,1]
    elif index==2:
        return [1,-1]
    else:
        return [1,1]

def get_vector_type(at,(vec_int_A,vec_int_B,vec_int_V)):
    if at[0]==1:
        return vec_int_A
    elif at[1]==1:
        return vec_int_B
    else:
        return vec_int_V

def get_neigh(dic,point,dim):
    l=[]
    neighs=[]
    for i in [-1,1]:
        for j in [-1,1]:
            neigh_coord=((i+point[0])%(2*dim),(j+point[1])%(2*dim))
            neighs.append(neigh_coord)
            l.append(dic[neigh_coord])
    return np.sum(l,axis=0),neighs

def get_energy(dic,p,dim,v_int_tup):
    at=dic[p]
    e=np.dot(get_vector_type(at,v_int_tup),get_neigh(dic,p,dim)[0])
    return e

def get_hipo_energy(dic,p,hipo_at_p,dim,v_int_tup):
    e=np.dot(get_vector_type(hipo_at_p,v_int_tup),get_neigh(dic,p,dim)[0])
    return e

def get_total_energy(dic,dim,v_int_tup):
    e=0
    for i in dic.keys():
        e+=get_energy(dic,i)
    return e/2

def get_vac_list(dic):
    vac_list=[]
    for k in dic.keys():
        val=dic[k]
        if val[2]==1:
            vac_list.append(k)
    return vac_list

def get_single_trans_prob(dic,vac_point,sum_prob,dim,T,v_int_tup):
    neigh_vect,l=get_neigh(dic,vac_point,dim)
    #print neigh_vect,l
    e_0_start=get_energy(dic,vac_point,dim,v_int_tup)
    #print e_0_start
    e_start_list=[]
    e_end_list=[]
    e_0_end_list=[]
    for p in l:
        e_start_list.append(get_energy(dic,p,dim,v_int_tup)+e_0_start)
        #print get_energy(dic,p),e_0_start
        e_end_list.append(get_hipo_energy(dic,p,vec_V,dim,v_int_tup))
        vec_p=dic[p]
        e_0_end_list.append(get_hipo_energy(dic,vac_point,vec_p,dim,v_int_tup))
        #print e_start_list
    delta_E_list=np.subtract(np.add(e_0_end_list,e_end_list),e_start_list)
    prob_list=[]
    for e in delta_E_list:
        prob=math.exp(-e/(2*kB*T))
        prob_list.append(prob+sum_prob)
        sum_prob+=prob
    return prob_list
    
def get_trans_probs(dic,dim,T,v_int_tup):
    match_dict={}
    vac_list=get_vac_list(dic)
    prob_list=[]
    group_prob_list=[]
    sum_prob=0
    #print vac_list
    index=0
    for v in vac_list:
        one_vac_probs=get_single_trans_prob(dic,v,sum_prob,dim,T,v_int_tup)
        prob_list+=one_vac_probs
        group_prob_list.append(one_vac_probs)
        sum_prob=one_vac_probs[-1]
        index+=1
    for i in range(0,len(vac_list)):
        for j in range(0,4):
            group_prob_list[i][j]/=sum_prob
    for i in range(0,len(vac_list)):
        match_dict[i]=vac_list[i]
    return match_dict,group_prob_list,sum_prob
    
def choose_trans(match_dict,group_prob_list):
    toss=random.random()
    last_ind=len(group_prob_list)-1
    guess=int(last_ind*toss)
    n=0
    flag=0
    while True:
        lv=group_prob_list[guess][0]
        uv=group_prob_list[guess][-1]
        #print toss,n,guess,lv,uv
        if lv<toss and uv>toss:
            break
        else:
            if lv>toss:
                guess-=1
                if guess<0:
                    guess=0
                    i=0
                    flag=1
                    break
                elif group_prob_list[guess][-1]<toss:
                    guess+=1
                    flag=1
                    i=0
                    break
            else:
                guess+=1
        n+=1
    if flag==0:
        good_list=group_prob_list[guess]
        for i in range(0,len(good_list)):
            #print i,toss,good_list
            if good_list[i]>toss:
                break         
        return match_dict[guess], i
    else:
        return match_dict[guess], i

def get_index_from_tuple(tup,dim):
    i=dim*tup[0]
    j=int(tup[1]/2)
    return i+j

def make_trans(dic,hop_vacancy,i,dim):
    hop_dir=index_move_translation(i)
    hop_site=((hop_vacancy[0]+hop_dir[0])%(2*dim),(hop_vacancy[1]+hop_dir[1])%(2*dim))
    dic[hop_site]
    dic[hop_vacancy]=dic[hop_site]
    dic[hop_site]=vec_V
    #return dic, hop_vacancy, hop_site
    return hop_vacancy, hop_site
    
def plot_conf_dic(dic,dim):    
    a=dic.keys()
    lx=[]
    ly=[]
    val_i=[]
    for i in a:
        lx.append(i[0])
        ly.append(i[1])
        val_i.append(dic[i])
    plt.scatter(lx,ly,c=val_i,s=3000/dim)
    #plt.show()
    
def live_ordering_plot(num_of_steps,d,dim,T,v_int_tup):
    #d=gen_random_lattice(dim,0.5,0.1) 
    n=0
    a=sorted(d.keys())
    #print a
    lx=[]
    ly=[]
    val_i=[]
    for i in a:
        lx.append(i[0])
        ly.append(i[1])
        val_i.append(d[i])
    plt.ion()
    plt.figure()
    plt.scatter(lx,ly,c=val_i,s=4000/dim)
    time_MC=0
    plt.title('t(MC)='+str(time_MC))
    while n<num_of_steps:
        match_dict,group_prob_list,sum_prob=get_trans_probs(d,dim,T,v_int_tup)
        hop_vacancy,i=choose_trans(match_dict,group_prob_list)
        hv,hs=make_trans(d,hop_vacancy,i,dim)
        print hv, hs
        hs_ind=get_index_from_tuple(hs,dim)
        color_site=val_i[hs_ind]
        val_i[get_index_from_tuple(hv,dim)]=color_site
        val_i[hs_ind]=vec_V
        plt.scatter(lx,ly,c=val_i,s=4000/dim)
        time_MC+=(1/sum_prob)
        plt.title('t(MC)='+str(time_MC))
        plt.pause(0.003*time_MC)
        plt.clf()
        n+=1


    