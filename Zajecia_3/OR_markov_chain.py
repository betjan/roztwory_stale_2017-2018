import random
import numpy as np
import pylab

prob_o=0.5

def test_seq(seq,lista):
    length_seq=len(seq)
    first_num=-length_seq
    flag=0
    if len(lista)>=len(seq):
        for s in seq:
            #print first_num
            #print lista[first_num], s
            if lista[first_num]==s:
                flag=1
            else:
                flag=0
                break
            first_num+=1 
    if flag==1:
         return True
    else:
         return False   
         
#print test_seq('ABB','AAABBBABB')



res_list=[]

coin_list=['A','B']

#print coin_list[-1]
seq1='AB'
seq2='BB'

#for s in seq2:
    #print s

win_1=0
win_2=0

num_of_tries=50000
n=0

mu1_list=[]
mu2_list=[]
'''
while n<num_of_tries:
    while True:
        res=random.choice(coin_list)
        res_list.append(res)
        if test_seq(seq1,res_list):
            win_1+=1
            mu1_list.append(len(res_list))
            break
        if test_seq(seq2,res_list):
            win_2+=1
            mu2_list.append(len(res_list))
            break
    res_list=[]
    n+=1

    
'''

#pylab.hist(mu1_list,bins=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5])
#pylab.hist(mu2_list,bins=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5])
#pylab.show()
'''
while n<num_of_tries:
    while True:
        res=random.choice(coin_list)
        res_list.append(res)
        if test_seq(seq1,res_list):
            win_1+=1
            mu1_list.append(len(res_list))
            break
    res_list=[]
    n+=1

n=0

while n<num_of_tries:
    while True:
        res=random.choice(coin_list)
        res_list.append(res)
        if test_seq(seq2,res_list):
            win_2+=1
            mu2_list.append(len(res_list))
            break
    res_list=[]
    n+=1
'''

pylab.figure()
pylab.hist(mu1_list,bins=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5,7.5,8.5,9.5,10.5,11.5])
pylab.figure()
pylab.hist(mu2_list,bins=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5,7.5,8.5,9.5,10.5,11.5])
pylab.show()

print win_1, win_2
print mu1_list,mu2_list
print np.mean(mu1_list), np.mean(mu2_list)
