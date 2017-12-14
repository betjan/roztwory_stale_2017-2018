# -*- coding: utf-8 -*-
import os

os.chdir("C:/Users/student/Desktop")

# Mój pierwszy program

#print "Hello world"

#z=range(0,10)

#print z

#for i in z:
    #print i

#To jest lista
at_list=['Ni','Al']
#To jest krotka
position=(0,0,0)

dim_x=3
dim_y=3
dim_z=3

Ni_list=[]
Al_list=[]

for i in range(0,dim_x):
    for j in range(0,dim_y):
        for k in range(0,dim_z):
            Ni_list.append((i,j,k))
            Al_list.append((i+0.5,j+0.5,k+0.5))

num_of_atoms=dim_x*dim_y*dim_z*2

with open('B2.xyz','w') as w:
    w.write(str(num_of_atoms)+'\n\n')
    for ni in Ni_list:
        w.write("%s %.3f %.3f %.3f\n"%('Ni',ni[0],ni[1],ni[2]))
    for al in Al_list:
        w.write("%s %.3f %.3f %.3f\n"%('Al',al[0],al[1],al[2]))
