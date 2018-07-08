# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 12:15:59 2018

@author: kgratsea
"""

# -*- coding: utf-8 -*-
"""
Quantum walk with GDT. 3 free parameters for each coin operator

function to be minimzed norm

@author: kgratsea
"""

import numpy as np
import cmath
import math
from scipy import optimize
global final
import random 

def tensor(vectorA,vectorB) :
    m = np.size(vectorA,0)
    n = np.size(vectorB,0)
    tens=np.zeros((m,n))
    for i in range(m) :
        for j in range(n) :
            tens[i][j] = vectorA[i]*vectorB[j]
    return (tens);

def func(z) :    
    n=5 #number of steps
    k=n+1 #number of sites at the final state
    
    initial = np.zeros((2*k,1),dtype=complex)
    #localised on one site
    initial[0][0]= 1.
    initial[1][0]= 1.5
    initial/= np.linalg.norm(initial)
    
    Initial = initial
    #print (Initial)
    
    f = open("test_GDT.txt","a+")
    f.write("Initial")
    f.close()
    with open('test_GDT.txt', 'a+') as f:
       print (Initial,file=f)
    f.close()
   
    
    #definition of invS
    invS = np.zeros((2*k,2*k),dtype=complex)
    matrixS = np.zeros((2*k,2*k),dtype=complex)
    for i in range (0,2*k,2) :
        invS[0+i][0+i] =1.
        matrixS[0+i][0+i] =  1.
        if (i+3)< 2*k :
            invS[1+i][3+i] = 1. #S-1
            matrixS[3+i][1+i] = 1.
    
    listSt = []
    listc = []
    listC = []

    listSt.append (initial)
    
    #Define coin operators with gdt
    
    l = 0 # for corresponding the correct coin parameters at each step n
    for j in range (0,n,+1) : 
        print ("n",j)
        c=np.zeros((2,2),dtype=complex)
        x=abs(z[0+l])
        y=z[1+l]
        v=z[2+l]
        while (1-x<0): 
            print ("here1")
            x /= 10.
            print (x)
        
        c[0][0]=   math.sqrt(x)
        c[0][1]= (math.sqrt(1-x)) * (math.cos(y*math.pi) + math.sin(y*math.pi)*1j) 
        c[1][0]= (math.sqrt(1-x)) * (math.cos(v*math.pi) + math.sin(v*math.pi)*1j)         
        c[1][1]= -(math.sqrt(x)) * (math.cos((y+v)*math.pi) + math.sin((y+v)*math.pi)*1j)  
        
        listc.append(c)
        matrixC = np.zeros((2*k,2*k),dtype=complex)
        print (c)
        
        for i in range (0,2*k,2):
            matrixC[0+i][0+i] = c[0][0]
            matrixC[1+i][1+i] = c[1][1]
            matrixC[0+i][1+i] = c[0][1]          
            matrixC[1+i][0+i] = c[1][0]   
         
        listC.append (matrixC)    
        
        
        m1 = np.dot(matrixC,initial)
        m2 = np.dot(matrixS,m1)   #next state
        print (m2)
        listSt.append (m2)
        initial = m2/np.linalg.norm(m2)
        l += 3 # moving to the next coin parameters
        
    Phi=initial    
    Phi_internal = np.delete(Phi,[0,1,2*k-2,2*k-1],None)
    print ("Phi_internal",Phi_internal)
    
    print ("Phi_internal",Phi_internal)
    Phi_internal_reshaped =np.zeros((2,4))
    for i in range(0,2,1):
        q=0
        for j in range(0,4,1): 
            Phi_internal_reshaped[i][j]= Phi_internal[i+q][0]
            q +=2

    print ("Phi_internal_reshaped",Phi_internal_reshaped)
    
    psiA, l, psiB = np.linalg.svd(Phi_internal_reshaped,full_matrices=1) #decomposition of initial matrix
    print ("l",l)
    
    NORM=0.0
    p=1.0 # p has to be larger or equal than 1 for the algorithm to work
    
    m = np.size(Phi_internal_reshaped,0)  #number of rows of initial matrix
    n = np.size(Phi_internal_reshaped,1)  #number of columns of initial matrix
    sum=np.zeros((m,n))
    
    for i in range(2) :
         tens= tensor(psiA[i],psiB[i])
         sum += math.pow(l[i],p-1)*tens
         NORM = NORM + math.pow(l[i],p) 

    NORM = math.pow(NORM,1./p)

    print (NORM)
    return (-NORM+math.sqrt(2))

    #max.ent.over sites i ={2,3}
    #Phi_target= np.array([[ 0.0066874],       [ 0.       ],       [ 0.6148   ],       [ 0.3492   ],       [ 0.3493   ],       [-0.6148   ],       [ 0.       ],       [ 0.0067874]])
    
    '''
    #max.ent.over sites i ={2,5}    
    #d = [0.48107031, 0.63099229]  #solution of orthogonality condition
    #Phi_target = np.array([ [0.4811], [0.52619533], [0.27592482], [ 0.26093188],  [0.2510644],[0 ], [0], [0.2290105] , [ 0.48391584], [ 0.3214155], [ 0.35487592],[0.6310]])
    Phi_target= np.array([  [0.4811], [0] ,[0.52619533],  [0.2290105] ,  [0.27592482],  [ 0.48391584], [ 0.26093188],  [ 0.3214155], [0.2510644], [ 0.35487592] ,[0 ], [0.6310] ])
    Phi_target /= np.linalg.norm(Phi_target)
    Phi_target_internal = np.delete(Phi_target,[0,1,2*k-2,2*k-1],None)
    print ("Phi_target",Phi_target) '''

   
   
#initial_coin_parameters=[1/2,0,math.pi,1/4,math.pi/2.,math.pi/2.,1/8,2*math.pi,0]  #n=3
initial_coin_parameters=[1/2,0,math.pi,1/4,math.pi/2.,math.pi/2.,1/8,2*math.pi,0,1/2,3*math.pi,math.pi/4,1/3,2*math.pi/3,math.pi/3] #n=5

'''
f=3
for i in range (1,3*f,3):    
    initial_coin_parameters.append(random.uniform(0,1))   
    initial_coin_parameters.append(random.uniform(0,2*math.pi))   
    initial_coin_parameters.append(random.uniform(0,2*math.pi))   
    '''
Initial_coin_par=initial_coin_parameters  
        
minimizer_kwargs = {"method": "BFGS"}
ret = optimize.basinhopping(func,initial_coin_parameters, minimizer_kwargs=minimizer_kwargs,niter=3, T=1.0, disp = True )
  
l=0
f=3
listc=[]
for j in range (0,f,+1) : 
        print ("j",j)
        c=np.zeros((2,2),dtype=complex)
        x=abs(ret.x[0+l])
        while (1-x<0): 
            print ("here1")
            x /= 10.
            print (x)
        y=ret.x[1+l]
        v=ret.x[2+l]
        c[0][0]=   math.sqrt(x)
        c[0][1]= (math.sqrt(1-x)) * (math.cos(y*math.pi) + math.sin(y*math.pi)*1j) 
        c[1][0]= (math.sqrt(1-x)) * (math.cos(v*math.pi) + math.sin(v*math.pi)*1j)         
        c[1][1]= -(math.sqrt(x)) * (math.cos((y+v)*math.pi) + math.sin((y+v)*math.pi)*1j)  
        
        listc.append(c)
        l+=3
