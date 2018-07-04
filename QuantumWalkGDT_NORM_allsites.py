
"""
Created on Fri Jun 29 12:15:59 2018
@author: kgratsea
"""

"""
Quantum walk with GDT. 3 free parameters for each coin operator
function to be minimzed norm
shift operator move to the left or to the right
@author: kgratsea
"""

import numpy as np
import cmath
import math
from scipy import optimize
global final

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
    Phi_reshaped =np.zeros((2,k),dtype=complex)
    for i in range(0,2,1):
        q=0
        for j in range(0,k,1): 
            Phi_reshaped[i][j] = Phi[i+q][0]
            q +=2
    #Phi_internal = np.delete(Phi,[0,1,2*k-2,2*k-1],None) #delete outer sites
    print ("Phi_reshaped",Phi_reshaped)
    
    psiA, l, psiB = np.linalg.svd(Phi_reshaped,full_matrices=1) #decomposition of initial matrix
    print ("l",l)
    
    NORM=0.0
    p=1.0 # p has to be larger or equal than 1 for the algorithm to work
    
    m = np.size(Phi_reshaped,0)  #number of rows of initial matrix
    n = np.size(Phi_reshaped,1)  #number of columns of initial matrix
    sum=np.zeros((m,n))
    
    for i in range(2) :
         tens= tensor(psiA[i],psiB[i])
         sum += math.pow(l[i],p-1)*tens
         NORM = NORM + math.pow(l[i],p) 

    NORM = math.pow(NORM,1./p)

    print (NORM)
    if ((-NORM+math.sqrt(2))<0.0001) : 
        f = open("GDT_NORM_allsites.txt","a+")
        f.write("Phi_reshaped")
        f.close()
        with open('GDT_NORM_allsites.txt', 'a+') as f:
           print (Phi_reshaped,file=f)
        f.close()
        f = open("GDT_NORM_allsites.txt","a+")
        f.write("l,NORM")
        f.close()
        with open('GDT_NORM_allsites.txt', 'a+') as f:
           print (l,file=f)
           print (NORM,file=f)
        f.close()
        f = open("GDT_NORM_allsites.txt","a+")
        f.write("coin parameters")
        f.close()
        with open('GDT_NORM_allsites.txt', 'a+') as f:
           print (z,file=f)
        f.close()
    return (-NORM+math.sqrt(2))
 
   
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
ret = optimize.basinhopping(func,initial_coin_parameters, minimizer_kwargs=minimizer_kwargs,niter=1, T=1.0, disp = True )
  
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