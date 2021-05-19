import numpy as np
import math 
import random
import matplotlib.pyplot as plt
from scipy.stats import poisson
import pandas as pandas
from scipy.stats import uniform
import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad

# Code to calculate Probability of survival for CHO cell at different values of CHI

# Ng is taken constant as Nbp is out of scope 
Ng=2500

# Nr and Lethality Criterion 
N_r=np.random.uniform(0.04,0.08)
Lambda=0.25

# Probability of producing lethal damage
Pl_r =Lambda*(1-poisson.pmf(1,N_r)-poisson.pmf(0,N_r)-poisson.pmf(2,N_r))
sig=0.371

doseCHO = np.arange(0,6,0.0001)

Se6=150   # Stopping potential for CHO
nsCHO=(4.26)*(10**(-4))  # Number density of complex damage sites on chromatin of CHO
NgCHO=((nsCHO*16*127*9.9)/math.pi)*1000  # Genome size of CHO

X0=0.35
X1=0.04

# Equation to find the number of lethal lesions for CHO (Equation 15)
Yl=((math.pi)*sig*Ng*doseCHO)/(16*Se6)

def Yield  (Ng,Se):
    return (math.pi*sig*Ng*doseCHO)/(16*Se)

YlCHO=Yield(NgCHO,Se6)

# CALCULATING THE PROBABILITY OF CELL SURVIVAL TAKING DIFFERENT CHI VALUES
def Different_chi10(a,b,YL):
    
    Pie=[]
    for i in Yl:
        
        if(i<(a/b)):
            Pie_surv=math.exp(-(1-a)*i - b*i*i)
            Pie.append(Pie_surv)
        else:
            Pie_surv=math.exp(-i)
            Pie.append(Pie_surv)
            
    return Pie

def chi_different(YlCHO,kai):
    
    piesur=[]
    for i in YlCHO:
        
        piesur.append(math.exp(-1*(1-kai)*(i)))
        
    return piesur

# GRAPHS
# Graph 1
# Plotting probability of survival at different CHI values
plt.title("Varying CHI")
plt.yscale("log")
a=chi_different(YlCHO,0.9)
b=chi_different(YlCHO,0.7)
c=chi_different(YlCHO,0.5)
d=chi_different(YlCHO,0.3)
e=chi_different(YlCHO,0.1)
plt.plot(doseCHO,a,color ='r',label ='X(chi)=0.9')
plt.plot(doseCHO,b,label ='X(chi)=0.7')
plt.plot(doseCHO,c,label ='X(chi)=0.5')
plt.plot(doseCHO,d,label ='X(chi)=0.3')
plt.plot(doseCHO,e,label ='X(chi)=0.1')
plt.xlabel('dose')
plt.ylabel('Probability of cell survival ')
plt.legend()
plt.show()

# Graph 2
plt.title("Varying CHI1 with CHI0 constant")
plt.yscale("log")
Pie1=Different_chi10(0.5,0.01,YlCHO)
Pie2=Different_chi10(0.5,0.11,YlCHO)
Pie3=Different_chi10(0.5,0.23,YlCHO)
Pie4=Different_chi10(0.5,0.45,YlCHO)
Pie5=Different_chi10(0.5,0.5,YlCHO)
plt.plot(doseCHO,Pie1,color ='r',label ='x1<<x0   | 0.04')
plt.plot(doseCHO,Pie2,color ='g',label ='x1<<x0   | 0.11')
plt.plot(doseCHO,Pie3,color ='b',label ='x<x0        | 0.23')
plt.plot(doseCHO,Pie4,color ='purple',label ='x1<x0      | 0.45')
plt.plot(doseCHO,Pie5,color ='orange',label ='x1=x0      | 0.5')
# plt.plot(doseCHO,Pie3)
plt.xlabel('dose')
plt.ylabel('Probability of cell survival ')
plt.legend()
plt.show()

# Graph 3
plt.title("Varying CHI0 with CHI1 constant")
plt.yscale("log")
Pie1=Different_chi10(0.05,0.05,YlCHO)
Pie2=Different_chi10(0.1,0.05,YlCHO)
Pie3=Different_chi10(0.5,0.05,YlCHO)
Pie4=Different_chi10(0.6,0.05,YlCHO)
Pie5=Different_chi10(0.9,0.05,YlCHO)
plt.plot(doseCHO,Pie1,color ='r',label ='x1=x0   | 0.04')
plt.plot(doseCHO,Pie2,color ='g',label ='x1<x0     | 0.1')
plt.plot(doseCHO,Pie3,color ='b',label ='x1<x0     | 0.5')
plt.plot(doseCHO,Pie4,color ='purple',label ='x1<<x0  | 0.6')
plt.plot(doseCHO,Pie5,color ='orange',label ='x1<<x0  | 0.9')
plt.xlabel('dose')
plt.ylabel('Probability of cell survival ')
# plt.plot(doseCHO,Pie3)
plt.legend()
plt.show()

# Graph 4
plt.title("Varying CHI1 and CHI0 keeping ratio of CHI0/CHI1 constant")
plt.yscale("log")
Pie1=Different_chi10(0.2,0.01,YlCHO)
Pie2=Different_chi10(0.4,0.02,YlCHO)
Pie3=Different_chi10(0.8,0.04,YlCHO)
plt.plot(doseCHO,Pie1,color ='r',label ='X0=0.2 | X1=0.01')
plt.plot(doseCHO,Pie2,color ='g',label ='X0=0.4 | X1=0.02')
plt.plot(doseCHO,Pie3,color ='b',label ='X0=0.8 | X1=0.04')
plt.xlabel('dose')
plt.ylabel('Probability of cell survival ')
plt.legend()
plt.show()
