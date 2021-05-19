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

# In this code, the dose is made time dependent and as a result, a significant decrease in relative cell survival is observed. 

# Calculating probability of survival for CHO cells
def Survival_Curves(doseCHO):
    
    Ng=2500
    sig=0.371
    Se6=150
    X0=0.35
    X1=0.04

    Yl=((math.pi)*sig*Ng*doseCHO)/(16*Se6)

    Pie=[]
    
    for i in Yl:
        
        if(i<(X0/X1)):
            Pie_surv=math.exp(-(1-X0)*i - X1*i*i)
            Pie.append(Pie_surv)
        else:
            Pie_surv=math.exp(-i)
            Pie.append(Pie_surv)

    return Pie

# Total time interval = 8mins which is dividen in 100 fractions
time=np.arange(0,8,0.1)

def time_dependent_dose(time,x,y):
    return x*time/y

# No. of doses
doseCHO1 = time_dependent_dose(time,3,3.1)#0.96
doseCHO2 = time_dependent_dose(time,3,3.2)#0.93
doseCHO3 = time_dependent_dose(time,3,3.3)#0.90
doseCHO4 = time_dependent_dose(time,3,3.4)#0.88
doseCHO5 = time_dependent_dose(time,3,3.5)#0.85

# Probability of survival for each dose
Pie=Survival_Curves(doseCHO1)
Pie=Survival_Curves(doseCHO2)
Pie=Survival_Curves(doseCHO3)
Pie=Survival_Curves(doseCHO4)
Pie=Survival_Curves(doseCHO5)

# Plotting the probability of survival
plt.yscale("log")
plt.plot(doseCHO1,Pie,label="Fraction Parameter = 0.96")
plt.plot(doseCHO2,Pie,label="Fraction Parameter = 0.93")
plt.plot(doseCHO3,Pie,label="Fraction Parameter = 0.90")
plt.plot(doseCHO4,Pie,label="Fraction Parameter = 0.88")
plt.plot(doseCHO5,Pie,label="Fraction Parameter = 0.85")
plt.xlabel('dose')
plt.ylabel('Probability of cell survival ')
plt.legend()
plt.show()
