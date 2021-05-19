import numpy as np
import math 
import random
import matplotlib.pyplot as plt
from scipy.stats import poisson
import pandas as pandas
from scipy.stats import uniform

# Two small innovations are performed in this code

#Small Innovations 1 - Get the Dose after which there will be no Cell Repair
# Select the Cell Line and Input the LET

# 1 

pi=3.142

X0=0.35
X1=0.04

sig=0.373

# ns is the number density of complex damage sites on chromatin of the respective cell
nsA549=(1.2)*(10**(-3))
nsAG1522=(4.2)*(10**(-4))
nsHela=(2.2)*(10**(-4))
nsNB1RGB=(3.2)*(10**(-4))
nsA1722=(2.4)*(10**(-4))
nsV79=(7.2)*(10**(-4))
nsCHO=(4.26)*(10**(-4))

# Stopping potential of the respective cells
Se1= 100 # A549
Se2=122# Ag1522
Se3=70 # HelA
Se4=54 #NB1RG
Se5=105 # A1722
Se6=120 #V79
Se7=150 # CHO

doseCHO=np.arange(0,6,0.0001)
doseCommon=np.arange(0,4,0.0001)

# Radius of respective cells
rA549=4.8
rAG1522=6.7
rHela=8.35
rNB1RGB=7.4
rA1722=8.15
rV79=5.3
rCHO=6.36

# Genome size of respective cells
NgA549=((nsA549*16*100*7.5)/math.pi)*1000
NgAG1522=((nsAG1522*16*144*10.6)/math.pi)*1000
NgHela=((nsHela*16*219*13)/math.pi)*1000
NgNB1RGB=((nsNB1RGB*16*172*11.6)/math.pi)*1000
NgA1722=(nsA1722*16*209*12.7)/math.pi
NgV79=((nsV79*16*88*8.2)/math.pi)*1000
NgCHO=((nsCHO*16*127*9.9)/math.pi)*1000

Seip=int(input("Input the Energy and it will show you the Dose at Maximal Dose at which there can be any cell reapir possible"))

CellSelect=int(input("A549:1, Ag1522:2, Hela :3, NB1RG:4, A1722:5,V79:5,CHO:7 "))

dosetest=np.arange(0,11,0.001)

# Function to find the number of lethal lesions i.e. Yl for the selected cell
def getYield(Ngi,Sei,sigi):
    return (pi*sigi*Ngi*dosetest)/(16*Sei)

def desiredYl(CellSelect):
    if (CellSelect==1):
        a=getYield(NgA549,Seip,sig)
    elif (CellSelect==2):
        a=getYield(NgAG1522,Seip,sig)
    elif (CellSelect==3):
        a=getYield(NgHela,Seip,sig)
    elif (CellSelect==4):
        a=getYield(NgNB1RGB,Seip,sig)
    elif (CellSelect==5):
        a=getYield(NgA1722,Seip,sig)
    elif (CellSelect==6):
        a=getYield(NgV79,Seip,sig)
    else:
        a=getYield(NgCHO,Seip,sig)
    return a


b=desiredYl(CellSelect)

# print(max(b)) # To get the max Yl

# "val" is the maximum value of dose to complete the cell repair
if (CellSelect==1):
    for i in b:
        if (i<X0/X1):
            pass
            # Still cell reapir is going on
        else:
            val=(i/((pi*0.371*NgA549)/(16*Seip))) 
            print(val)
            break
        
elif  (CellSelect==2):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir is going on
        else:
            val=(i/((pi*0.371*NgAG1522)/(16*Seip)))
            print(val)
            break
        
elif (CellSelect==3):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir is going on
        else:
            val=(i/((pi*0.371*NgHela)/(16*Seip)))
            print(val)
            break
        
elif (CellSelect==4):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir is going on
        else:
            val=(i/((pi*0.371*NgNB1RGB)/(16*Seip)))
            print(val)
            break
        
elif (CellSelect==4):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir is going on
        else:
            val=(i/((pi*0.371*NgA1722)/(16*Seip)))
            print(val)
            break

elif (CellSelect==6):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir is going on
        else:
            val=(i/((pi*0.371*NgV79)/(16*Seip)))
            print(val)
            break
        
else:
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir is going on
        else:
            val=(i/((pi*0.371*NgCHO)/(16*Seip)))
            print(val)
            break


# 2. Small Innovation of Finding the Number of DSB induced Lesions 


print(" For CHO Cells " )
Narray=[]

for i in doseCHO:
    N=i*pi*(sig*NgCHO*i)/(16*Se1) # Dose for CHO
    Narray.append(N)

plt.plot(doseCHO,Narray,label='The number of Lesions per Cell')
plt.xlabel('dose')
plt.ylabel('number of Lesions per Cell ')
plt.legend()
plt.show()


inter=pi*(rCHO**2)*15*10000000/(Se7*1.602)
print(inter)

larr=[] # No. of losions per particle
narr=[] # No. of particles
### Lambda ###

# For calculating lesions per particle
for i in doseCHO:
    N=i*pi*(sig*NgCHO*i)/(16*Se1)
    n=i*inter ### Number of Particles
    narr.append(n)
    l=N/n  ### Lesions Per Particle 
    larr.append(l)

plt.plot(doseCHO,larr,label='The DSBs per Particle Irradiated')
plt.xlabel('dose')
plt.ylabel('number of Lesions per Cell ')
plt.legend()
plt.show()

print("The Maximum  number of Particles are")
print(max(narr))











