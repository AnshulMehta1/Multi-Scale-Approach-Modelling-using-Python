import numpy as np
import math 
import random
import matplotlib.pyplot as plt
from scipy.stats import poisson
import pandas as pandas
from scipy.stats import uniform

### Small Innovations - Get the Dose after which there will be no Cell Repair
### Select the Cell Line and Input the LET

pi=3.142


X0=0.35
X1=0.04

sig=0.373
nsA549=(1.2)*(10**(-3))
nsAG1522=(4.2)*(10**(-4))
nsHela=(2.2)*(10**(-4))
nsNB1RGB=(3.2)*(10**(-4))
nsA1722=(2.4)*(10**(-4))
nsV79=(7.2)*(10**(-4))
nsCHO=(4.26)*(10**(-4))

Se1= 100 # A549
Se2=122# Ag1522
Se3=70 # HelA
Se4=54 #NB1RG
Se5=105 # A1722
Se6=120 #V79
Se7=150 # CHO

doseCHO=np.arange(0,6,0.0001)
doseCommon=np.arange(0,4,0.0001)

rA549=4.8
rAG1522=6.7
rHela=8.35
rNB1RGB=7.4
rA1722=8.15
rV79=5.3
rCHO=6.36

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

print(max(b))
 
if (CellSelect==1):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir gooing on
        else:
            val=(i/((pi*0.371*NgA549)/(16*Seip)))
            print(val)
            break
elif  (CellSelect==2):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir gooing on
        else:
            val=(i/((pi*0.371*NgAG1522)/(16*Seip)))
            print(val)
            break
elif (CellSelect==3):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir gooing on
        else:
            val=(i/((pi*0.371*NgHela)/(16*Seip)))
            print(val)
            break
elif (CellSelect==4):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir gooing on
        else:
            val=(i/((pi*0.371*NgNB1RGB)/(16*Seip)))
            print(val)
            break
elif (CellSelect==4):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir gooing on
        else:
            val=(i/((pi*0.371*NgA1722)/(16*Seip)))
            print(val)
            break

elif (CellSelect==6):
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir gooing on
        else:
            val=(i/((pi*0.371*NgV79)/(16*Seip)))
            print(val)
            break
else:
    for i in b:
        if (i<X0/X1):
            pass
            #  Still cell reapir gooing on
        else:
            val=(i/((pi*0.371*NgCHO)/(16*Seip)))
            print(val)
            break


### Small Innovation of Finding the Number of DSB induced Lesions ###


print(" For CHO Cells " )
Narray=[]

for i in doseCHO:
    N=i*pi*(sig*NgCHO*i)/(16*Se1)
    Narray.append(N)

plt.plot(doseCHO,Narray,label='The number of Lesions per Cell')
plt.show()

inter=pi*(rCHO**2)*15*10000000/(Se7*1.602)
print(inter)
larr=[]
narr=[]
### Lambda ###


for i in doseCHO:
    N=i*pi*(sig*NgCHO*i)/(16*Se1)
    n=i*inter ### Number of Particles
    narr.append(n)
    l=N/n  ### Lesions Per Particle 
    larr.append(l)

plt.plot(doseCHO,larr,label='The DSBs per Particle Irradiated')
plt.show()
print("The Maximum  number of Particles are")
print(max(narr))











