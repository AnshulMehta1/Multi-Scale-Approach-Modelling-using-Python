import numpy as np
import math
import random
import matplotlib.pyplot as plt
from scipy.stats import poisson
import pandas as pandas
from scipy.stats import uniform

rA549=4.8
rAG1522=6.7
rHela=8.35
rNB1RGB=7.4
rA1722=8.15
rV79=5.3
rCHO=6.36

X0=0.35
X1=0.04

dose = np.arange(0,6,0.0001)
doseCHO= np.arange(0,6,0.0001)
doseCommon=np.arange(0,4,0.0001)

pi=3.14
Ng=2500
## Stopping Potential
Se=150
sig=0.371
Yl=(pi*sig*Ng*dose)/(16*Se)

Se1= 100
Se1H=25
Se2=122
Se2H=17
Se3=70 # HelA
Se3H=15
Se4=54 #NB1RG
Se4H=13
Se5=105 # A1722
Se5H=78
Se6=120 #V79
Se6H=100
Se7=150 # CHO
Se7H=72

nsA549=(1.2)*(10**(-3))
nsAG1522=(4.2)*(10**(-4))
nsHela=(2.2)*(10**(-4))
nsNB1RGB=(3.2)*(10**(-4))
nsA1722=(2.4)*(10**(-4))
nsV79=(7.2)*(10**(-4))
nsCHO=(4.26)*(10**(-4))


NgA549=((nsA549*16*100*7.5)/math.pi)*1000
NgAG1522=((nsAG1522*16*144*10.6)/math.pi)*1000
NgHela=((nsHela*16*219*13)/math.pi)*1000
NgNB1RGB=((nsNB1RGB*16*172*11.6)/math.pi)*1000
NgA1722=(nsA1722*16*209*12.7)/math.pi
NgV79=((nsV79*16*88*8.2)/math.pi)*1000
NgCHO=((nsCHO*16*127*9.9)/math.pi)*1000


# Linear Alpha
#########

### A549 ###
lalphaA549 = ( 3.14 * 0.371 * NgA549)/ (16.0 * Se1)

lalphahypoxicA549 = ( 3.14 * 0.0371* NgA549)/ (16.0 * Se1H)

### AG1522 ###
lalphaAG1522 = (3.14* 0.371*NgAG1522)/ (16.0 * Se2)

lalphahypoxicAG1522 = (3.14* 0.0371*NgAG1522)/ (16.0 * Se2H)


### Hela ###

lalphaHela = ( 3.14 *0.371* NgHela)/ (16.0 *Se3)


lalphahypoxicHela = (3.14 *0.0371* NgHela)/ (16.0 * Se3H)

### NB1RG ###
lalphaNB1RG = ( 3.14 * 0.371*NgNB1RGB)/ (16.0 * Se4)

lalphahypoxicNB1RG = (3.14 * 0.0371* NgNB1RGB)/ (16.0 * Se4H)


### A1722 ###
lalphaA1722 = (3.14 * 0.371* NgA1722)/ (16.0 * Se5)

lalphahypoxicA1722 = (3.14* 0.0371*Ng)/ (16.0 * Se5H)


### V79 ####
lalphaV79 = ( 3.14* 0.371* NgV79)/ (16.0 * Se6)

lalphahypoxicV79 = (3.14 * 0.0371*NgV79)/ (16.0 * Se6H)


### CHO CeLLAs###

lalphaCHO = ( 3.14 *0.371* NgCHO)/ (16.0 *Se7)


lalphahypoxicCHO = (3.14 *0.0371* NgCHO)/ (16.0 * Se7H)



####  Cells-CHO ###
def func (n):
    return math.exp(-lalphaCHO*n)
y = list(map(func,doseCHO))
def func (n):
    return math.exp(-lalphahypoxicCHO*n)
c = list(map(func,doseCHO))

plt.yscale('log')
p = plt.plot(dose,y,label='High LEt')
q=plt.plot(doseCHO,c,label='Low LET')
plt.legend()
plt.xlabel('Dose (Gy)')
plt.ylabel('Survival Fraction (pi)')
plt.title('Figure 1 CHO')
plt.show()


####  Cells-CHO ###
def func (n):
    return math.exp(-lalphaHela*n)
y = list(map(func,doseCommon))
def func (n):
    return math.exp(-lalphahypoxicHela*n)
c = list(map(func,doseCommon))

plt.yscale('log')
p = plt.plot(doseCommon,y,label='High LEt')
q=plt.plot(doseCommon,c,label='Low LET')
plt.legend()
plt.xlabel('Dose (Gy)')
plt.ylabel('Survival Fraction (pi)')
plt.title('Figure 1 Hela')
plt.show()


#### A549 Cells ###
def func (n):
    return math.exp(-lalphaA549*n)
y = list(map(func,doseCommon))
def func (n):
    return math.exp(-lalphahypoxicA549*n)
c = list(map(func,doseCommon))

plt.yscale('log')
p = plt.plot(doseCommon,y,label='High LEt')
q=plt.plot(doseCommon,c,label='Low LET')
plt.legend()
plt.xlabel('Dose (Gy)')
plt.ylabel('Survival Fraction (pi)')
plt.title('Figure 1 A549 ')
plt.show()



#### AG1522 Cells ###
def func (n):
    return math.exp(-lalphaAG1522*n)
y = list(map(func,doseCommon))
def func (n):
    return math.exp(-lalphahypoxicAG1522*n)
c = list(map(func,doseCommon))

plt.yscale('log')
p = plt.plot(doseCommon,y,label='High LEt')
q=plt.plot(doseCommon,c,label='Low LET')
plt.legend()
plt.xlabel('Dose (Gy)')
plt.ylabel('Survival Fraction (pi)')
plt.title('Figure 1 AG1522')
plt.show()

#### NB1RGB Cells ###
def func (n):
    return math.exp(-lalphaNB1RG*n)
y = list(map(func,doseCommon))
def func (n):
    return math.exp(-lalphahypoxicNB1RG*n)
c = list(map(func,doseCommon))

plt.yscale('log')
p = plt.plot(doseCommon,y,label='High LEt')
q=plt.plot(doseCommon,c,label='Low LET')
plt.legend()
plt.xlabel('Dose (Gy)')
plt.ylabel('Survival Fraction (pi)')
plt.title('Figure 1 NB1RGB')
plt.show()

#### A17222 Cells ###
def func (n):
    return math.exp(-lalphaA1722*n)
y = list(map(func,doseCommon))
def func (n):
    return math.exp(-lalphahypoxicA1722*n)
c = list(map(func,doseCommon))

# plt.yscale('log')
# p = plt.plot(doseCommon,y,label='High LEt')
# q=plt.plot(doseCommon,c,label='Low LET')
# plt.legend()
# plt.xlabel('Dose (Gy)')
# plt.ylabel('Survival Fraction (pi)')
# plt.title('Figure-1 A1722')
# plt.show()

#### For V79 Cells ###
def func (n):
    return math.exp(-lalphaV79*n)
y = list(map(func,doseCommon))
def func (n):
    return math.exp(-lalphahypoxicV79*n)
c = list(map(func,doseCommon))

plt.yscale('log')
p = plt.plot(doseCommon,y,label='High LEt')
q=plt.plot(doseCommon,c,label='Low LET')
plt.legend()
plt.xlabel('Dose (Gy)')
plt.ylabel('Survival Fraction (pi)')
plt.title('Figure 1 V79')
plt.show()
