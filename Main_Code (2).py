import numpy as np
import matplotlib.pyplot as plt
from random import randint
import math
import scipy. stats as ss
from math import exp


dose = np.arange(0,8,0.001)
x0 = 0.35
x1 = 0.04
sig = 0.371
#Ng = 2500
Se = 150
lamda=0.15



Dn_A549= 9.6
Dn_AG1522= 13.4
Dn_Hela=16.7
Dn_NB1RGB=14.8
Dn_A1722=16.3
Dn_V79=10.6
Dn_CHO = 12.7

ns_A549=(1.2)*(10**(-3))
ns_AG1522=(4.2)*(10**(-4))
ns_Hela=(2.2)*(10**(-4))
ns_NB1RGB=(3.2)*(10**(-4))
ns_A1722=(2.4)*(10**(-4))
ns_V79=(7.2)*(10**(-4))
ns_CHO=(4.26)*(10**(-4))


rA549=4.8
rAG1522=6.7
rHela=8.35
rNB1RGB=7.4
rA1722=8.15
rV79=5.3
rCHO=6.36

Se1= 100
Se1H=25
Se2=122
Se2H=17
Se3=70
Se3H=15
Se4=54
Se4H=13
Se5=105
Se5H=78
Se6=120
Se6H=100
Se7=150
Se7H=72

#gradNumOfLesion = Ns*sig

def getZ(Dn):
    return (math.pi*Dn)/4
    
z_A549= getZ(Dn_A549)
z_AG1522= getZ(Dn_AG1522)
z_Hela=getZ(Dn_Hela)
z_NB1RGB=getZ(Dn_NB1RGB)
z_A1722=getZ(Dn_A1722)
z_V79=getZ(Dn_V79)
z_CHO = getZ(Dn_CHO)

def getNbp(Ns,An,z):
    return (160*An*z*Ns)/(3*math.pi)

Nbp_A549= getNbp(ns_A549,144,z_A549)
Nbp_AG1522= getNbp(ns_AG1522,100,z_AG1522)
Nbp_Hela=getNbp(ns_Hela,219,z_Hela)
Nbp_NB1RGB=getNbp(ns_NB1RGB,172,z_NB1RGB)
Nbp_A1722=getNbp(ns_A1722,209,z_A1722)
Nbp_V79=getNbp(ns_V79,88,z_V79)
Nbp_CHO = getNbp(ns_CHO,127,z_CHO)


def getNg(Nbp):
    return Nbp/3.33

Ng_A549= getNg(Nbp_A549)*1000
Ng_AG1522= getNg(Nbp_AG1522)*1000
Ng_Hela=getNg(Nbp_Hela)*1000
Ng_NB1RGB=getNg(Nbp_NB1RGB)*1000
Ng_A1722=getNg(Nbp_A1722)*1000
Ng_V79=getNg(Nbp_V79)*1000
Ng_CHO = getNg(Nbp_CHO)*1000

def Yield  (Ng,Se):
    return (math.pi*sig*Ng*dose)/(16*Se)

YlA549 = Yield(Ng_A549,Se1)
YlAG1522 = Yield(Ng_AG1522,Se2)
YlHela = Yield(Ng_Hela,Se3)
YlNB1RGB = Yield(Ng_NB1RGB,Se4)
YlA1722 = Yield(Ng_A1722,Se5)
YlV79 = Yield(Ng_V79,Se6)
YlCHO = Yield(Ng_CHO,Se7)


# N_r=randint(1,10)
N_r=np.random.uniform(0.04,0.08)


Pl_r= (1- (ss.poisson.pmf(0, N_r) + ss.poisson.pmf(1, N_r) + ss.poisson.pmf(2, N_r)))
Pl_r=Pl_r*lamda

#actual cell survival curve
def expected_survival(n):
    cell_surv1=[]
    for i in n:
        cell_surv=math.exp(-i)
        cell_surv1.append(cell_surv)
    return cell_surv1

Pie=[]
for i in YlCHO:
    if(i<(x0/x1)):
        Pie_surv=math.exp(-(1-x0)*i - x1*i*i)
        Pie.append(Pie_surv)
    else:
        Pie_surv=math.exp(-i)
        Pie.append(Pie_surv)



Pie_expected = expected_survival(YlCHO)
plt.yscale('log')
plt.title('Transition from Linear Quadratic to Linear regime')
plt.plot(dose,Pie,color='orange',label='  Cell survival curve ')
plt.legend()
plt.show()

plt.yscale('log')
plt.plot(dose,Pie,label=' actual Cell survival curve with repair ')
plt.plot(dose,Pie_expected,color ='g',label ='expected survival curve in absence of repair factor of')
plt.xlabel('dose')
plt.ylabel('Probability of cell survival ')
plt.legend()
plt.show()


damage=[]
for i in Pie:
    val=1-i
    damage.append(val)


# plt.xlabel('dose')
# plt.ylabel('Probability of damage')
# plt.plot(dose,damage)
# plt.show()

dose=np.arange(0,2,0.0001)

def MSA_to_LQ (n,a,b):
    d=[]
    for i in n:
        d.append(exp(-a*i  - b*pow(i,2)))
    return d


def alpha(X0,X1,Ng,Se):
    return(((1-X0)*math.pi*sig*Ng/(16*Se)))

def beta(X0,X1,Ng,Se):
    return(X1*pow((math.pi*sig*Ng)/(16*Se),2))

Alpha_A549=alpha(x0,x1,Ng_A549,Se1)
Alpha_AG1522=alpha(x0,x1,Ng_AG1522,Se2)
Alpha_Hela=alpha(x0,x1,Ng_Hela,Se3)
Alpha_NB1RGB=alpha(x0,x1,Ng_NB1RGB,Se4)
Alpha_A1722=alpha(x0,x1,Ng_A1722,Se5)
Alpha_V79=alpha(x0,x1,Ng_V79,Se6)
Alpha_CHO=alpha(x0,x1,Ng_CHO,Se7)

Beta_A549=beta(x0,x1,Ng_A549,Se1)
Beta_AG1522=beta(x0,x1,Ng_AG1522,Se2)
Beta_Hela=beta(x0,x1,Ng_Hela,Se3)
Beta_NB1RGB=beta(x0,x1,Ng_NB1RGB,Se4)
Beta_A1722=beta(x0,x1,Ng_A1722,Se5)
Beta_V79=beta(x0,x1,Ng_V79,Se6)
Beta_CHO=beta(x0,x1,Ng_CHO,Se7)


plt.yscale("log")
plt.xlabel('Dose (Gy)')
plt.ylabel('Survival Fraction(pi)')
plt.plot(dose,MSA_to_LQ(dose,Alpha_A549,Beta_A549),label="LQ equation of A546 Cell line")
plt.plot(dose,MSA_to_LQ(dose,Alpha_AG1522,Beta_AG1522),label="LQ equation of AG1522 Cell line")
plt.plot(dose,MSA_to_LQ(dose,Alpha_Hela,Beta_Hela),label="LQ equation of HeLa Cell line")
plt.plot(dose,MSA_to_LQ(dose,Alpha_NB1RGB,Beta_NB1RGB),label="LQ equation of NB1RGB Cell line")
plt.plot(dose,MSA_to_LQ(dose,Alpha_A1722,Beta_A1722),label="LQ equation of A172 Cell line")
plt.plot(dose,MSA_to_LQ(dose,Alpha_V79,Beta_V79),label="LQ equation of V79 Cell line")
plt.plot(dose,MSA_to_LQ(dose,Alpha_CHO,Beta_CHO),label="LQ equation of CHO Cell line")

plt.legend()
plt.show()