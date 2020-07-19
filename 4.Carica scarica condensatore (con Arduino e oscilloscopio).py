import pylab
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mplt

###CARICA

print('\n\n\nCARICA\n\n\n')

x,y=pylab.loadtxt(r'C:\Users\Satellite\Desktop\Università 2\Laboratorio 2\_ESAME\lab2 Giulio\esercitazione 4\dati\DatiCond1_C.txt',unpack=True) #file Arduino:prima colonna tempo,seconda voltaggio
# ATTENZIONE ALLE UNITA' DI MIAURA
dx=4/1000
dy=np.array([1]*len(y))
x=x/1000

#Funzione carica
def carica(x, q, t):
    return q*(1-np.exp(-x/t))

init=(1023,15) #q è il max valore in digit di delta V (1023), t è tau=RC
#UNITA' DI MISURAAAAA

pylab.figure(0) #Serve per vedere se gli init sono giusti
pylab.errorbar(x,y,dy,dx,linestyle = '', color = 'black', marker = '.')
xx=np.linspace(min(x),max(x),2000)
pylab.plot(xx,carica(xx,*init), color='red')

#FIT

sigma=dy
w=1/sigma**2
pars,covm=curve_fit(carica,x,y,init,sigma)
chi2 = ((w*(y-carica(x,*pars))**2)).sum()
ndof=len(x)-len(init)

print('pars:',pars)
print('covm:\n',covm)
dq, dt=np.sqrt(covm.diagonal())
print('\ndq=%f , dt=%f' %(dq, dt))
print ('chi2, ndof:',chi2, ndof)
covnorm=np.zeros((2,2))
#covarianza normalizzata
for i in range(2):
    for j in range(2):
        covnorm[i][j]=covm[i][j]/(np.sqrt(covm[i][i]*covm[j][j]))
print('covnorm:\n', covnorm)



#FIT CON ERRORE EFFICACE

print('\n\nCARICA CON ERROE EFFICACE\n\n')

Deff=np.sqrt(dy**2+((pars[0]/pars[1])*dx*np.exp(-x/pars[1])**2))

#errore efficace
for i in range (2):
    pars,covm=curve_fit(carica,x,y,pars,Deff)
    Deff=np.sqrt(dy**2+((pars[0]/pars[1])*dx*np.exp(-x/pars[1])**2))

sigma=Deff
w=1/sigma**2
parse,covme=curve_fit(carica,x,y,init,sigma)
chi2 = ((w*(y-carica(x,*parse))**2)).sum()
ndof=len(x)-len(init)

print('pars:',parse)
print('covm:\n',covme)
dq, dt=np.sqrt(covme.diagonal())
print('\nda=%f , db=%f' %(dq, dt))
print ('chi2, ndof:',chi2, ndof)
covnorm=np.zeros((2,2))
#covarianza normalizzata
for i in range(2):
    for j in range(2):
        covnorm[i][j]=covme[i][j]/(np.sqrt(covme[i][i]*covme[j][j]))
print('covnorm:\n', covnorm)

#GRAFICO CARICA
pylab.figure(1)
pylab.suptitle('Carica e Scarica condensatore')

pylab.subplot(2,3,1)
pylab.errorbar(x,y,dy,dx,linestyle = '', color = 'black', marker = '.')
pylab.rc('font',size=18)
pylab.ylabel('X  [digit]')
pylab.minorticks_on()
pylab.plot(xx,carica(xx,*parse), color='red')
#Residui
pylab.subplot(2,3,4)
r = (y-carica(x,*parse))/sigma
# bellurie
pylab.rc('font',size=18)
pylab.ylabel('Norm. res.')
pylab.xlabel('t [ms]')
pylab.minorticks_on()
#pylab.ylim((-.9,.9))
pylab.plot(x,r,linestyle="--",color='blue',marker='o')



##SCARICA

x,y=pylab.loadtxt(r'C:\Users\Satellite\Desktop\Università 2\Laboratorio 2\_ESAME\lab2 Giulio\esercitazione 4\dati\DatiCond1_S.txt',unpack=True)

dx=4/1000
dy=np.array([1]*len(y))
x=x/1000

#Funzione scarica due parametri
def scarica2(x, q, t):
    return q*np.e**(-x/t)

#Funzione carica tre parametri
def scarica3(x, q, t, c):
    return q*np.e**(-x/t)+c


init2=(1023,15) #Gli init sono quelli dell carica
init3=(1023, 15, 0.001)

#Guardiamo se gli init sono giusti
pylab.figure(2)
pylab.errorbar(x,y,dy,dx,linestyle = '', color = 'black', marker = '.')
xx=np.linspace(min(x),max(x),2000)
pylab.plot(xx,scarica2(xx,*init2), color='blue')

print('\n\n\nSCARICA 2 PARAMETRI\n\n\n')

#FIT 2

sigma=dy
w=1/sigma**2
pars2,covm2=curve_fit(scarica2,x,y,init2,sigma)
chi2 = ((w*(y-scarica2(x,*pars2))**2)).sum()
ndof=len(x)-len(init2)

print('\n Senza errore efficace \n')
print('pars2:',pars2)
print('covm2:\n',covm2)
dq, dt=np.sqrt(covm2.diagonal())
print('\ndq=%f , dt=%f' %(dq, dt))
print ('chi22, ndof:',chi2, ndof)

#FIT CON ERRORE EFFICACE

print('\n\nERRORE EFFICACE 2 scarica\n\n')
Deff=np.sqrt(sigma**2+((pars2[0]/pars2[1])*dx*np.exp(-x/pars2[1])**2))

for i in range(2):
    pars2,covm2=curve_fit(carica,x,y,pars2,Deff)
    Deff=np.sqrt(sigma**2+((pars2[0]/pars2[1])*dx*np.exp(-x/pars2[1])**2))

sigma=Deff
w=1/sigma**2
pars2,covm2=curve_fit(scarica2,x,y,init2,sigma)
chi2 =((w*(y-scarica2(x,*pars2))**2)).sum()
ndof=len(x)-len(init2)

print('pars:',pars2)
print('covm:\n',covm2)
dq, dr =np.sqrt(covm2.diagonal())
print('\ndq=%f , dr=%f' %(dq, dr))
print ('chi2, ndof:',chi2, ndof)

covnorm=np.zeros((2,2))
#covarianza normalizzata
for i in range(2):
    for j in range(2):
        covnorm[i][j]=covm1[i][j]/(np.sqrt(covm1[i][i]*covm1[j][j]))
print('covnorm:\n', covnorm)


pylab.figure(1)
#Plot
pylab.subplot(2,3,2)
pylab.errorbar(x,y,dy,dx,linestyle = '', color = 'black', marker = '.')
pylab.rc('font',size=18)
pylab.ylabel('X  [digit]')
pylab.minorticks_on()
pylab.plot(xx,scarica2(xx,*pars2), color='red')

#Residui
pylab.subplot(2,3,5)
r = (y-scarica2(x,*pars2))/sigma
# bellurie
pylab.rc('font',size=18)
pylab.ylabel('Norm. res.')
pylab.xlabel('t [ms]')
pylab.minorticks_on()
#pylab.ylim((-.9,.9))
pylab.plot(x,r,linestyle="--",color='blue',marker='o')

print('\n\n\nSCARICA 3 PARAMETRI\n\n\n')
#FIT 3

sigma=dy
w=1/sigma**2
pars3,covm3=curve_fit(scarica3,x,y,init3,sigma)
chi2 = ((w*(y-scarica3(x,*pars3))**2)).sum()
ndof=len(x)-len(init2)

print('\n Senza errore efficace \n')
print('pars3:',pars3)
print('covm3:\n',covm3)
dq, dt, dc=np.sqrt(covm3.diagonal())
print('\ndq=%f , dt=%f' %(dq, dt))
print ('chi23, ndof:',chi2, ndof)

#FIT CON ERRORE EFFICACE

print('\n\nERRORE EFFICACE 3 scarica\n\n')
Deff=np.sqrt(sigma**2+((pars3[0]/pars3[1])*dx*np.exp(-x/pars3[1])**2))

for i in range(2):
    pars3,covm3=curve_fit(scarica3,x,y,pars3,Deff)
    Deff=np.sqrt(sigma**2+((pars3[0]/pars3[1])*dx*np.exp(-x/pars3[1])**2))

sigma=Deff
w=1/sigma**2
pars3,covm3=curve_fit(scarica3,x,y,init3,sigma)
chi2 =((w*(y-scarica3(x,*pars3))**2)).sum()
ndof=len(x)-len(init3)

print('pars:',pars3)
print('covm:\n',covm3)
dq, dr, dc =np.sqrt(covm3.diagonal())
print('\ndc=%f , dr=%f, dc=%f' %(dq, dr, dc))
print ('chi2, ndof:',chi2, ndof)

covnorm=np.zeros((3,3))
#covarianza normalizzata
for i in range(2):
    for j in range(2):
        covnorm[i][j]=covm1[i][j]/(np.sqrt(covm3[i][i]*covm3[j][j]))
print('covnorm:\n', covnorm)


#Plot
pylab.subplot(2,3,3)
pylab.errorbar(x,y,dy,dx,linestyle = '', color = 'black', marker = '.')
pylab.rc('font',size=18)
pylab.ylabel('X  [digit]')
pylab.minorticks_on()
pylab.plot(xx,scarica3(xx,*pars3), color='red')

#Residui
pylab.subplot(2,3,6)
r = (y-scarica3(x,*pars3))/sigma
# bellurie
pylab.rc('font',size=18)
pylab.ylabel('Norm. res.')
pylab.xlabel('t [ms]')
pylab.minorticks_on()
#pylab.ylim((-.9,.9))
pylab.plot(x,r,linestyle="--",color='blue',marker='o')

pylab.show()