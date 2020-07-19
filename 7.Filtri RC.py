#esperienza 7
import numpy as np
import pylab
import matplotlib as plt
from scipy.optimize import curve_fit

#acquisizione dati
f, df, A, dA=np.loadtxt(r'C:\Users\Satellite\Desktop\Università 2\Laboratorio 2\_ESAME\lab2 Giulio\esercitazione 7\dati7.txt', unpack=True)

'''Se ho i dati del voltaggio
f, df, Vin, dVin, Vout, dVout=np.loadtxt('filtro dati di merda.txt', unpack=True)

#calcolo guadagno
A=Vin/Vout
dA=np.sqrt((1/Vout*dVin)**2+(Vin/(Vout**2)*dVout)**2)
'''
#Guadagno per il passa basso
def G(f,ft,a):
    return a/np.sqrt(1+(f/ft)**2)

'''
#Funzione di trasferimento per il passa alto
def T(f,ft,a):
    return a/(1-1j(ft/f))
#Guadagno per il passa alto
def G(f,ft,a):
    return a/np.sqrt(1+(ft/f)**2)'''
    
#in decibel
GD= 20*np.log10(A)
dGD=(20/np.log(10))*(dA)/A

init=[382,1]  #frequenza di taglio utlizzata

#Prima figura con dati e funzione calcolata con gli init
pylab.figure(0)
pylab.title('Prova init')
pylab.ylabel('A [UA]')
pylab.xlabel('f [Hz]')
pylab.grid(color = 'gray')
pylab.errorbar(f,A,dA,df,'.', color='orange')
pylab.plot(f, G(f,382,1), color='orange', label='Sbarra 1' )

#Fit
pars, covm=curve_fit(G,f,A,init,dA,absolute_sigma=False)
ft,a= pars
dft,da = pylab.sqrt(covm.diagonal())

pylab.figure(1)
pylab.subplot(211)
pylab.title('Filtro passa')#cambia nome
pylab.ylabel('A [UA]')
pylab.xlabel('f [Hz]')
pylab.grid(color = 'gray')
pylab.errorbar(f,A,dA,df,'.', color='orange')
pylab.loglog(f, G(f,ft,a), color='orange', label='Sbarra 1' )

#residui
pylab.subplot(212)
pylab.title('Residui normalizzati')
pylab.errorbar(f, (A-G(f,ft,a))/dA, dA,df, '.',ls='--', color='orange' )
pylab.xscale('log')
#chi2 e gradi di libertà
chisq = sum(((A-G(f,ft,a))/dA)**2)
ndof=len(f)-2

#Stampa
print('chiquadro/ndof: %f' %(chisq/ndof))
print('frequenza di taglio: %f' %(ft))
print('Coefficente: %f' %(a))
covnorm=np.zeros((2,2))
#covarianza normalizzata
for i in range(2):
    for j in range(2):
        covnorm[i][j]=covm[i][j]/(np.sqrt(covm[i][i]*covm[j][j]))
print('covnorm \n', covnorm)
#Diagramma di Bode
pylab.figure(2)
pylab.ylabel('A [dB]')
pylab.xlabel('f [Hz]')
pylab.errorbar(f,GD,dGD,df ,'.', color='blue')
pylab.xscale('log')

b=1j
#Diagramma di NYQUIST
pylab.figure(3)
fp=np.logspace(-2,5,1000)
#Funzione di trasferimento per il passa basso
T=1/(1+1j*(fp/ft))
'''#Funzione di trasferimento per il passa alto
T=1/(1-1j*(ft/f))'''
pylab.ylabel('Im(T)')
pylab.xlabel('Re(T)')
tr=T.real
ti=T.imag
pylab.plot(tr,ti)
pylab.show()

















