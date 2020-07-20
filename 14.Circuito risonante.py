#CIRCUITO RLC IN RISONANZA

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#SCRIVI TUTTO CIO' CHE HAI
C = 0.10e-6	# capacita' nominale del condensatore [microFarad]
R = 354.3	# resistenza inserita tra quelle del banco [ohm]
#Inizializzazione del fit
a = 3.4e-5 #RC
b = 4.0e-5 #(R+r)C
f0 = 717. #[Hz] Frequenza di risonanza
r = 41.1	# resistenza dell'induttanza [ohm]

f, df, Vin, dVin, Vout, dVout = np.loadtxt(r'C:\Users\Satellite\Desktop\Università 2\Laboratorio 2\_ESAME\PR ESAME finale\14.Circuito risonante\data.txt', unpack=True, dtype=float)


A = Vout / Vin	# attenuazione
dA = np.sqrt((dVout/Vin)**2 + (dVin*A/Vin)**2)

#Primo plot con solo dati
plt.figure(0)
plt.title('Dati')
plt.xlabel('$f$ [Hz]', size=14)
plt.ylabel("$A(f)$ [dB]", size=14)
plt.minorticks_on()
plt.xlim(0, np.max(f)+1e2)
plt.ylim(0, 1)
plt.errorbar(f, A, xerr=df, yerr=dA, marker='.', linestyle='none')
plt.show()

df_fwhm_att = 2*np.pi*np.sqrt(3)*(R+r)*C*f0**2	# larghezza a semialtezza [Hz]

#Dichiarazione funzione di fit
def Af(x, a, b, f0):	# attenuazione A(f) con 3 parametri liberi
    return (2*np.pi*a*x)/( np.sqrt( (2*np.pi*b*x)**2 + (1-(x/f0)**2)**2 ) )

X = np.linspace(0.1, max(f)+1e2, 1e5)

# plot curva di risonanza con i parametri iniziali attesi
plt.figure(1)
plt.xlabel('$f$ [Hz]', size=14)
plt.ylabel("$A(f)$", size=14)
plt.minorticks_on()
plt.ylim(0, 1)
plt.plot(X, Af(X, a, b, f0), color='black')
plt.errorbar(f, A, xerr=df, yerr=dA, marker='.', linestyle='none', label='dati')
plt.legend(fontsize=14)


#FIT 
init = [a, b, f0]

sigma = np.sqrt(dA**2 + (a*df*(1 - ((((f/f0)**2-1)*2/(f0**2) + b**2)*(f/np.sqrt((b*f)**2 + (1 - (f/f0)**2)**2))**2)))**2)

popt, pcov = curve_fit(Af, f, A, init, sigma, absolute_sigma = True)

a_3par, b_3par, f0_3par = popt
da_3par, db_3par, df0_3par = np.sqrt(np.diag(pcov))

#Residui normalizzati
res_3par = (A - Af(f, popt[0],popt[1], popt[2]))/sigma

print('Absolute_sigma = True')
print('a = (%.3f +- %.3f)*10e-3 Hz' %(a_3par*1e3, da_3par*1e3))
print('b = (%.3f +- %.3f)*10e-3 Hz' %(b_3par*1e3, db_3par*1e3))
print('f0 = (%.3f +- %.3f) Hz' %(f0_3par, df0_3par))
print('Delta fwhm att= (%.3f)'%(df_fwhm_att))

chi2 = (res_3par**2).sum()
cov = pcov[0,1]/(np.sqrt(pcov[0,0]*pcov[1,1]))

print('\nchi2/ndof = %.1f/%d' %(chi2, len(f)-len(popt)))

covnorm=np.zeros((3,3))
for i in range (3):
    for j in range(3):
        covnorm[i][j]=pcov[i][j]/(np.sqrt(pcov[i][i]*pcov[j][j]))
print('Covarianze normalizzate')

print(covnorm[0][1])
print(covnorm[0][2])
print(covnorm[1][2])

hm_3 = np.max(Af(X, a_3par, b_3par, f0_3par))/2	

# plot curva di risonanza con best-fit
plt.figure(2)
plt.subplot(2,1,1)
plt.title('Best-fit', size=18)
plt.ylabel("$A(f)$", size=14)
plt.minorticks_on()
plt.xlim(0, np.max(f)+1e2)
plt.axhline(y=hm_3, color='black', linestyle='--')
plt.ylim(0, 1)
plt.plot(X, Af(X, a_3par, b_3par, f0_3par), color='black')
plt.errorbar(f, A, xerr=df, yerr=dA, marker='.', linestyle='none')
plt.grid()


#subplot residui
res_2par = (A - g(f, popt[0],popt[1]))/sigma
plt.subplot(2,1,2)
plt.xlabel('$f$ [Hz]', size=14)
plt.ylabel("Norm. res.", size=14)
plt.minorticks_on()
plt.grid()
plt.xlim(0, np.max(f)+1e2)
plt.plot(f, res_3par, color='black', marker='o', linestyle='--')

'''
# plot curva di risonanza in dB
plt.figure(3)
plt.title('Plot $A(f)$ [dB]$', size=18)
plt.xlabel('$f$ [Hz]', size=14)
plt.ylabel("A(f) [dB]", size=14)
plt.minorticks_on()
plt.grid()
plt.xlim(0, np.max(f)+1e2)
plt.ylim(-25, 0)
plt.errorbar(f, 20*np.log10(A), xerr=df, yerr=20*dA/(np.log(10)*A), marker='.', linestyle='none')
plt.plot(X, 20*np.log10(Af(X, a_3par,b_3par, f0_2par)), color='black')
plt.axhline(y= (20*np.log10(hm_3*2)-3), color='black', linestyle='--', label='$A_{max}$-3dB (3 par.)' )

'''


plt.show()

