import numpy as np
import matplotlib as mlpt
import matplotlib.pyplot as plt
import pylab 
from scipy.optimize import curve_fit

t,V=pylab.loadtxt(r'C:\Users\aministratore\Documents\Università\Secondo anno\Laboratorio 2\File elearning\Esperienza 12\Dati normalista\data1.txt',unpack=True)


DV=np.array([1]*len(V))
Dt=4e-6
pylab.figure(1)
pylab.errorbar(t,V,DV,Dt,linestyle = '', color = 'black', marker = '.') #per stimare init


C=0.1e-6
#A e p vanno stimati dal grafico (A è l'ampiezza, p la trovi a t=0)

def funz(t, w, tt, A, p, c):
    
    return A*np.exp(-t/tt)*np.cos(w*t+p)+c

init=(0.0021,60000, 400, np.pi/2, 450)
xx=np.linspace(min(t),max(t),2000)

pylab.plot(xx,funz(xx,*init), color='blue')  #plot con pars

# set the error
#sigma=DV #prima chiamata
sigma=pars[3]*np.exp(-t/pars[1])*((1/pars[1])*np.cos(pars[0]*t+pars[4])+pars[0]*np.sin(pars[0]*t+pars[4])) #seconda
w=1/sigma**2

# call the minimization routine
pars,covm=curve_fit(funz,t,V,init,sigma)

# calculate the chisquare for the best-fit function
chi2 = ((w*(V-funz(t,*pars))**2)).sum()

# determine the ndof
ndof=len(V)-len(init)

# print results on the console
print('pars:',pars)
print('covm:\n',covm)
print ('chi2, ndof:',chi2, ndof)
pylab.figure(2)
pylab.subplot(2,1,1)
# plot the best fit curve
pylab.plot(xx,funz(xx,*pars), color='red')
pylab.errorbar(t,V,DV,Dt,linestyle = '', color = 'black', marker = '.') #per vedere come sono
'''
inutile perché dV>>dt e non ha senso propagare dt su dv

## errore efficace

Deff=np.sqrt(sigma**2+(init[2]*np.exp(-t/init[1])*((1/init[1])*np.cos(init[0]*t+init[3])+init[0]*np.sin(init[0]*t+init[3]))*Dt)**2)

#errore efficace
for i in range (2):
    pars,covm=curve_fit(funz,t,V,init,sigma)
    Deff=np.sqrt(sigma**2+(init[2]*np.exp(-t/init[1])*((1/init[1])*np.cos(init[0]*t+init[3])+init[0]*np.sin(init[0]*t+init[3]))*Dt)**2)

sigma=Deff
w=1/sigma**2
parse,covme=curve_fit(funz,t,V,init,sigma)
chi2 = ((w*(V-funz(t,*parse))**2)).sum()
ndof=len(V)-len(init)

print('\n\n ERRORE EFFICACE')
print('pars:',parse)
dw, dtt, dA, dp, dc=np.sqrt(covme.diagonal())
print('\ndw=%f , dtt=%f, dA=%f, dp=%f, dc=%f' %(dw, dtt, dA, dp, dc))
print('covm:\n',covme)
print ('chi2, ndof:',chi2, ndof)
covnorm=np.zeros((5,5))
#covarianza normalizzata
for i in range(5):
    for j in range(5):
        covnorm[i][j]=covme[i][j]/(np.sqrt(covme[i][i]*covme[j][j]))


print('covnorm:\n', covnorm)
print('\n\n')
plt.figure(3)
pylab.plot(xx,funz(xx,*pars), color='red')
pylab.errorbar(t,V,Deff,Dt,linestyle = '', color = 'black', marker = '.') #per vedere come son
'''
##residui

# switch to the residual plot
pylab.figure(2)
pylab.subplot(2,1,2)
pylab.grid()
# build the array of the normalized residuals
r = (V-funz(t,*pars))/sigma

# bellurie
pylab.rc('font',size=18)
pylab.ylabel('Norm. res.')
pylab.minorticks_on()
# set the vertical range for the norm res
#pylab.ylim((-.9,.9))

# plot residuals as a scatter plot with connecting dashed lines
pylab.plot(t,r,linestyle="--",color='blue',marker='o')

# show the plot

pylab.show()

##outliers


xx2=np.array([])
yy2=np.array([])
Dx2=np.array([])
Dy2=np.array([])

outx=np.array([])
outy=np.array([])
Dxout=np.array([])
Dyout=np.array([])

j=0
k=0

for i in range(len(V)):
    if (np.abs(V[i]-funz(t[i], *pars))<5*DV[i]):
        xx2=np.insert(xx2, j, t[i])
        yy2=np.insert(yy2, j, V[i])
        Dx2=np.insert(Dx2, j, Dt)
        Dy2=np.insert(Dy2, j, DV[i])
        j+=1
    else:
        outx=np.insert(outx, k, t[i])
        outy=np.insert(outy, k, V[i])
        Dxout=np.insert(Dxout, k, Dt)
        Dyout=np.insert(Dyout, k, DV[i])
        k+=1
    
sigma=Dy2
pars,covm=curve_fit(funz,xx2,yy2,init,sigma)

w=1/sigma**2

# calculate the chisquare for the best-fit function
chi2 = ((w*(yy2-funz(xx2,*pars))**2)).sum()

# determine the ndof
ndof=len(yy2)-len(init)

# print results on the console
print('pars:',pars)
print('covm:\n',covm)
print ('chi2, ndof:',chi2, ndof)
pylab.figure(4)
pylab.subplot(2,1,1)
# plot the best fit curve
pylab.plot(xx,funz(xx,*pars), color='blue')
pylab.errorbar(xx2,yy2,Dy2,Dx2,linestyle = '', color = 'black', marker = '.')
pylab.errorbar(outx,outy,Dyout,Dxout,linestyle = '', color = 'red', marker = '.')

pylab.subplot(2,1,2)
pylab.grid()
# build the array of the normalized residuals
r = (yy2-funz(xx2,*pars))/sigma
rout=(outy-funz(outx,*pars))/Dyout
# bellurie
pylab.rc('font',size=18)
pylab.ylabel('Norm. res.')
pylab.minorticks_on()


# plot residuals as a scatter plot with connecting dashed lines
pylab.plot(xx2,r,linestyle="--",color='blue',marker='o')
pylab.plot(outx,rout,linestyle="",color='red',marker='o')

# show the plot
pylab.show()























