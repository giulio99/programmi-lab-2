import pylab
import numpy 
from scipy.optimize import curve_fit

# data load
x,y=pylab.loadtxt(r'C:\Users\GIULIO\Desktop\università\aa2019-2020\Laboratorio 2\esercitazione 5\ondegiulio0.txt',unpack=True)

#converto i dati da microsecondi a secondi in questa maniera mi evito eventuali cambi di unità di misura più avanti nel fit, soprattuto nel passaggio da tempi a frequenze
x=x/100000

dx=numpy.array([4/100000]*len(x))
dy=numpy.array([1]*len(y))

pylab.figure(1)
pylab.subplot(2,1,1)

# scatter plot with error bars
pylab.errorbar(x,y,dy,dx,linestyle = '--', color = 'black', marker = '.')

# bellurie 
pylab.rc('font',size=18)
pylab.ylabel('$V$  [digit]')
pylab.minorticks_on()

pylab.show()

#definisci la funzione di x con 4 parametri, rispettivamente ampiezza, frequenza(ANGOLARE), sfasamento e offset
def ff(x, a, w,d, c):
    return a*numpy.sin(w*x+d)+c

#dalla figura vedi l'ampiezza come metà della distanza tra max e min, frequenza come 2pi*(numero di picchi in un secondo)
#(approssimata, fai il calcolo), fase va bene zero e offset come valore medio dei punti in figura uno

init=(400, 250, -0.8,500)

# prepare a dummy xx array (with 2000 linearly spaced points)
xx=numpy.linspace(min(x),max(x),2000)

# plot the fitting curve computed with initial values
# AT THE SECOND ATTEMPT THE FOLLOWING LINE MUST BE COMMENTED 



pylab.figure(2)
pylab.errorbar(x,y,dy,dx,linestyle = '--', color = 'black', marker = '.')
pylab.plot(xx,ff(xx,*init), color='blue') 

pylab.show()

# set the error
sigma=dy


# call the minimization routine
pars,covm=curve_fit(ff,x,y,init,sigma)

# calculate the chisquare for the best-fit function
w=1/sigma**2
chi2 = ((w*(y-ff(x,*pars))**2)).sum()

# determine the ndof
ndof=len(x)-len(init)

#determina la frequenza dell'onda dalla frequenza angolare del fit e l'errore
f=pars[1]/2*numpy.pi
df=numpy.sqrt(covm[1][1])
# print results on the console
print('f= %f +- %f' %(f, df))
print('pars:',pars)
print('covm:',covm)
print ('chi2, ndof:',chi2, ndof)

covnorm=np.zeros((4,4))
#covarianza normalizzata
for i in range(4):
    for j in range(4):
        covnorm[i][j]=covm[i][j]/(np.sqrt(covm[i][i]*covm[j][j]))
print('covnorm \n', covnorm)
        

pylab.figure(1)
# plot the best fit curve
pylab.subplot(2,1,1)
pylab.errorbar(x,y,dy,dx,linestyle = '--', color = 'black', marker = '.')
pylab.plot(xx,ff(xx,*pars), color='red') 

# switch to the residual plot
pylab.subplot(2,1,2)

#calcolo l'array dei residui normalizzati
r = (y-ff(x,*pars))/sigma

# bellurie 
pylab.rc('font',size=18)
pylab.ylabel('Norm. res.')
pylab.xlabel('$t$  [ms]')
pylab.minorticks_on()

#plotta i residui
pylab.plot(x,r,linestyle="--",color='blue',marker='o')


# show the plot
pylab.show()


