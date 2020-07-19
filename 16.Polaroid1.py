import pylab
import numpy
from scipy.optimize import curve_fit

# data load
x,y=pylab.loadtxt(r'C:\Users\Satellite\Desktop\Università 2\Laboratorio 2\_ESAME\PR ESAME finale\16.Polarizzazione luce\data1.txt',unpack=True)
Dx=0.01
Dy=np.array([0.1]*len(y))

# define the linear function
#nota che la tua x e' RJ+rA (nel terzo fit)
def ff(x,A,B,C):
    x = x/180.*np.pi
    C = C/180.*np.pi
    return A*(np.cos(x-C))**2+B #Togli la fase se
    
# make the array with initial values
init=(60,-5, 60) #ampiezza, dove tocca l'asse y (fase verticale), fase (fase orizzontale)

#Per controllare i parametri iniziali
pylab.figure(0)
pylab.errorbar(x,y,Dy,Dx,linestyle = '', color = 'black', marker = '+')
xx=numpy.linspace(min(x)-5.,max(x)+5.,1000)
pylab.plot(xx,ff(xx,*init), color='red')


#Fit
pylab.figure(1)
pylab.subplot(2,1,1)
pylab.errorbar(x,y,Dy,Dx,linestyle = '', color = 'black', marker = '+')
pylab.plot(xx,ff(xx,pars[0],pars[1],pars[2]), color='red')

# bellurie 
pylab.rc('font',size=16)
pylab.xlabel(' $\Theta_{0} $ [degree]')
pylab.ylabel('$I_{ph}$ [$\mu$A]') 
pylab.title('Data plot w numerical fit')
pylab.minorticks_on()


# set the error
#nota che prima sigma=Dy per trovare i parametri iniziali e poi fare sigma aggiungendo in quadratura l'errore su x
sigma=Dy
w=1/sigma**2

# call the routine (NOTE THE ABSOLUTE SIGMA PARAMETER!)
pars,covm=curve_fit(ff,x,y,init,sigma)

# calculate the chisquare for the best-fit funtion
# note the indexing of the pars array elements
chi2 = ((w*(y-ff(x,pars[0],pars[1],pars[2]))**2)).sum()

# determine the ndof
ndof=len(x)-len(init)

# print results on the console
print(pars)
print(covm)
print (chi2, ndof)

# print the same in a slightly more readible version (not very useful, yet)
print('a = ', pars[0], '+/-', numpy.sqrt(covm[0,0]))
print('b = ', pars[1], '+/-', numpy.sqrt(covm[1,1]))
print('c= ',pars[2],'+/-',numpy.sqrt(covm[2,2]))
print('norm cov = ', covm[0,1]/(numpy.sqrt(covm[0,0]*covm[1,1])))
print('norm cov2 =',covm[1,2]/(numpy.sqrt(covm[1,1]*covm[2,2])))
print('norm cov3 =',covm[0,2]/(numpy.sqrt(covm[0,0]*covm[2,2])))


# prepare a dummy xx array (with 100 linearly spaced points)
xx=numpy.linspace(min(x)-5.,max(x)+5.,1000)

#residui normalizzati
pylab.subplot(2,1,2)
r=(y-ff(x,pars[0],pars[1],pars[2]))/sigma
pylab.ylabel('Norm. res.')
pylab.minorticks_on()
pylab.plot(x,r,linestyle="--",color='blue',marker='o')
pylab.show()

# show the plot
pylab.show()