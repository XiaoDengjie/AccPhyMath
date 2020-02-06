import matplotlib.pyplot as plt
import numpy as np
from pylab import mpl
import math

"""
  Gaussian Elimination
"""

t=? #Ratio of A over B
m=? #Truncation
A=0.003*t
B=0.003
l=float(np.sqrt(A**2+B**2))
u=np.arccosh(A/l)
sigma0=4.228e7
tao=8.0055e-12
Z0=377
c=2.997925e8
epsilon=8.8e-12
u1=v1=0
#Def function
sigma=sigma0/(1-I*omega*tao)
k=omega/c
W=i*k*l**2
lambda0=np.sprt(I*k*Z0*sigma0)
Ru0=I*l*lambda0*np.cosh(u)
d[0]=-W/4*np.sinh(2u)-(I*k/(lambda0**2)+I/k)*Ru0
z[0]=W/8*np.sinh(2u)
t[0]=-1/(2Pi*epsilon*c)
s[0]=W/4*np.cosh(2un
s[n]:=W/(8*(2n+1))*(np.sinh(2n*u)/np.sinh((2n+2)u))
z[n]:=W/(8*(2n+1))*(np.sinh((2n+2)u)/np.sinh(2n*u)+np.cosh((2n+2)u)/np.cosh(2n*u)
t[n]:=1/(Pi*epsilon*c)np.cosh(2n*u1)*np.cos(2n*v1)(np.tanh(2n*u)-np.coth(2n*u))
d[n]:=-W(np.sinh((2n+2)u)/((16n+8)np.sinh(2n*u))+np.cosh((2n+2)u)/((16n+8)np.cosh(2n*u))
         +np.sinh((2n-2)u)/((16n-8)np.sinh(2n*u))+np.cosh((2n-2)u)/((16n-8)np.cosh(2n*u)))
        (4nI)/k-I/k*Ru0*np.coth(2n*u)-(I*k)/(lambda0**2)*(Ru0*(np.tanh(2n*u)+np.coth(2n*u))-4n)
                
