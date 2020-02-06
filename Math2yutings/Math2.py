#coding=utf-8
"""
filename:       Math2.py
Description:
Author:         Yuting Wang, Dengjie Xiao
IDE:            PyCharm
Change:         2020/2/6  下午1:04    Dengjie Xiao        Create


"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class GaussElim(object):

    def __init__(self, t, m, omega):
        self.t = t
        self.m = m
        self.omega = omega
        self.i = 1j
        self.A = 0.003*self.t
        self.B = 0.003
        self.l = float(np.sqrt(self.A ** 2 - self.B ** 2))
        self.u = np.arccosh(self.A/self.l)
        self.sigma0 = 4.228e7
        self.tao = 8.0055e-12
        self.Z0 = 377
        self.c = 2.997925e8
        self.epsilon = 8.8e-12
        self.u1 = 0
        self.v1 = self.u1
        # Def function
        self.sigma = self.sigma0 / (1 - self.i * self.omega * self.tao)
        self.k = omega / self.c
        self.W = self.i * self.k * self.l ** 2
        self.lambda0 = np.sqrt(self.i * self.k * self.Z0 * self.sigma0)
        self.Ru0 = self.i * self.l * self.lambda0 * np.cosh(self.u)
        self.d0 = -self.W / 4 * np.sinh(2*self.u)-(self.i * self.k / (self.lambda0 ** 2) + self.i / self.k) * self.Ru0
        self.z0 = self.W / 8 * np.sinh(2*self.u)
        self.t0 = -1 / (2*np.pi * self.epsilon * self.c)
        self.s0 = self.W / 4 * np.cosh(2*self.u)

    def ss(self,n):
        if n == 0:
            ss = self.s0
        else:
            ss = self.W / (8 * (2*n+1)) * (np.sinh(2*n* self.u) / np.sinh((2*n+2)*self.u))
        return (ss)

    def zz(self,n):
        if n == 0:
            zz = self.z0
        else:
            zz = self.W/(8*(2*n+1))*(np.sinh((2*n+2)*self.u)/np.sinh(2*n*self.u)+np.cosh((2*n+2)*self.u)/np.cosh(2*n*self.u))
        return (zz)

    def tt(self, n):
        if n == 0:
            tt = self.t0
        else:
            tt = 1/(np.pi * self.epsilon*self.c)*np.cosh(2*n*self.u1)*np.cos(2*n*self.v1)*(np.tanh(2*n*self.u)-1/np.tanh(2*n*self.u))
        return (tt)

    def dd(self,n):
        if n == 0:
            dd = self.d0
        else:
            dd = -self.W*(np.sinh((2*n+2)*self.u)/((16*n+8)*np.sinh(2*n*self.u))+np.cosh((2*n+2)*self.u)/((16*n+8)*np.cosh(2*n*self.u))   \
                   +np.sinh((2*n-2)*self.u)/((16*n-8)*np.sinh(2*n*self.u))+np.cosh((2*n-2)*self.u)/((16*n-8)*np.cosh(2*n*self.u)))*      \
                      (4*n*self.i)/self.k-self.i/self.k*self.Ru0*(1/np.tanh(2*n*self.u))-(self.i*self.k)/(self.lambda0**2)*(self.Ru0*(np.tanh(2*n*self.u)+1/np.tanh(2*n*self.u))-4*n)
        return (dd)

    def CC(self,n):
        if n == 0:
            CC = self.dd(0)
        else:
            CC = self.dd(n) - (self.ss(n-1)*self.zz(n-1)/self.CC(n-1))
        return (CC)

    def DD(self,n):
        DD = self.ss(n-1)/self.CC(n-1)
        return (DD)

    def T(self,n):
        if n == 0 :
            T = self.tt(0)
        else:
            T = self.tt(n) - self.T(n-1)*self.DD(n)
        return (T)

    def X(self,n):
        if n == self.m :
            X = self.T(self.m)/self.CC(self.m)
        else :
            X = (self.T(n) - self.zz(n)*self.X(n+1))/self.CC(n)
        return (X)

    def Coe(self, ome=None):
        if ome != None:
            omega = ome
            self.__init__(t,m,omega)
        Coe  = 0
        for i in range(self.m+1):
            Coe = Coe + (-1)**i*self.X(i)
        return (Coe)

    def plot(self,OmeRange):
        """
        plot and save data
        """
        data = pd.DataFrame({})
        for i in np.arange(10, OmeRange, 10):
            ome = 10+i
            data.loc[i,'ome'] = ome
            data.loc[i, 'Coe'] =  np.abs(self.Coe(ome))
        data.to_csv("result.csv",index=None)
        plt.plot(np.log(data['ome']),data['Coe'])
        # plt.plot(data['ome'], data['Coe'])
        plt.xlabel('omega')
        plt.ylabel('Coe')
        plt.show()
        return (data)



t = 1.01
m = 10
omega = 10
test = GaussElim(1.01,10,10)
OmegaRange = 1000
test.plot(OmegaRange)




