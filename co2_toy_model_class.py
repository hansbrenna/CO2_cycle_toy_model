# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 14:06:53 2016

@author: hanbre
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

mco2 = 44.0

def flux(cc1,cc2,k,c):
    f = abs(cc1-cc2)*k*c
    return f

class reservoir():
    def __init__(self,total_m,m_c,mr):    
        self.total_mass=total_m
        self.mass_c0=m_c
        self.mass_c=0
        self.mres=mr
        self.c0=self.concentration(self.mass_c0,self.total_mass,mco2,self.mres)
        self.c = 0
        self.cl = [self.c0]
        self.ml = [self.mass_c0]
        
    def concentration(self,mass_c,total_mass,mco2,mres):
        conc=mass_c/total_mass*mco2/mres
        return conc
        
    def store_previous(self):
        self.cl.append(self.c)
        self.ml.append(self.mass_c)
        self.c0=self.c
        self.mass_c0 = self.mass_c
        
    def update(self,F,dt):
        self.mass_c = self.mass_c0+sum(F)*dt
        self.c=self.concentration(self.mass_c,self.total_mass,mco2,self.mres)
        self.store_previous()

Mocn = 1.4e18*1000
    
atm = reservoir(4.5e18,589e12*3.665,29)
socn = reservoir(0.025*Mocn,900e12*3.665,18)
docn = reservoir(0.975*Mocn,37100e12*3.665,18)
btm = reservoir(1.0e30,1e25,1)

kao=1
koa=1
ksd=1
kds=1
kdb=1
kbd=0.0

Fao = flux(atm.c0,socn.c0,kao,atm.mass_c0)
Foa = flux(atm.c0,socn.c0,koa,socn.mass_c0)
Fsd = flux(socn.c0,docn.c0,ksd,socn.mass_c0)
Fds = flux(socn.c0,docn.c0,kds,docn.mass_c0)
Fbd = flux(btm.c0,docn.c0,kbd,btm.mass_c0)
Fdb = flux(btm.c0,docn.c0,kdb,docn.mass_c0)

Fao1 = [Fao]
Foa1 = [Foa]
Fsd1 = [Fsd]
Fds1 = [Fds]
Fbd1 = [Fbd]
Fdb1 = [Fdb]

print(Fao,Foa,Fsd,Fds,Fdb,Fbd)

t0=0
t=t0
dt = 1
tstop = 2000000

i=0
while t < tstop:
    #if i == np.floor(tstop/dt/2.0):
    #    atm.mass_c0+=5000e12
    #    print('here')
    t+=dt
    atm.update([-Fao,Foa],dt)
    socn.update([-Foa,Fao,-Fsd,Fds],dt)
    docn.update([Fsd,-Fds,Fbd,-Fdb],dt)
    btm.update([Fdb,-Fbd],dt)
    Fao = flux(atm.c,socn.c,kao,atm.mass_c)
    Foa = flux(atm.c,socn.c,koa,socn.mass_c)
    Fsd = flux(socn.c,docn.c,ksd,socn.mass_c)
    Fds = flux(socn.c,docn.c,kds,docn.mass_c)
    Fbd = flux(btm.c,docn.c,kbd,btm.mass_c)
    Fdb = flux(btm.c,docn.c,kdb,docn.mass_c)
    Fao1.append(Fao)
    Foa1.append(Foa)
    Fsd1.append(Fsd)
    Fds1.append(Fds)
    Fdb1.append(Fdb)
    Fbd1.append(Fbd)
    i+=1
   
print(Fao1[-1],Foa1[-1],Fsd1[-1],Fds1[-1],Fdb1[-1],Fbd1[-1])
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
ax1.semilogy(np.array(atm.ml))
ax2.semilogy(np.array(socn.ml))
ax3.semilogy(np.array(docn.ml))
ax4.semilogy(np.array(atm.cl))
plt.show()

print(Fao1[-1]/atm.ml[-1])
