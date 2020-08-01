#!/usr/bin/env python
import time
import numpy as np
# from sol_exact import *
from scipy.sparse import diags
import matplotlib.pyplot as plt

#Solution exacte
def sol_exacte(x, t):
    if t < 1.:
        if x <= t:
            u = 1.
        else:
            if x >= 1.:
                u = 0.
            else:
                u = (1. - x) / (1. - t)
    else:
        if x < (t + 1.) / 2.:
            u = 1.
        else:
            u = 0.

    return u
    
    
#Solution initiale
def sol_t0(x):
    return sol_exacte(x, 0)


#Génération du maillage
def gen_mesh(xmin, xmax, Nx):
    dx = np.abs(xmax - xmin) / Nx 
    x = np.zeros(Nx + 2)
    
    for i in range(Nx + 2):
        x[i] = xmin + (i - 0.5) * dx
        
    return x, dx

#Flux physique
def flux(u):
    return u * u / 2.


#Dérivée du flux physique
def dflux(u):
    return u 


#Flux de Rusanov
def flux_rusanov(uL, uR):
    vmax = max(abs(dflux(uL)), abs(dflux(uR)))
    return (flux(uL) + flux(uR)) / 2.0  - vmax * (uR - uL) / 2.0

 
 


#Résolution VF-Rusanov 
def sol_rusanov(x, dx, Nx, tmax, cfl):
    #Initialisation
    tnow = 0
    u = np.asarray([sol_t0(x_i) for x_i in x])
    up1 = np.zeros(Nx + 2)
    
    while tnow < tmax:
        #Recherche de la vitesse max sur tout le domaine
        vmax = 0
        for i in range(Nx + 2):
            vloc = abs(dflux(u[i]))
            if vloc > vmax :
                vmax = vloc

        dt = cfl * dx / vmax

        #Résolution en espace
        for i in range(1, Nx + 1):
            up1[i] = u[i] - dt / dx * (flux_rusanov(u[i], u[i+1]) - flux_rusanov(u[i-1], u[i])) 
        
        #Recopie 
        for i in range(1, Nx + 1):
            u[i] = up1[i]
            
        tnow += dt
            
    return u


xmin = - 1.
xmax = 2.
Nx = 256
tmax = 1
cfl = 0.99

#Discrétisation en espace
x, dx = gen_mesh(xmin, xmax, Nx)


#Solution exacte
u_rus  = sol_rusanov(x, dx, Nx, tmax, cfl)
u_exacte   = [sol_exacte(xi, tmax) for xi in x]

#Tracé des solution
plt.plot(x, u_exacte, label='Sol. exacte')
plt.plot(x, u_rus, label='Sol. Rusanov')
plt.legend(loc="upper right", fontsize=18)
plt.show()


