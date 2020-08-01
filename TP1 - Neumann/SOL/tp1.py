#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from sol_exact import *

#Paramètres du calcul
Np = 100
L  = 1.
t  = 0
N  = 100
xp = np.linspace(0., L, Np)

#Génération des solutions (discrète). Liste par compréhension
u_t0    = [sol_t0(x) for x in xp]
u_exact = [sol_exacte(x, t, N, L) for x in xp]

plt.plot(xp, u_t0, label='Condition initiale')
plt.plot(xp, u_exact, label='Série tronquée')
plt.legend(loc="upper right", fontsize=18)
plt.show()
