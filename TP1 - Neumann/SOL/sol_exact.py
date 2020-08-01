#!/usr/bin/env python
import numpy as np

# Condition initiale à t=0
def sol_t0(x):
    if x >= 3./8 and x <= 5./8:
        u = 1
    else :
        u = 0
    return u
    
# Solution exacte tronquée
def sol_exacte(x, t, N, L):
    u = 1./4
    for i in range (1, N + 1):
        ipi = i * np.pi
        c = 2. / ipi * (np.sin(5 * ipi / 8) - np.sin(3 * ipi / 8))

        a  = (i**2 * np.pi**2 / L**2)
        m  = i * np.pi / L
        u += c * np.exp(-a * t) * np.cos(m * x)
    return u
