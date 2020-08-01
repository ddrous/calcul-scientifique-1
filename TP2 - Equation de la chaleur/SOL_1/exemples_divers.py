import time
import numpy as np
from scipy.sparse import diags, linalg as sla
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve, lu_solve, lu_factor

#0. Construction d'une matrice tridiagonale

#Méthode 1, utilisant np.diag()

#Construction d'une matrice tri-diagonale de taille N
N = 4 
#Diagonale principale
diag_main      = 2 * np.ones(N)
diag_main[0]   = 1
diag_main[N-1] = 1

#Diagonales inférieure et supérieure
diag_off = -1 * np.ones(N-1)
diag_sup = -1 * np.ones(N-1) #On en a pas vraiment besoin car diag_off = diag_sup

#Assemblage de A
A  =  np.diag(diag_main, k=0)
A +=  np.diag(diag_off,  k=-1)
A +=  np.diag(diag_off,  k=1)
print("Méthode 1 : utilisation de diag")
print(A)

#Méthode 2, élément par élément
#Création de la matrice nulle de taille NxN
A = np.zeros((N,N))


#Modification des première entrées car diag_main et diag_off démarrent en position  0,0 et 0,1 

A[0][0] = diag_main[0]
A[0][1] = diag_off[0]

for i in range(1,N-1): #Boucle de 1 à N-2 inclus
    A[i][i] = diag_main[i]
    A[i][i-1] = diag_off[i]
    A[i][i+1] = diag_sup[i]

#Modification des dernières entrées car diag_main et diag_off vont jusqu'à (N-1, N-1) et (N-1,N-2)
A[N-1][N-1] = diag_main[N-1]
A[N-1][N-2] = diag_sup[N-2]
print("Méthode 2 : élement par élement")
print(A)


#Méthode 3. Utilisant le format matrice creuse 'sparce' de scipy. Ici on ne stocke pas les 0.
diag_main      = 2 * np.ones(N).ravel()
diag_main[0]   = 1
diag_main[N-1] = 1

#Diagonales inférieure et supérieure
diag_off = -1 * np.ones(N-1).ravel()
diag_sup = -1 * np.ones(N-1).ravel() #On en a pas vraiment besoin car diag_off = diag_sup

#Assemblage de A
diagonals = [diag_off, diag_main, diag_sup]
A = diags(diagonals, [-1,0,1], shape=(N, N),format='csc')
print("Méthode 3 : utilisant scipy sparse")
print(A.toarray())

#1. Exemple d'un produit matrice vecteur u = A.v
print("Produit matrice vecteur u = A.v")
v = np.random.randint(100, size=N) #Création d'un vecteur aléatoire pour les tests
u = A.dot(v)
print(u)

# 2. Exemple de résolution d'un système linéaire
print("Résolution d'un système linéaire u = A.v")
#Construction d'une matrice non creuse
diag_main = np.random.randint(100, size=N)
diag_off  = np.random.randint(100, size=N-1)
diag_sup  = np.random.randint(100, size=N-1)

#Assemblage de A
A  =  np.diag(diag_main, k=0)
A +=  np.diag(diag_off,  k=-1)
A +=  np.diag(diag_off,  k=1)

print("Matrice pleine")
#Méthode frontale
x = solve(A, v)
print(x)
#Méthode LU
LU = lu_factor(A)
x =  lu_solve(LU, v)

print(x)


#Méthode LU si vous avez utilisé des matrices creuses
#Construction d'un matrice M creuse tridiagonale
diag_main      = np.random.randint(100, size=N).ravel()
diag_off       = np.random.randint(100, size=N-1).ravel()
diag_sup       = np.random.randint(100, size=N-1).ravel()
diagonals = [diag_off, diag_main, diag_sup]
M = diags(diagonals, [-1,0,1], shape=(N, N),format='csc')
print("Matrice creuse M")
print(M.toarray())

#Méthode frontale pour matrice creuse
x = spsolve(M,u)
print(x)


#Méthode LU pour matrice creuse
LU = sla.splu(M)
x = LU.solve(u)
print(x)







