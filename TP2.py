# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 08:39:43 2021

@author: elebr

TP n°2 Méthode de Cholesky pour la résolution de systèmes linéaires
Jimmy Hoarau - Emeric Lebreton
"""
import numpy as np
import time as t
import matplotlib.pyplot as plt
from math import *

# 1- Décomposition de Cholesky

# Question 1

#A doit être inversible, définie et positive
#On cherche L une matrice telle que A = L(LT)
"""
def Cholesky(A):
    #rend la matrice L de la décomposition de Cholesky de A, une matrice symétrique déﬁnie positive (avec donc A = LL)
    n, m = A.shape
    L = np.zeros((n,n))
    for i in range (0,n):
        for k in range (0,n):
            if i==k : 
                #Calcul des coefficients sur  la Diagonale
                somme1 = 0
                for j in range(0,i+1):
                    somme1 = somme1 + (L[i,j])**2
                L[i,k] = (A[i,j] - somme1)**0.5
            if i > k :
                #Calcul coefficients hors diagonaux
                somme2 = 0
                for j in range (0,i+1) :
                    somme2 = somme2 + (L[i,j]*L[k,j])
                L[i,k] = ((A[i,k] - somme2) / L[k,k])
    #print (L)
    return L

#A = np.array([[4,-2,-4], [-2,10,5], [-4,5,6]])   
#Cholesky (A)
"""
def Cholesky(A):
    n,n = np.shape(A)
    L = np.zeros((n,n))

    for i in range(0, n):
        S = 0

        for j in range(0, i):
            S += L[i, j] ** 2
        L[i,i] = sqrt(A[i, i] - S)

        for k in range(i + 1, n):
            S = 0

            for j in range(0, i):
                S += L[k, j] * L[i, j]

            L[k, i] = (A[k, i] - S) / L[i, i]

    return L
"""
def  npcholesky(A) :
    L = np.linalg.cholesky(A)
    #print("On obtient la matrice L = ",L)
    return(L)
"""
#A = np.array([[4,-2,-4], [-2,10,5], [-4,5,6]])
#npcholesky(A)

# 2- Résolution de systèmes à l’aide de la décomposition de Cholesky

# Question 1

#On cherche AX= B
#On pose A = L(LT)
#Puis on pose Y = (LT)X
#Donc LY = B
#On peut calculer le déterminant de A facilement

def ResolutionSystTriInf(Taug):
    n,m=Taug.shape
    X=[]
    X.append(1)
    for i in range(n-1): #création de la matrice colonne X solution (contenant que des 0 pour le moment)
        X.append(0)
       
    for i in range(0,n):#balayage ligne 0 à N
        for j in range(0,n):  #balayage terme de la ligne en question
            if i > j : # si le coeff X[i] solution est connu,
                A = Taug[i,j] * X[j]        #il est multiplié au terme correspondant dans la matrice
                Taug[i,m-1] = Taug[i,m-1] - A # puis soustrait à B[i]
            if i == j :     #nous sommes sur la diagonale,il suffit donc de diviser le B[i] par ce terme
                Taug[i][m-1] = Taug[i][m-1] / Taug[i,i]
        X[i] = Taug[i,m-1]
    return(X)
    
    
    
def ResolutionSystTriSup(Taug):
    n,m=Taug.shape
    X=[]
    for i in range(n-1): #création de la matrice colonne X solution (contenant que des 0 pour le moment)
        X.append(0)
    X.append(1)
    for i in range(n-1,-1,-1):#balayage ligne N à 0
        for j in range(n-1,-1,-1):  #balayage terme de la ligne en question
            if i < j : # si le coeff X[i] solution est connu,
                A = Taug[i,j] * X[j]        #il est multiplié au terme correspondant dans la matrice
                Taug[i,m-1] = Taug[i,m-1] - A # puis soustrait à B[i]
            if i == j :     #nous sommes sur la diagonale,il suffit donc de diviser le B[i] par ce terme
                Taug[i][m-1] = Taug[i][m-1] / Taug[i,i]
        X[i] = Taug[i,m-1]
    return(X)
"""

def ResolCholesky(L,B) :
    
    #On résout LY = B
    L = Cholesky(A)
    n, m = A.shape
    Laug = np.insert(L,n,B.T,axis = 1)
    Y = ResolutionSystTriInf(Laug)
    #print("On obtient la matrice Y = ",Y)
    
    #On résout (LT)X = Y
    Lt = np.transpose(L)
    Ltaug=np.insert(Lt,n,Y, axis=1)
    X = ResolutionSystTriSup(Ltaug)
    #print("On obtient la matrice X = ",X)
    return(X)

def ResolCholesky(L,B) :
    L = Cholesky(A)
    n, m = A.shape
    M1 = np.column_stack([L, B.T])
    Y = ResolutionSystTriInf(M1)
    N = np.column_stack([L.T, Y])
    X = ResolutionSystTriSup(N)
    #print(X)
    return (X)
#A = np.array([[4,-2,-4], [-2,10,5], [-4,5,6]])   
#L = np.array([[2,0,0], [-1,3,0], [-2,1,1]])
#B = np.array([6,-9,-7])
#ResolCholesky(L,B)

def ResolnpCholesky(L,B) :
    
    #On résout LY = B
    L = np.linalg.cholesky(A)
    n, m = A.shape
    Laug = np.insert(L,n,B,axis = 1)
    Y = ResolutionSystTriInf(Laug)
    #print("On obtient la mtrice Y = ",Y)
    
    #On résout (LT)X = Y
    Lt = np.transpose(L)
    Ltaug=np.insert(Lt,n,Y, axis=1)
    X = ResolutionSystTriSup(Ltaug)
    #print("On obtient la mtrice X = ",X)
    return(X)
"""
#ResolnpCholesky(L,B)


def ResolCholesky(A,B):
    n, n = np.shape(A)
    B = B.reshape(n,1)
    L = Cholesky(A)
    Lt = np.transpose(L)
    x = np.zeros(n)
    y = []

    #Pour: Ly = b
    for i in range(0,n):
        y.append(B[i])

        for k in range(0,i):
            y[i] = y[i] - L[i, k] * y[k]
        y[i] = y[i] / L[i,i]

    #Pour: LTx = y
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(Lt[i, i + 1:], x[i + 1:])) / Lt[i, i]

    return x

def ResolnpCholesky(A,B):
    n, n = np.shape(A)
    B = B.reshape(n,1)
    L = np.linalg.cholesky(A)
    Lt = np.transpose(L)
    x = np.zeros(n)
    y = []

    #Pour: Ly = b
    for i in range(0,n):
        y.append(B[i])

        for k in range(0,i):
            y[i] = y[i] - L[i, k] * y[k]
        y[i] = y[i] / L[i,i]

    #Pour: LTx = y
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(Lt[i, i + 1:], x[i + 1:])) / Lt[i, i]

    return x

def  Resolution_solve(A,B) :
    X = np.linalg.solve(A,B)
    #print("On obtient la mtrice X = ",X)
    return(X)

#A = np.array([[4,-2,-4], [-2,10,5], [-4,5,6]])
#B = np.array([6,-9,-7])
#Resolution_solve(A,B)


"""

def choleskyalternative ( A ):
    n,m = A.shape
    L = np.eye(n) # On initialise L comme identité pour avoir les termes diagonaux .
    D = np.eye(n)
    for k in range (n):
    # D’abord le calcul de dkk en commençant par la somme.
        S =0.
        for j in range (k):
            S = S + L [k,j]**2* D[j,j]
        D[k,k]= A[k,k] - S
        # Puis le calcul de la colonne k de L .
        # On calcule lik pour chaque i > k .
        for i in range (k+1,n):
            S = 0
            for j in range (k):
                S = S + L[i,j]*L[k,j]*D[j,j]
            L[i,k]=(A[i,k]-S)/D[k,k]
    print("On obtient L :",L,"\n et D :",D)
    return (L,D)
    


#A = np.array([[4,-2,-4], [-2,10,5], [-4,5,6]])
#B = np.array([6,-9,-7])
#choleskyalternative(A)
"""

def  ReductionGauss(Aaug) :
    m,n = Aaug.shape                                 
    if m != n-1:
        print("On ne peux pas effectuer la fonction. La matrice doit être augmentée pour obtenir un résultat")
        return
    A = np.copy(Aaug)                                
    for i in range(0,m-1):
        l = A[i,i]                                   
        if (l == 0):                                
            print("La réduction de A n'existe pas")
            return
        else:
            for j in range(i+1,m):
                g = A[j,i]/l                         
                A[j,:] -= g*A[i,:]                   
    return(A)
    
def Gauss(A,B):
    m1,n1 = A.shape
    m2,n2 = B.shape
    if m1 != n1:
        print("La matrice A n'est pas carrée, le calcul n'est donc pas possible")
        return
    if n2 != 1 or m2 != m1:
        print("La matrice B n'est pas aux bonnes dimensions, le calcul n'est pas possible")
        return
    Aaug = np.concatenate((A,B),axis=1)
    Taug =  ReductionGauss(Aaug)
    X = ResolutionSystTriSup(Taug)
    return(X)

#A = np.array([[4,-2,-4], [-2,10,5], [-4,5,6]])
#B = np.array([6,-9,-7])
#Gauss(A,B)


#Vérification de la solution :
#np.allclose(np.dot(A, X), B)

# 3- Expérimentation des méthodes

# Question 1



Taillematrice = []
for i in range(10,300,10):
    Taillematrice.append(i)

len_Taillematrice = len(Taillematrice)

Erreur_Cholesky = []
Erreur_GaussLU = []
Erreur_solve = []
#Erreur_CholeskyAlternatif = []
Erreur_npcholesky = []

Temps_Cholesky = []
Temps_GaussLU = []
Temps_solve = []
#Temps_CholeskyAlternatif = []
Temps_npcholesky = []

for i in range(len_Taillematrice):
    n = Taillematrice[i]
    M = np.random.rand(n,n)
    B = np.random.rand(n,1)
    A = np.dot(M.T,M)
    A1 = np.copy(A)
    A2 = np.copy(A)
    A3 = np.copy(A)
    A4 = np.copy(A)
    
    start = t.process_time()
    X = ResolCholesky(A1,B)
    stop = t.process_time()
    #E = np.linalg.norm(A@X-np.ravel(B))
    E = np.linalg.norm(np.dot(A1,X)- np.ravel(B))
    Erreur_Cholesky.append(E)
    Temps_Cholesky.append(stop - start)
    
    start = t.process_time()
    X = Resolution_solve(A2,B)
    stop = t.process_time()
    #E = np.linalg.norm(A@X-np.ravel(B))
    E = np.linalg.norm(np.dot(A2,X)-B)
    Erreur_solve.append(E)
    Temps_solve.append(stop - start)
    
    start = t.process_time()
    X = ResolnpCholesky(A3,B)
    stop = t.process_time()
    #E = np.linalg.norm(A@X-np.ravel(B))
    E = np.linalg.norm(np.dot(A3,X)-np.ravel(B))
    Erreur_npcholesky.append(E)
    Temps_npcholesky.append(stop - start)
    
    start = t.process_time()
    X = Gauss(A4,B)
    stop = t.process_time()
    #E = np.linalg.norm(A@X-np.ravel(B))
    E = np.linalg.norm(np.dot(A4,X)-np.ravel(B))
    Erreur_GaussLU.append(E)
    Temps_GaussLU.append(stop - start)
    
    
plt.plot(Taillematrice, Temps_Cholesky, label = "Cholesky")
plt.plot(Taillematrice, Temps_npcholesky, label = "npcholesky")
plt.plot(Taillematrice, Temps_GaussLU, label ="GaussLU")
plt.plot(Taillematrice, Temps_solve, label ="np.linalg.solve")
#plt.plot(Taillematrice, Temps_CholeskyAlternatif, label ="Cholesky Alternatif")
plt.title("Temps de calcul CPU en fonction de la taille de la matrice")
plt.xlabel("Taille de la matrice")
plt.ylabel("Temps de calcul en (s)")
plt.legend()
plt.show()
    
# Question 2

plt.plot(Taillematrice, Erreur_Cholesky, label = "Cholesky")
plt.plot(Taillematrice, Erreur_npcholesky, label = "npcholesky")
plt.plot(Taillematrice, Erreur_GaussLU, label ="GaussLU")
plt.plot(Taillematrice, Erreur_solve, label ="np.linalg.solve")
#plt.plot(Taillematrice, Erreur_CholeskyAlternatif, label ="Cholesky Alternatif")
plt.title("Erreurs de calcul en fonction de la taille de la matrice", fontsize = 10)
plt.xlabel("Taille de la matrice")
plt.ylabel("||AX - B||")
plt.legend()
plt.show()
