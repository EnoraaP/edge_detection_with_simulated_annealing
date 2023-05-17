# -*- coding: utf-8 -*-
"""
Created on Thu May 20 16:11:34 2021

@author: Enora
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import special

Px = [-2, 0, 5, 8]
Py = [-4, 2, 4, -3]

Px2 = [2,-2, 1, 6, 8, 6]
Py2 = [-4, -1, 6, 5, 2, -2]

NB_POINTS = 100

couleurs = ['red', 'darkorchid','orange', 'limegreen', 'royalblue', 'darkorchid','lightseagreen', 'deepskyblue']
#couleurs = ['navy', 'slateblue', 'darkslategrey', 'purple','dodgerblue',  'steelblue', 'darkturquoise','indigo' ]
# Interpolation de Lagrange

def Lagrange(Px, Py):
    n = len(Px)
    X = np.linspace(np.amin(Px), np.amax(Px), 50)
    P = 0
    for i in range(n):
        L = 1
        for j in range(0, i):
            L = L * (X - Px[j]) / (Px[i] - Px[j])
        for j in range(i+1, n):
            L = L * (X - Px[j]) / (Px[i] - Px[j])
        P = P + L * Py[i]
        print(Py[i])
        
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlim(-3,9)
    ax.set_ylim(-5, 7)
    ax.plot(X, P, color = 'orange')
    ax.plot(Px, Py, '+', color = 'dimgrey')


# Algorithme de Castelnau
    
def Castelnau(Px, Py, t):

    j = 0
    
    while len(Px) != 1:
        plt.plot(Px, Py, '+', color = couleurs[j])
        plt.plot(Px, Py, color = couleurs[j])
        print(Px, Py)
        
        
        Mx, My = [], []
        for i in range(0, len(Px)-1):
            Mx.append((1-t)*Px[i] + t*Px[i+1])
            My.append((1-t)*Py[i] + t*Py[i+1])
        Px = Mx.copy()
        Py = My.copy()
        
        j += 1
        
    plt.plot(Px, Py, 'o', color = 'k')
    print(Px, Py)
        
        

# Courbe de Bézier
    
def Bezier(Px, Py):
    t = np.linspace(0, 1)
    
    X = np.zeros(len(t))
    Y = np.zeros(len(t))
    
    n = len(Px) - 1
    
    for i in range(n+1):
        X = X + special.binom(n, i)*(1-t)**i*t**(n-i)*Px[i]
        Y = Y + special.binom(n, i)*(1-t)**i*t**(n-i)*Py[i]
        
        
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlim(-3,9)
    ax.set_ylim(-5, 7)
    
    ax.plot(Px[0:2],Py[0:2], '-', color = 'darkgrey')
    ax.plot(Px[2:4],Py[2:4], '-', color = 'darkgrey')
    ax.plot(Px, Py, '.', color = 'dimgrey')
    ax.plot(X,Y, color = 'darkorchid')
    #plt.plot(Px, Py)
    
def Bezier_coeff(n):
    t = np.linspace(0, 1)
    for i in range(n+1):
        plt.plot(t, special.binom(n, i)*(1-t)**i*t**(n-i), color = couleurs[i])
        
        
    #ne passe pas par les points initiaux
def B_splines_coeff(n, k): #affiche les courbes des coefficients
    C = list(range(0,n+k)) #les coefficients
    
    T = np.linspace (0,n+k,NB_POINTS)
    
    #bi, 1 #initialisation
    f0 = lambda i : lambda t : 1 if C[i]<=t<C[i+1] else 0
    S0 = [[f0(i)(t) for t in T] for i in range (0,n+k-1)]
    
    
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlim(0,n+k-1)
    
    
    for j in range (0,n+k-1):
        ax.plot(T, S0[j], couleurs[0])
    plt.show()
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlim(0,n+k-1)
    
    S = [S0]
    T = np.array([t for t in T])
    for p in range(1,k):
        Sp = [((T - C[i])*S[p-1][i])/ (C[i+p]-C[i]) + ((C[i+p+1] - T)*S[p-1][i+1]) / (C[i+p+1]-C[i+1]) 
             for i in range (0,n+k-(p+1))]
        
        for i in range (0,len(Sp)):
            ax.plot(T, Sp[i], color = couleurs[p+1])
            
        plt.show()
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        ax.set_xlim(0,n+k-1)
        S.append(Sp)


def B_splines_coeff_fonctionnel(n, p):  #renvoit les coefficients
    
    C = [0]*(p-1)+list(range(0,n-p+2))+[n-p+1]*(p-1)
    # C = [0,0,0,1,2,3,4,4,4] #les coefficients pour 6
    
    T = np.linspace (0,max(C),NB_POINTS)
    
    #bi, 1 #initialisation
    f0 = lambda i : lambda t : 1 if C[i]<=t<=C[i+1] else 0
    S0 = [[f0(i)(t) for t in T] for i in range (0,len(C)-1)]
    
    
    S = [S0]

    
    def aux (x,y): #permet de contrer le cas ou deux noeuds ont egaux
        if x == y:
            return 0
        return 1/(x-y)
    
    for k in range(1,p):
        Sk = [(T - C[i])*S[k-1][i] * aux(C[i+k],C[i]) + (C[i+k+1] - T)*S[k-1][i+1] * aux(C[i+k+1],C[i+1]) 
              for i in range (0,len(C)-(k+1))]
                
        S.append(Sk)
    return Sk

def B_splines_coeff_fonctionnel2(n, k):  #renvoit les coefficients
    
    C = list(range(0,n+k)) #les coefficients
    
    T = np.linspace (0,n+k,NB_POINTS)
    
    #bi, 1 #initialisation
    f0 = lambda i : lambda t : 1 if C[i]<=t<C[i+1] else 0
    S0 = [[f0(i)(t) for t in T] for i in range (0,n+k-1)]

    
    S = [S0]
    for p in range(1,k):
        Sp = [((T - C[i])*S[p-1][i])/ (C[i+p]-C[i]) + ((C[i+p+1] - T)*S[p-1][i+1]) / (C[i+p+1]-C[i+1]) 
             for i in range (0,n+k-(p+1))]
    
        S.append(Sp)
    return(T, Sp)
    
def B_splines(Px, Py, p, c):

    n = len(Px)
    X = np.zeros(NB_POINTS)
    Y = np.zeros(NB_POINTS)
    
    T, S = B_splines_coeff_fonctionnel2(n,p)

    
    for i in range(0, len(S)):
        Nx = Px[i]*S[i] #on multiplie le tableau numpy par le coefficient
        Ny = Py[i]*S[i]
        X = X+Nx
        Y = Y+Ny
    
    xmin = -1
    xmax = -1
    for k in range(len(T)):
        if T[k]>p-1 and xmin == -1:
            xmin = k
        elif T[k]>n and xmax == -1:
            xmax = k
    
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlim(-3,9)
    ax.set_ylim(-5, 7)
    ax.plot(X[xmin:xmax],Y[xmin:xmax], '+', color = couleurs[c])
    ax.plot(Px,Py, '--', color = 'grey')
    ax.plot(Px, Py, '+', color = 'dimgrey')
        

        
            
    