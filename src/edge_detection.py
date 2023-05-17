# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 16:41:45 2020
****************** Détection de contours *****************
@author: Gulliver
"""

import math as m
import random as rd
from PIL import Image, ImageFilter, ImageOps, ImageDraw
import matplotlib.pyplot as plt
import numpy as np


irm3 = Image.open('Images/irm 3.png')
irm4 = Image.open('Images/irm 4.jpg')
irm5 = Image.open('Images/irm 5.jpg')

# Détècte un contour horizontal pour une image PNG !

def contours(nom_image):
    (liste_bandes, largeur, hauteur) = decoupage_et_pretraitement(nom_image)
    # On récupère une version simplifiée de l'information issue de l'image,
    # Sous forme d'une liste de listes qui représentent chacune une tranche verticale,
    # de 8 pixels de large et ses contrastes
    X = []
    Y = []
    abscisse = 0
    
    for liste in liste_bandes:
        # On ajoute à X l'abscisse de la tranche nommée "liste" correspondante
        # Et on ajoute à Y la position du minimum de cette liste:
        X.append(abscisse)
        Y.append(recuit2(liste))
        abscisse += 1
    TraceX,TraceY = trace2(X,Y)
    
    fusion_images(TraceX,TraceY,nom_image)






def decoupage_et_pretraitement(nom):
    #On applique un filtre matriciel à l'image qui fait ressortir les discontinuités,
    #Puis on la convertit en niveaux de gris
    img = ImageOps.grayscale(filtre_prewitt(nom))
    
    largeur,hauteur = img.size
    tabpixel = img.load()
    liste_bande = []
    # On établit ensuite une liste qui contient chaque "tranche" d'image, chaque tranche étant
    # une liste contenant la valeur moyenne sur 8 pixels de large 
    for i in range(0,int(largeur/8)):
        liste = []
        for j in range(1,hauteur-1):
            total = 0
            for k in range(0,8):
                total += tabpixel[i*8+k,j]
            liste.append(total / 8)
        liste_bande.append(liste)
    return(liste_bande, largeur, hauteur)

def filtre_prewitt(nom):
    img = Image.open(nom)
    #On applique un filtre matriciel à l'image qui fait ressortir les discontinuités,
    #Puis on la convertit en niveaux de gris
    img = ImageOps.grayscale(img)
    largeur,hauteur = img.size
    Result = Image.new("RGB",(largeur - 2, hauteur - 2), "white")

    for i in range(1,hauteur-1):
        for j in range(1,largeur-1):
            p1 = img.getpixel((j-1,i))
            p2 = img.getpixel((j,i-1))
            p3 = img.getpixel((j+1,i))
            p4 = img.getpixel((j,i+1))
            n = 255 - int(Norme(p1,p2,p3,p4))
            Result.putpixel((j-1,i-1),(n,n,n))
    Result.show()
    return(Result)


def Norme(p1,p2,p3,p4):
    n = m.sqrt((p1-p3)*(p1-p3) + (p2-p4)*(p2-p4))    
    return n


def b_splines_coeff2(n): #cubic #nombre de points
    
    C = [0,0]+list(range(0,n-1))+[n-2,n-2]
   # C = [0,0,0,1,2,3,4,4,4] #les coefficients pour 6
    
    T = np.linspace (0,max(C),200)
    
    #bi, 1 #initialisation
    f0 = lambda i : lambda t : 1 if C[i]<=t<=C[i+1] else 0
    S0 = [[f0(i)(t) for t in T] for i in range (0,len(C)-1)]
    
    
    S = [S0]
    T = np.array([t for t in T])
    
    def aux (x,y): #permet de contrer le cas ou deux noeuds ont egaux
        if x == y:
            return 0
        return 1/(x-y)
    
    for p in range(1,3):
        Sp = [(T - C[i])*S[p-1][i] * aux(C[i+p],C[i]) + (C[i+p+1] - T)*S[p-1][i+1] * aux(C[i+p+1],C[i+1]) 
             for i in range (0,len(C)-(p+1))]
                
        S.append(Sp)
    return Sp


def trace2(Px, Py):
    
    plt.plot(Px, Py,'.')
    
    X = np.zeros(200)
    Y = np.zeros(200)
    
    S = b_splines_coeff2(len(Px))
    
    for i in range(0, len(S)):
        Nx = Px[i]*S[i] #on multiplie le tableau numpy par le coefficient
        Ny = Py[i]*S[i]
        X += Nx
        Y += Ny
    plt.plot(X,Y)
    return(X,Y)



def recuit2(liste):
    hauteur = len(liste)
    mini,maxi = 0,hauteur -1
    T0,tseuil,tau = 2,1000,10000
    # 
    def f(x):
        pos = int(x)
        return(liste[pos]*(x - pos) + liste[pos+1]*(pos + 1 - x))
    
    globalcursor = (mini + maxi)/2
    globalmin = f((mini + maxi)/3)
    lastf = f((mini + maxi)/2)
    #Initialisations de la température et du temps:
    t = 0
    T = T0
    while t < tseuil:
        #On définit ici une fonction qui rapproche les points aléatoires
        # sur le minimum conjecturé aux itérations précedentes
        if t%25 == 0:
            expo = 1+(1-m.exp(-t/1000))*t/300
            def densite(decoup_segm,globalcursor):
                param = rd.uniform(0,1)
                if param > decoup_segm:
                    return((maxi - globalcursor)*((param - decoup_segm)*(1/(1- decoup_segm)))**expo + globalcursor)
                else:
                    return(-(globalcursor - mini)*((decoup_segm - param)*(1/decoup_segm))**expo + globalcursor)

        decoup_segm = abs((globalcursor- mini)/(maxi - mini))  #On découpe le segment [0,1] 
        #pour donner une densité de probabilité d'apparition du point autour du dernier minimum enregistré
        newcursor = densite(decoup_segm,globalcursor)
        f_new = f(newcursor)
        if f_new < lastf:
            lastf = f_new
        elif rd.uniform(0,1) < m.exp(-1/T):
            lastf = f_new
        if f_new < globalmin:
            globalmin = f_new
            globalcursor = newcursor
        t += 1
        T = T0*m.exp(-t/tau)
    return(hauteur - globalcursor)


def fusion_images(TraceX,TraceY,image):
    img = Image.open(image)
    width, height = img.size
    RED = (255, 52, 19)
    draw = ImageDraw.Draw(img)
  
    for i in range(0,len(TraceX)-1):
        draw.line((8*TraceX[i], height - TraceY[i],8*TraceX[i+1], height - TraceY[i+1]), RED)
    img.show()
    img.save("Détection_contours.png","PNG")












