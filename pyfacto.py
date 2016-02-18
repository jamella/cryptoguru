#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016 Amaury Behague <amaury.behague@gmail.com>
#
# This file is part of cryptoguru.
#
# cryptopwn is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

__author__ = "Amaury Behague"
__copyright__ = "Copyright 2016, Amaury Behague"
__license__ = "GPL"
__version__ = "3"
__email__ = "amaury.behague@gmail.com"
__status__ = "Beta"

"""
.. module:: pyfacto
	:platform: Unix
	:synopsis: Useful factoring tools.
	
.. moduleauthor:: Amaury Behague <amaury.behague@gmail.com>
"""

import itools, random, time
from multiprocessing import Process, Queue

random.seed()

#Retourne un facteur de n
#pré-requis : n impair
def facto_fermat(n, verbose=False):
    if(n%2 == 0):
        print("Erreur factoFermat : n pair")
        return 0
    aux = itools.isqrt(n)
    if aux*aux == n:
        if(verbose): print("n est un carré")
    aux += 1
    res = aux*aux - n
    if(verbose): print(res)
    r = itools.isqrt(res)
    e = r*r - res
    while(e != 0):
        aux += 1
        res = aux*aux - n
        r = itools.isqrt(res)
        e = r*r - res
    res = aux + itools.isqrt(res)
    if(vernose): print("factoFermat :\n",n,"=",res,"*",n//res)
    return res

#n = itools.mersenne(17) * itools.mersenne(19)
#itools.rabin_miller(n,10)
#print(facto_fermat(n,True))

#retourne un petit facteur de n
#pré-requis : n composé
def rho_pollard(n,verbose=False):
	x = random.randint(1,n)
	y = x
	g = 1
	while(g==1):
		x = (x*x + 1)%n
		y = (y*y + 1)%n
		y = (y*y + 1)%n
		g = itools.gcd(n,x-y)
	if(verbose): print("rho_pollard :\n",n,"=",g,"x",n//g)
	return g
	
#n = itools.rand_prime(2**35,2**36,20)*itools.rand_prime(2**35,2**36,20)
#rho_pollard(n,True)


#Version améliorée par Brent : http://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf
def rho_pollard_brent(n,verbose=False):
	x = random.randint(1,n)
	c = random.randint(1,n)
	m = random.randint(1,n)
	y = x
	g = 1
	r = 1
	q = 1
	while(g==1):
		x = y
		for i in range(r):
			y = (y*y + c)%n
		k = 0
		while( k<r and g==1):
			ys = y
			for i in range(min(m,r-k)):
				y = (y*y + c)%n
				q = q*(x-y) % n
			g = itools.gcd(q,n)
			k += m
			
		r *= 2
	
	if(g==n):
		g = 1
		while(g==1):
			ys = (ys*ys + c)%n
			g = itools.gcd(x-ys,n)
	
	if(verbose): print("rho_pollard_brent :\n",n,"=",g,"x",n//g)
	return g

#n = itools.rand_prime(2**48,2**49,20)*itools.rand_prime(2**48,2**49,20)
#rho_pollard_brent(n,True)

#version parallélisée
#d'intérêt limité, permet juste d'éviter de souffrir d'un mauvais tirage
#augmente quand même les chances de vite tomber sur une relation puisque l'on a plusieurs points de départ (racine du nombre de procs)
def rho_pollard_brent_p(n,jobs=8,verbose=False):
	def core(n, output_queue): output_queue.put(rho_pollard_brent(n,False))
	queue = Queue()
	procs = []
	for j in range(jobs):
		procs.append(Process(target=core, args=(n, queue)))
	for p in procs:
		p.start()
	g = queue.get()
	for p in procs:
		p.terminate()
	if(verbose): print("rho_pollard_brent_p :\n",n,"=",g,"x",n//g)
	return g
	
#rho_pollard_brent_p(n,8,True)

#retourne un facteur de n
#pré-requis : n possède un facteur p tel que p-1 soit B-lisse
#nbRM : nombre de tours de Rabin-Miller pour générer la liste des premiers
#B = 1000000 permet de trouver 1/4 des facteurs de 10 chiffres, et 1/27 de ceux de 18 chiffres
def pm1_pollard(n,B,nbRM=20,verbose=False):
	lB = itools.ilog(B,2)
	primes = itools.get_primes(1,B,nbRM)
	a = random.randint(1,n)
	for q in primes:
		e = lB//itools.ilog(q,2)
		a = itools.exp_mod(a,q**e,n) #q**e ~ B
	g = itools.gcd(a-1,n)
	
	if g>1 and g<n:
		if(verbose): print("pm1_pollard :\n",n,"=",g,"x",n//g)
		return g
	else:
		if(verbose): print("pm1_pollard failed with", B, ":(")
		return 0
		

#n = itools.rand_prime(2**45,2**46,20)*itools.rand_prime(2**45,2**46,20)
#pm1_pollard(n,10000,10,True)

def pm1_pollard_auto(n,Bmax,verbose=False):
	#paramètres ajustés pour une complexité raisonnable
	B = 16000
	nbRM = 10
	alternate = False
	g = 1
	while(g==1 and B<=Bmax):
		lB = itools.ilog(B,2)
		primes = itools.get_primes(1,B,nbRM)
		g = n
		while(g==n):
			a = random.randint(1,n)
			for q in primes:
				e = lB//itools.ilog(q,2)
				a = itools.exp_mod(a,q**e,n) #q**e ~ B
			g = itools.gcd(a-1,n)
			if(g==n and verbose):
				print("N reached, resuming...")
		if(g==1):
			c = a
			primes = itools.get_primes(B,2*B,nbRM)
			d = primes[0]
			a = itools.exp_mod(c,d,n)
			for i in range(1,len(primes)):
				g = itools.gcd(a-1,n)
				if(g!=1): 
					alternate=True
					break
				d = primes[i]-primes[i-1]
				a *= itools.exp_mod(c,d,n) #peut être amélioré en mémorisant les c**(2k) dans une table
			if(g==1):
				B *= 4
	if(g!=1 and g!=n):
		if(verbose): 
			if(alternate): print("pm1_pollard_auto ( B' =",2*B,") :\n",n,"=",g,"x",n//g)
			else: print("pm1_pollard_auto ( B =",B,") :\n",n,"=",g,"x",n//g)
		return g
	else:
		if(verbose): print("pm1_pollard_auto failed with", Bmax, ":(")
		return 0
	
#pm1_pollard_auto(n,256000,True)

#calcule V[mq] % n quand v = V[m]
#V est la suite de Lucas définie par V[i] = A*V[i-1] - V[i-2]
#on utilise les formules suivantes :
#V[2n] = V[n]*V[n] - 2
#V[m+n] = V[m]*V[n] - V[m-n]
#on calcule à la fois V[qm] et V[(q+1)m] (pour résoudre la dépendance de l'addition)
def lucas_mul(v,q,n):
	x = v
	y = (v*v - 2) % n
	for bit in bin(q)[3:]:
		if(bit=='1'):
			x=(x*y-v) % n
			y=(y*y-2) % n
		else:
			y=(x*y-v) % n
			x=(x*x-2) % n
	return x

#retourne un facteur de n
#pré-requis : n possède un facteur p tel que p+1 soit B-lisse
#nbRM : nombre de tours de Rabin-Miller pour générer la liste des premiers
def pp1_williams(n,B,nbRM=20,verbose=False):
	lB = itools.ilog(B,2)
	primes = itools.get_primes(3,B,nbRM)
	a = random.randint(1,n)
	v = a
	i=0
	for q in primes:
		e = lB//itools.ilog(q,2)
		#plus efficace que :  for i in range(e): v = lucal_mul(v,q,n)
		#évite de calculer e fois la décomposition binaire de q et d'appeler e fonctions, au prix d'une exponentiation entière
		v = lucas_mul(v,q**e,n) #q**e ~ B
		g = itools.gcd(v-2,n)
		if(g!=1 and g!=n): break
				
	if(g>1 and g<n):
		if(verbose): print("pp1_williams :\n",n,"=",g,"x",n//g)
		return g
	else:
		if(verbose): print("pp1_williams failed with", B, ":(")
		return 0
		
def pp1_williams_auto(n,Bmax,verbose=False):
	B = 1000
	nbRM = 10
	alternate = False
	g = 1
	while(g==1 and B<=Bmax):
		lB = itools.ilog(B,2)
		B2 = (10 ** (itools.ilog(B,10)//5))*B #tiré de resultats expérimentaux : http://www.loria.fr/~zimmerma/records/Pplus1.html
		primes = itools.get_primes(3,B,nbRM)
		g = n
		
		#Step 1
		while(g==n):
			a = random.randint(1,n)
			v = a
			for q in primes:
				e = lB//itools.ilog(q,2)
				#plus efficace que :  for i in range(e): v = lucal_mul(v,q,n)
				#évite de calculer e fois la décomposition binaire de q et d'appeler e fonctions, au prix d'une exponentiation entière
				v = lucas_mul(v,q**e,n) #q**e ~ B
			#les valeurs atteintes par g sont croissantes et stagent à n lorsqu'atteints
			#=> pas besoin de les évaluer à chaque multiplication de l'indice
			#si g == n -> on a dépassé les facteurs, il faut recommencer avec la même borne (voire plus petite ?)
			#le tirage de a sera déterminant
			#si g==1 la borne n'est pas suffisante, on passe au Step 2 puis on augmente
			g = itools.gcd(v-2,n)
			if(g==n and verbose):
				print("N reached, resuming...")
		
		#Step 2 TODO : optimal? L'idée des multiples de 6 peut-elle être reprise pour p-1 Pollard ?
		if(g==1):
			c = (B//6)*6
			v6 = lucas_mul(v,6,n)
			Vl = []
			Vl.append(lucas_mul(v,c,n))
			Vl.append(lucas_mul(v,c+6,n))
			i = 2
			c = c+12
			while(c<B2):
				Vl.append( (Vl[i-1]*v6 - Vl[i-2]) % n )
				c += 6
				i += 1
			produit = 1
			for vi in Vl:
				produit = ( produit * (vi-v) ) % n
			g = itools.gcd(produit,n)
			if(g==1):
				B *= 10
			else:
				alternate = True
			
	if(g!=1 and g!=n):
		if(verbose): 
			if(alternate): print("pp1_williams_auto ( B' =",B2,") :\n",n,"=",g,"x",n//g)
			else: print("pp1_williams_auto ( B =",B,") :\n",n,"=",g,"x",n//g)
		return g
	else:
		if(verbose): print("pp1_williams_auto failed with", Bmax, ":(")
		return 0

#n = itools.rand_prime(2**25,2**36,20)*itools.rand_prime(2**25,2**36,20)
#pm1_pollard_auto(n,1024000,True)
#pp1_williams_auto(n,1000000,True)

def edwards_add(P1,P2,Edwards):
	x1,y1 = P1
	x2,y2 = P2
	a,d,n = Edwards
	x1y1 = (x1*y1)%n
	x2y2 = (x2*y2)%n
	p = ( x1y1 + x2y2 ) % n
	q = ( (y1*y2 % n) + (a*x1*x2 % n) ) % n
	x = p//q
	y=-1
	p = ( x1y1 - x2y2 + n) % n
	q = ( (x1*y2 % n) - (y1*x2 % n) + n) % n
	if(q==0): print("edwards_add ERROR : division par 0")
	else: y = p//q
	return (x,y)

def edwards_double(P1,Edwards):
	x,y = P
	a,d,n = Edwards
	xy = (x*y)%n
	x2 = (x*x)%n
	y2 = (y*y)%n
	return ( ((2*xy)%n) // ((1 + d*((xy*xy)%n))%n) ), ( ((y2 - a*x2)%n) // ((n+1 - d*((xy*xy)%n))%n) )  
	
#retourne un facteur de n
#algorithme le plus efficace pour trouver des facteurs de moins de 25 chiffres
#version avec courbes de d'Edwards (plus rapide que les courbes de Montgomery)
#E(a,d) : ax² + y² = 1 + dx²y²  (a et d non nuls)
#On considère le cas où a = 1 et d != 0,1
#TODO
def eecm(n):
	u = random.randint(2,n-1)
	u2 = (u*u)%n
	u2p1 = (u*u+1)%n
	x = ((u2-1)%n)//u2p1
	y = (((u-1)*(u-1))%n)//u2p1
	d1 = ( ((u2p1*u2p1)%n)*(u2p1) ) % n
	d1 = ( d1*(u2-4*u+1) ) % n
	d2 = ( (u-1)*(u-1) ) % n
	d2 = ( ((d2*d2)%n)*d2 ) % n
	d2 = (d2 * ((u+1)*(u+1))%n ) % n
	d = d1 // d2
	

	
#TODO : version avec Courbes d'Edward


#TODO : Multiple Polynomial Quadratic Sieve (MPQS) => cf COURS CRYPTO (GNFS trop compliqué)
#MPQS : méthode la plus rapide pour trouver des facteurs entre 25 et 100 chiffres

#TODO : fonction de facorisation appelant d'abord p-1,p+1, puis rho & ECM puis MPQS

def dico_add_to_first(d1,d2):
	for p in d2.keys():
		if p in d1.keys():
			d1[p] += d2[p]
		else:
			d1[p] = d2[p]
			

def easy_facto(n,verbose=False):
	dico = {}
	p = n
	while(itools.rabin_miller(p,10) == False):
		q = rho_pollard_brent_p(p,8)
		if(q != 1):
			if(itools.rabin_miller(q,10) == True):
				if(q in dico.keys()):
					dico[q] += 1
				else:
					dico[q] = 1
			else:
				dico_add_to_first(dico, easy_facto(q))
			p = p//q

	if(p != 1):
		if(p in dico.keys()):
			dico[p] += 1
		else:
			dico[p] = 1
	if(verbose): print("easy_facto :\nN =",n,"\nfactorisation :", dico)
	return dico

#print(easy_facto(344))	
#print(easy_facto(random.randint(1,100000000)))

