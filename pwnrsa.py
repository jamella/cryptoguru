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
.. module:: pwnrsa
	:platform: Unix
	:synopsis: Efficient attacks on RSA.
	
.. moduleauthor:: Amaury Behague <amaury.behague@gmail.com>
"""


import itools, math
from multiprocessing import Pool

# Generation des listes de convergents de a et b
# pré-requis : a,b rationnels (pour que les listes soient finies)
def gen_convergents(a, b, verbose=False, denom_only=True):
    conv = []
    if a >= b:
        r0 = a
        r1 = b
    else:
        r0 = b
        r1 = a
    u0 = 1
    v0 = 0
    u1 = 0
    v1 = 1
    if verbose:
        print("q\tr\tu\tv")
        print("0\t{0}\t{1}\t{2}".format(r0, u0, v0))
    fini = False
    while not fini:
        q = r0 // r1
        temp = r0 % r1
        r0 = r1
        r1 = temp
        temp = u0 - q * u1
        u0 = u1
        u1 = temp
        temp = v0 - q * v1
        v0 = v1
        v1 = temp
        if denom_only:
            conv.append(abs(v0))
        else:
            conv.append([abs(u0),abs(v0)])
        if verbose:
            print("{0}\t{1}\t{2}\t{3}".format(q, r0, u0, v0))
        if r1 == 0:
            fini = True
            if verbose:
                print("Algo termine")
                print("pgcd({0},{1}) = {2} = {3}x{4} + {5}x{6}".format(a, b, r0, u0, a, v0, b))
    return conv

#Conv = gen_convergents(60728973,160523347,False,False)
#print(Conv)

#Retourne p et q tels que n = pq étant donné n et l'indicatrice d'Euler.
#p et q sont réels ici, ils ne sont entiers que si phi est bien l'indicatrice d'Euler de n.
def get_pq(n,phi,verbose=False):
	a = n - phi +1
	delta = a*a - 4*n
	if(delta<0) :
		if(verbose): print("get_pq : Error! delta < 0")
		return 0,0
	p = (a - math.sqrt(delta))/2
	q = (a + math.sqrt(delta))/2
	return p,q


# n : modulo, e : exposant chiffrement, m : message, c : chiffré (m^e % n) (m et c ne servent qu'à vérifier d)
# d (exposant de déchiffrement) doit vérifier :
# d < (1/3)n^(1/4)
# d est alors le dénominateur d'une fraction réduite de e/n
# ici phi(n) est assimilé à n
def wiener1(n,e,m,c,verbose=False):
    conv = gen_convergents(n,e)
    for d in conv:
    	if(verbose): print("d prob =",d)
    	p = itools.exp_mod(c,d,n)
    	if(verbose): print("pow =",p)
    	if  p == m%n:
            if(verbose): print("\nwiener : success!!!\nSecret exponent :", d, "\n")
            return d
    if(verbose): print("\nwiener failed :(\n")
    return 0
    
#wiener1(160523347, 60728973, 313, 75454098, True)
    
#Wiener sans utiliser de message, on vérifie d en calculant l'indicatrice d'Euler
def wiener2(n,e,verbose=False):
	conv = gen_convergents(n,e,False,False)
	for f in conv:
		k = f[0]
		d = f[1]
		if (k!=0):
			t = (e*d-1) % k
			if(t==0):
				if(verbose):print("d prob :",d)
				phi = (e*d-1)//k
				p,q = get_pq(n,phi,verbose)
				if(verbose): print("p,q prob:", p, q)
				if(int(p)*int(q) == n):
					if(verbose): print("\nwiener : success!!!\nSecret exponent :", d, "\n")
					return d
	if(verbose): print("\nwiener failed :(\n")
    
#wiener2(160523347, 60728973, True)


# n : modulo, e : exposant chiffrement, m : message, c : chiffré (m^e % n) (m et c ne servent qu'à vérifier d)
# d (exposant de déchiffrement) doit vérifier :
# d < (n^(3/4))/|p-q|
# p et q doivent être proches
# d est alors le dénominateur d'une fration réduite de e / (n+1-2*sqrt(n))
# ici phi(n) est assimilé à (n+1-2*sqrt(n)) (car p et q proches de sqrt(n))
def weger1(n,e,m,c,verbose=False):
	conv = gen_convergents(n+1-2*itools.isqrt(n), e)
	for d in conv:
		if(verbose): print("d prob =",d)
		p = itools.exp_mod(c,d,n)
		if(verbose): print("pow =",p)
		if  p == m%n:
			if(verbose): print("\nweger : success!!!\nSecret exponent :", d, "\n")
			return d
	if(verbose): print("\nweger failed :(\n")
	return 0
	
#Weger sans utiliser de message, on vérifie d en calculant l'indicatrice d'Euler
def weger2(n,e,verbose=False):
	conv = gen_convergents(n+1-2*itools.isqrt(n), e, False, False)
	for f in conv:
		k = f[0]
		d = f[1]
		if (k!=0):
			t = (e*d-1) % k
			if(t==0):
				if(verbose):print("d prob :",d)
				phi = (e*d-1)//k
				p,q = get_pq(n,phi,verbose)
				if(verbose): print("p,q prob:", p, q)
				if(int(p)*int(q) == n):
					if(verbose): print("\nweger : success!!!\nSecret exponent :", d, "\n")
					return d
	if(verbose): print("\nweger failed :(\n")
    
#instance résolue par Weger sur laquelle Wiener échoue
"""
p = 10037
q = 10039
n = p*q
phi = (p-1)*(q-1)
d = 509
e = itools.inversion_modulaire(d,phi)
borne_wiener = (1/3) * (n**(0.25))
borne_weger = (n**0.75)/abs(p-q)
print("Borne Wiener :", borne_wiener,"\nBorne Weger :", borne_weger)
weger2(n,e,True)
wiener2(n,e,True)
"""

#Weger étendue, parallélisée
# n : modulo, e : exposant chiffrement
# n = pq avec p et q vérifiant :
# q < p < 2q
# p/q ~ 1 + a/b
# a,b dans [0,B]
# B borne définie par l'utilisateur, la complexité étant proportionelle à B²

#sous-routine de Weger dédiée à la parallélisation
def sub_weger(args):
	n = args[0]
	e = args[1]
	amin = args[2]
	amax = args[3]
	B = args[4]
	for a in range(amin,amax):
		for b in range(1,B):
			x = a/b
			F = n+1 - ((2+x)/math.sqrt(1+x))*math.sqrt(n)		
			conv = gen_convergents(F, e, False, False)
			for f in conv:
				k = int(f[0])
				d = int(f[1])
				if (k!=0):
					t = (e*d-1) % k
					if(t==0):
						phi = (e*d-1)//k
						p,q = get_pq(n,phi)
						if(int(p)*int(q) == n):
							return d
	return 0

def weger_ex(n,e,B,jobs=8,verbose=False):
	map_args = []
	for j in range(jobs):
		amin = int((j/jobs)*B+1)
		amax = int(((j+1)/jobs)*B+1)
		map_args.append([n,e,amin,amax,B])
	
	pool = Pool(processes=jobs)
	Ld = pool.map(sub_weger, map_args)
	for d in Ld:
		if(d!=0):
			if(verbose): print("\nweger extended : success!!!\nSecret exponent :", d, "\n")
			return d
	if(verbose): print("\nweger extended failed :(\n")

	
#instances aléatoires dont la taille fait que Wiener et Weger échouent alors que Weger étendue fonctionne (la plupart du temps)

#q = itools.rand_prime(100000000000,200000000000,20)
#p = itools.rand_prime(100000000000,200000000000,20)
#n = p*q
#print("\nn =",n)
#phi = (p-1)*(q-1)
#d = 50077069
#print("d = n^", math.log(d)/math.log(n))
#print("(p-q)d = n^", math.log(abs(p-q)*d)/math.log(n),"\n\n")
#e = itools.inversion_modulaire(d,phi)

#wiener2(n,e,True)
#weger2(n,e,True)
#weger_ex(n,e,800,8,True)



