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


import random
random.seed()

#pgcd
def gcd(a, b):
    while(b!=0):
    	 a, b = b, a%b
    return abs(a)

def euclide_extended(a, b, verbose=False):
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
        if verbose:
            print("{0}\t{1}\t{2}\t{3}".format(q, r0, u0, v0))
        if r1 == 0:
            fini = True
            if verbose:
                print("Algo termine")
                print("pgcd({0},{1}) = {2} = {3}x{4} + {5}x{6}".format(a, b, r0, u0, a, v0, b))
    return r0, u0, v0

#print(euclide_extended(255, 45))



#1/a mod p
def inversion_modulaire(a, p):
    pgcd, u, v = euclide_extended(a%p, p, False)
    if pgcd != 1:
        print("{0} non inversible modulo {1}".format(a, p))
        return 0
    return (v+p)%p

#print inversion_modulaire(5,23)



def inversibles_ZnZ(n):
    liste = []
    for i in xrange(1, n):
        pgcd, u, v = euclide_extended(i, n, False)
        if pgcd == 1:
            liste.append(i)
    return liste

#print(inversibles_ZnZ(13))
#print(inversibles_ZnZ(18))

#reconstruit modulo N (produit des ni) x = ai mod ni où les La est la liste des ai et Ln celle des ni.
#pré-requis : La et Ln de même taille
def crt(La, Ln, N=0):
	#si N n'est pas fourni
	if(N==0):
		N = 1
		for ni in Ln :
			N *= ni
	x = 0	
	for i in range(len(Ln)):
		p = N // Ln[i]
		v = inversion_modulaire(p, Ln[i])
		x = (x + La[i]*v*p) % N
	return x
	
#print(crt([2,3,2],[3,5,7],105))

def ind_euler(n):
	phi = 0
	for i in range(1, n):
		pgcd, u, v = euclide_extended(i, n, False)
		if pgcd == 1:
			phi += 1
	return phi


# a^n mod p (fastexp)
def exp_mod(a, n, p):
    r1 = 1
    i = 0
    while i < n:
        r = a
        j = 1
        while 2 * j <= n - i:
            r = (r * r) % p
            j = 2 * j
        i = i + j
        r1 = (r1 * r) % p
    return r1
    
#racine carrée entière
def isqrt(n):
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

def ilog(x, b): # plus grand entier l tel que b**l <= x.
    r = x
    l = 0
    while r >= b:
        r = r // b
        l += 1
    return l

def temoin_miller(a,n):
    #n>2, a>1
    r=0
    s=0
    d = n-1
    while(r==0):
        r = d % 2
        if(r==0):
            d = d//2
            s += 1

    x = exp_mod(a,d,n)
    if (x==1 or x==n-1):
        return False

    while(s>1):
        x = x*x % n
        if(x == n-1):
            return False
        s = s-1

    return True
    
    
def rabin_miller(n,k=1,verbose=False):
	if(n<4):
		return True
	for i in range(k):
		a = random.randint(2,n-2)
		if(temoin_miller(a,n)):
			if(verbose): print("Rabin-Miller :",n,"is not prime")
			return False
			
	if(verbose): print("rabin_miller :", n, "\nhas", 100*(1-1/(2**k)), "% chances to be prime")
	return True

def mersenne(p): #Mersenne numbers (not all primes)
    return 2**p - 1
    
#rabin_miller(mersenne(107),10,True)

#retourne la liste des nombres premiers de [a,b[ testés avec k tours de Rabin-Miller
def get_primes(a,b,k,verbose=False):
	primes = []
	if(a<5):
		start = 5
		primes = [2,3]
	else:
		start = a
	for n in range(start,b):
		#rabin_miller ne fait les k tours que si nécessaire
		if(rabin_miller(n,k)):
			primes.append(n)
	if(verbose): print("get_primes : those numbers have", 100*(1-1/(2**k)), "% chances to be prime")
	return primes
	
#print(get_primes(10000,10100,20,True))
		
#génère un nombre premier entre a et b, testé avec k tours de Rabin-Miller
def rand_prime(a,b,k,verbose=False):
	r = random.randint(a,b)
	while(not rabin_miller(r,k)):
		r = random.randint(a,b)
	if(verbose): print("rand_prime :", r, "\nhas", 100*(1-1/(2**k)), "% chances to be prime")
	return r
	

#retourne 0 si aucun generateur n'est trouve : p n'est pas premier
def trouver_generateur(p):
	for a in range(2, p):
		r = a
		i = 1
		while r != 1:
			i += 1
			r = (r * a) % p
		if i == p-1:
			return a
	return 0
    
#print trouver_generateur(13)
    
#trouve un générateur du groupe cyclique Z/pZ* par random
def rand_gen(p):
 	i = 1
 	while(i < p-1):
 		a = random.randint(2,p-1)
 		r = a
 		i = 1
 		while(r != 1):
 			i += 1
 			r = (r*a)%p
 	return a
 	
def get_group(n, verbose=False):
	if(verbose): print("\n====== Group generator ======")
	if(verbose): print("prime seed :",n)
	k = 2
	p = k*n + 1
	while(rabin_miller(p,10)==False):
		k += 1
		p = k*n + 1
	g = 1
	if(verbose): print(p,"=",k,"x",n,"+ 1")
	while(g==1):
		h = random.randint(2,p-1)
		g = exp_mod(h,k,p)
	return g,p
	
#n = rand_prime(1000000,2000000,10)
#g1 = rand_gen(n)
#g2 = trouver_generateur(n)
#print("generateurs pour N =",n,":",g1,g2)
	
