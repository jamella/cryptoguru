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
.. module:: itools
	:platform: Unix
	:synopsis: Useful integers functions
	
.. moduleauthor:: Amaury Behague <amaury.behague@gmail.com>
"""

import random
random.seed()


def gcd(a, b):
	"""Computes the greatest common divisor.
	
	Args:
		a (int): some integer
		b (int): some integer
	
	Returns:
		(int) The GCD.
	"""
	
	while(b!=0):
		a, b = b, a%b
	return abs(a)

def euclide_extended(a, b, verbose=False):
	"""Extended version of Euclide's algorithm.
	
	Args:
		a (int): some integer
		b (int): some integer
	
	Optional args:
		verbose (bool): Set to True to get a display
		
	Returns:
		(int,int,int) Three integers, the first one is the gcd, the others are Bezout coefs.
	
	It doesn't matter if a<b, the algorithm switches a and b if necessary.
	However, if a < b, the first coefficient returned will be the one associated to b (gcd,coef_b,coef_a).
	"""
	
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



def inversion_modulaire(a, p):
	"""Inverts a mod p.
	
	Args:
		a (int): some integer
		p (int): some integer, greater than a
		
	Returns:
		(int) 1/a mod p if gcd(a,p) = 1
			  0 if gcd(a,p) > 1
	"""
	
	pgcd, u, v = euclide_extended(a%p, p, False)
	if pgcd != 1:
		print("{0} non inversible modulo {1}".format(a, p))
		return 0
	return (v+p)%p

#print inversion_modulaire(5,23)



#reconstruit modulo N (produit des ni) x = ai mod ni où les La est la liste des ai et Ln celle des ni.
#pré-requis : La et Ln de même taille
def crt(La, Ln, N=0):
	"""Uses the Chinese Remainder Theorem to deduct X = La[i] mod Ln[i] for all i < len(La) = len(Ln)
	X is computed modulo N which is the product of the elements of Ln.
	
	Args:
		La (List): a list of integers
		Ln (List): a list of integers, which size must be len(La).
		
	Optional args:
		N (int): the final modulo. Will be computed if omitted.
		
	Returns:
		(int) the result of the reconstruction modulo N.
		
	For all i, the following statement should be true :
	La[i] < Ln[i]
	
	"""
	
	#if N is omitted
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


def exp_mod(a, n, p):
	"""An efficient function to compute a^n mod p.
	
	Args:
		a (int): an integer
		n (int): an integer
		p (int): an integer > a
		
	Returns:
		int. a^n mod p
	
	Actually just a custom implementation of fastexp. Uses the binary decomposition of n.
	"""
		
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
	
	
def isqrt(n):
	"""Integer Square Root.
	
	Args:
		n (int): an integer
		
	Returns:
		int. The greatest int s such that s*s <= n.
		
	Uses Newton's iterative method.
	"""
	
	
	x = n
	y = (x + 1) // 2
	while y < x:
		x = y
		y = (x + n // x) // 2
	return x

def ilog(x, b):
	"""Integer logarithm in base b.
	
	Args:
		x (int): an integer
		b (int): a base
	
	Returns:
		int. The greatest integer l such that b**l <= x.
		
	"""
	
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
	"""Rabin-Miller primality test.
	
	Args:
		n (int) : the integer which primality you wish to test
	
	Optional Args:
		k (int) : number of repetitions of the test
		verbose (boot) : set to True if you want a display
		
	If this test returns False, you are 100% sure that n isn't prime.
	However, if this test returns True, there's a probability of 1/(2**k) that n isn't prime.
	"""
	
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


def get_primes(a,b,k,verbose=False):
	"""Computes the list of primes between two integers.
	
	Args:
		a (int): an integer
		b (int): an integer such as b > a
		k (int): the number of repetitions Rabin-Miller primality test.
		
	Optional Args:
		verbose (bool): set to True if you want a display.
		
	Returns:
		(List) the list of pseudo-prime integers in [a,b[
	
	Why "pseudo-prime"? Because you can never be sure that they are prime with RM primality test.
	"""
		
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
	"""Generates a random prime number between two integers.
	
	Args:
		a (int): an integer
		b (int): an integer such as b > a
		k (int): the number of repetitions Rabin-Miller primality test.
		
	Optional Args:
		verbose (bool): set to True if you want a display.
		
	Returns:
		(int) a random pseudo-prime between a and b.
		
	Why "pseudo-prime"? Because you can never be sure that they are prime with RM primality test.
	"""
	
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
	"""Generates a "safe" group for the DLP
	
	Args:
		n (int): a prime integer : security modulo / prime seed
		
	Optional Args:
		verbose (bool): set to True if you want a display.
		
	Returns:
		(int,int) g,p with p a prime integer such that p = k*n with some small integer k. And g is a generator of Z/pZ*
		
	"""
	
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
	
