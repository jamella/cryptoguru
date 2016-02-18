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
.. module:: pwndlp
	:platform: Unix
	:synopsis: Efficient Attacks on the DLP.
	
.. moduleauthor:: Amaury Behague <amaury.behague@gmail.com>
"""


import itools, random, pyfacto, time
from multiprocessing import Process, Queue

random.seed()


def rho_pollard_dlp(g, h, p, n, verbose=False):
	"""Pollard's Rho applied to DLP.
	
	Args:
		- *g (int)*: a generator
		- *h (int)*: an integer in <g>
		- *p (int)*: the order of <g>, must be prime (else call Pohlig-Hellman first).
		- *n (int)*: the modulo
		
	Optional args:
		- *verbose (bool)*: set to True if you want a display.
		
	Returns:
		- *(int)*: x such that [x]g = h mod n
	
	"""
	
	#Bruteforce to start
	if(verbose): print("Brute forcing..")
	x = g
	for i in range(1,min(p+1,100)):
		x = (x*g)%n
		if(x==h):
			if(verbose): print("rho_pollard_dlp :\n[",i+1,"]",g,"=",h)
			return i+1
	
	#then standard rho
	S0 = n//3
	S1 = 2*(n//3)
	x,y = 0,0
	while(y==0):
		if(verbose): print("starting new walk..")
		gx, hx = random.randint(0,p), random.randint(0,p)
		x = ( itools.exp_mod(g,gx,n)*itools.exp_mod(h,hx,n) ) % n
		gy, hy = gx, hx
		y = x
	
		loop = True
	
		while(loop):
			#x[i]
			if(x <= S0):
				x = (x*h)%n
				hx = (hx+1)%p
			elif(x <= S1):
				x = (x*x)%n
				gx, hx = (2*gx)%p, (2*hx)%p
			else:
				x = (x*g)%n
				gx = (gx+1)%p
		
			#x[2i]
			for _ in range(2):
				if(y <= S0):
					y = (y*h)%n
					hy = (hy+1)%p
				elif(y <= S1):
					y = (y*y)%n
					gy, hy = (2*gy)%p, (2*hy)%p
				else:
					y = (y*g)%n
					gy = (gy+1)%p
				
			if(x==y): loop = False
	
		x = (gy - gx +p)%p
		y = (hx - hy +p)%p
	
	y = itools.inversion_modulaire(y,p)
	if(verbose): print("rho_pollard_dlp :\n[",(x*y)%p,"]",g,"=",h)
	return (x*y)%p


#TODO: use mixed walks (instead of only adding walks).
def rho_pollard_dlp_adv(g, h, p, n, b, k, verbose=False):
	"""Improved version of Pollard's Rho applied to DLP. It uses k-adding walks.
	
	Args:
		- *g (int)*: a generator
		- *h (int)*: an integer in <g>
		- *p (int)*: the order of <g>, must be prime (else call Pohlig-Hellman first)
		- *n (int)*: the modulo
		- *b (int)*: exponent bound used to generate random walks (p seems to be the best)
		- *k (int)*: number of partitions (20 or more is advised)
		
	Optional args:
		- *verbose (bool)*: set to True if you want a display.
		
	Returns:
		- *(int)*: x such that [x]g = h mod n
	
	"""
	#Bruteforce to start
	x = g
	if(verbose): print("Brute Forcing..")
	for i in range(1,min(p+1,100)):
		x = (x*g)%n
		if(x==h):
			if(verbose): print("rho_pollard_dlp_adv :\n[",i+1,"]",g,"=",h)
			return i+1
			
	if(verbose): print("computing f..")
	puissances = []
	coefficients = []
	for s in range(k):
		ms,ns = 0,0
		while(ms==0 and ns==0):
			ms,ns = random.randint(0,b), random.randint(0,b)
		puissances.append((ms,ns))
		coefficients.append( (itools.exp_mod(g,ms,n) * itools.exp_mod(h,ns,n)) % n)
	
	x,y = 0,0
	while(y==0):
		if(verbose): print("starting new walk..")
		gx, hx = random.randint(0,p-1), random.randint(0,p-1)
		x = ( itools.exp_mod(g,gx,n)*itools.exp_mod(h,hx,n) ) % n
		gy, hy = gx, hx
		y = x
	
		loop = True
	
		while(loop):
			#x[i]
			i = x % k
			ms,ns = puissances[i]
			coef = coefficients[i]			
			x = (x*coef) % n
			gx, hx = (gx+ms)%p, (hx+ns)%p
		
			#x[2i]
			for _ in range(2):
				i = y % k
				ms,ns = puissances[i]			
				coef = coefficients[i]			
				y = (y*coef) % n
				gy, hy = (gy+ms)%p, (hy+ns)%p
				
			if(x==y): loop = False
	
		x = (gy - gx +p)%p
		y = (hx - hy +p)%p

	y = itools.inversion_modulaire(y,p)
	if(verbose): print("rho_pollard_dlp_adv :\n[",(x*y)%p,"]",g,"=",h)
	return (x*y)%p
	
#sub-function for parallelized Pollard's Rho
#TODO: set distinguished points criteria dynamically.
def sub_rho(g, h, p, n, coefficients, puissances, k, dico):
	
	if(itools.ilog(n,2) > 40):
		critere = 2**(itools.ilog(n,2)-14)
	elif(itools.ilog(n,2) > 20):
		critere = 2**(itools.ilog(n,2)-10)
	else:
		critere = 2**(itools.ilog(n,2)//2)
	
	points = 0
	y = 0
	while(y==0):
	
		gx, hx = random.randint(0,p-1), random.randint(0,p-1)
		x = ( itools.exp_mod(g,gx,n)*itools.exp_mod(h,hx,n) ) % n
	
		loop = True
	
		while(loop):	
			i = x % k
			ms,ns = puissances[i]
			coef = coefficients[i]			
			x = (x*coef) % n
			gx, hx = (gx+ms)%p, (hx+ns)%p
		
			if(x % critere == x):
				#print("distinguished point hit", points)
				points += 1
				if(x in dico.keys()):
					loop = False
					gy,hy = dico[x]
					x = (gy - gx +p)%p
					y = (hx - hy +p)%p
				else:
					dico[x] = (gx,hx)
		
	y = itools.inversion_modulaire(y,p)
	return (x*y)%p
			
	

def rho_pollard_dlp_par(g, h, p, n, b, k, jobs=8, verbose=False):
	"""Improved parallelized version of Pollard's Rho applied to DLP. Uses distinguished points for optimal efficiency.
	
	Args:
		- *g (int)*: a generator
		- *h (int)*: an integer in <g>
		- *p (int)*: the order of <g>, must be prime (else call Pohlig-Hellman first)
		- *n (int)*: the modulo
		- *b (int)*: exponent bound used to generate random walks (p seems to be the best)
		- *k (int)*: number of partitions (20 or more is advised)
		
	Optional args:
		- *jobs (int)*: number of threads to launch. Should be your number of virtual cores.
		- *verbose (bool)*: set to True if you want a display.
		
	Returns:
		- *(int)*: x such that [x]g = h mod n
	"""
	
	#Bruteforce to start
	x = g
	if(verbose): print("Brute Forcing..")
	for i in range(1,min(p+1,100)):
		x = (x*g)%n
		if(x==h):
			if(verbose): print("rho_pollard_dlp_par :\n[",i+1,"]",g,"=",h)
			return i+1
			
	if(verbose): print("computing f..")
	puissances = []
	coefficients = []
	for s in range(k):
		ms,ns = 0,0
		while(ms==0 and ns==0):
			ms,ns = random.randint(0,b), random.randint(0,b)
		puissances.append((ms,ns))
		coefficients.append( (itools.exp_mod(g,ms,n) * itools.exp_mod(h,ns,n)) % n)
		
	if(verbose): print("launching distributed attack")
	dico = {}
		
	def core(g, h, p, n, coefficients, puissances, k, dico, output_queue): output_queue.put(sub_rho(g, h, p, n, coefficients, puissances, k, dico))
	queue = Queue()
	procs = []
	for j in range(jobs):
		procs.append(Process(target=core, args=(g, h, p, n, coefficients, puissances, k, dico, queue)))
	for p in procs:
		p.start()
	x = queue.get()
	for p in procs:
		p.terminate()
	if(verbose): print("rho_pollard_dlp_par :\n[",x,"]",g,"=",h)
	return x
		

def pohlig_hellman(g, h, n, log_file, verbose=False):
	"""Pohlig-Hellman's algorithm to solve DLP.

	Args:
		- *g (int)*: a generator
		- *h (int)*: an integer in <g>
		- *n (int)*: the modulo
		- *log_file (File)*: an opened file in which results will be written
		
	Optional args:
		- *verbose (bool)*: set to True if you want a display.
		
	Returns:
		- *(int)*: x such that [x mod n-1]g = h mod n
	
	This functions factors n-1 in prime integers in order to be able to call Pollard's Rho on the subgroups
	and then reconstructs the result with the CRT.
	"""
	
	
	p = n-1
	if(verbose):print("\n==== POHLIG_HELLMAN ====")
	facteurs = pyfacto.easy_facto(n-1,verbose)
	La = []
	Ln = []
	for pi in facteurs.keys():
		if(verbose):print("\n=================\nFactor :",pi,"\n=================")
		ei = facteurs[pi]
		xi = 0
		y = h
		b = itools.inversion_modulaire(g,n)
		q = p//pi
		gi = itools.exp_mod(g,q,n)
		for i in range(ei):
			w = itools.exp_mod(y,q,n)
			t0 = time.monotonic()
			xj = rho_pollard_dlp_par(gi, w, pi, n, pi, 50, 8, verbose)
			t1 = time.monotonic()
			log_file.write(str(t1-t0) + " s\n")
			print(t1-t0, "s")
			xi += xj*(pi**i)
			y = (y*itools.exp_mod(b,xj,n)) % n
			b = itools.exp_mod(b,pi,n)
			q = q//pi
		La.append(xi)
		Ln.append(pi**ei)
	x = itools.crt(La,Ln,p)
	if(verbose): print("pohlig_hellman :\n[",x,"]",g,"=",itools.exp_mod(g,x,n))
	return x


