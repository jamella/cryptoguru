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


import itools, random, pyfacto, time
from multiprocessing import Process, Queue

random.seed()

#retourne x tel que [x]g = h mod n
#p : ordre du groupe dont g est le générateur
#pré-requis : p premier, h dans <g>
def rho_pollard_dlp(g, h, p, n, verbose=False):
	
	#brute force pour commencer
	if(verbose): print("Brute forcing..")
	x = g
	for i in range(1,min(p+1,100)):
		x = (x*g)%n
		if(x==h):
			if(verbose): print("rho_pollard_dlp :\n[",i+1,"]",g,"=",h)
			return i+1
	
	#puis rho standard
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


#retourne x tel que [x]g = h mod n
#p : ordre du groupe dont g est le générateur
#b : borne maximale pour les puissances de la fonction de marche aléatoire, p semble plus efficace que les petites valeurs
#k : nombre de parties
#pré-requis : p premier, h dans <g>
def rho_pollard_dlp_adv(g, h, p, n, b, k, verbose=False):
	#brute force pour commencer
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
	
#sous-routine pour rho pollard parallélisé
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
			
	

#retourne x tel que [x]g = h mod n
#p : ordre du groupe dont g est le générateur
#b : borne maximale pour les puissances de la fonction de marche aléatoire, p semble plus efficace que les petites valeurs
#k : nombre de parties
#jobs : nombre de sous-routines à lancer dans la parallélisation (à ajuster selon le nombre de coeurs logiques)
#pré-requis : p premier, h dans <g>
def rho_pollard_dlp_par(g, h, p, n, b, k, jobs=8, verbose=False):
	#brute force pour commencer
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
		

#Résout un DLP dans Z/nZ*
#g un générateur, on cherche x tq [x mod n-1]g = h mod n
def pohlig_hellman(g, h, n, log_file, verbose=False):
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

log_file = open("../results/rho_dlp_log.txt", "w+")

for l in range(30,70,5):
	print("\n\n ==============\n Taille",l,"bits\n ==============\n")
	log_file.write("\n\n====== Taille de l'entrée : " + str(l) + " bits =======\n")
	p = itools.rand_prime(2**(l-1),2**l,10)
	g,n = itools.get_group(p, True)
	print("N =",n,"générateur trouvé :",g)
	x = random.randint(2,n-1)
	h = itools.exp_mod(g,x,n)
	print("N =",n,", g =",g,", x =",x,", h =",h)
	if(l<60):
		for _ in range(20):
			x2 = pohlig_hellman(g,h,n,log_file,True)
	else:
		for _ in range(10):
			x2 = pohlig_hellman(g,h,n,log_file,True)

log_file.close()


