C
C   Copyright (c) 1997 Silvano Bonazzola
C
C    This file is part of LORENE.
C
C    LORENE is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    LORENE is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with LORENE; if not, write to the Free Software
C    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C


	subroutine citxms(n,n64,index,f,cc,cs)
c
c## routine modifiee le 08.10.1993: ajout des save
c				     et suppression des variables declarees 
c				     mais non utilisees
c

	implicit double precision(a-h,o-z)

c
c **** inversion de la tf de tchebytchev de n64 fonctions echantillonees
c **** rarement a l'origine. (intervalle d'echantillonage [0-1]
c **** dimension minimale des tableaux=n+1
c **** routine distruptive
c
c		arguments de la routine:
c
c	n	= nombre des degres de liberte'
c	n64	=nombre des fonctions a inverser.
c	index	= parametre: il doit etre =0 si la fonction d'entree
c		  est paire, =1 si impaire.
c	f	= imput: tableau contenant les coefficients de tche-
c		bitchev des fonctions a inverser.
c		 f doit contenir les coefficients de n64 fonctions 'en paral-
c		 lel' (cfr. la routine tfnms). si le nombre des fonctions
c		 a inverser est un multiple de 8, n64 doit etre, pour des rai-
c		 sons de craytinisation = le nombre des fonctions=1.
c
c	cc	= tableau de travail.
c	cs	= output.
c
c			les dimensions des tableaux doivent etre gt.eq.
c			(n64+1)*(n+2)
c			
c		routine testee le 4/2/1987.
c
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:31:27  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:33:57  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	dimension cc(*),f1(257),fn(257),cs(*),f(*),sen(1025),cose(1025)

	data ncontr/0/,nfonc/0/
c##
	save ncontr, nfonc, m1025, n1, n11, n2, n21, n3, x0, pi
	save pi2, sen, cose, n257, n65, n63, n66, nm65, nm651, nm66
	save n2651, n265, n6565, n633, n365
c
	if(ncontr.eq.n) go to 999
c
	m1025=257
	if(n.lt.m1025) go to 9
	print 200,n
  200	format(10x,'DIMENSIONS INSUFFISANTES DANS CITXMS: N=',i5)
	print*,'EXECUTION INTERROMPUE'
	call exit
  9	continue
c
	n1=n+1
	n11=n+2
	n2=n/2
	n21=n2+1
	n3=n2-1
	x0=0
	pi=2.*acos(x0)/n
	pi2=pi/2
	sen(1)=.25/sin((n2-1)*pi)+.5
	do 1 l=2,n3
	sen(l)=.25/sin((n2-l)*pi)+.5
	cose(l)=1/sin((l-1)*pi2)
  1	continue
	do 2 l=n21,n1
	sen(l)=.25/sin((n2-l)*pi)+.5
	cose(l)=1/sin((l-1)*pi2)
  2	continue
	cose(n2)=1/sin(n3*pi2)
c
  999	continue
c
		if(n.ne.ncontr.or.n64.ne.nfonc) then
c
	n257=257
	if(n64.gt.257) then
	print*,'dimension insuffisantes dans la sub. cirxms, n64=',n64
	print*,'execution interrompue'
	call exit
	endif
c
	ncontr=n
	nfonc=n64
c
	n65=n64
	if((n64/8)*8.eq.n64)n65=n64+1
	n63=n64-1
	n66=n65+1
	nm65=n*n65
c	nm651=nm65+1
	nm66=nm65+1
	n2651=n21*n65
	n265=n2651-n65
	n6565=n65+n65
	n633=n65-1
	n365=n2651-n6565
	endif
c
c
	if(index.eq.0) go to 11
c
	do 3 m=1,n64
	cs(m)=-f(m)
  3	continue

	do 4 m=nm65-n633,nm65
	cc(m)=-f(m)
  4	continue
c
	do 5 m=nm66,nm65+n65
	f(m)=cc(m-n65)
  5	continue
c
	jm1=n65+1
	jm2=jm1+n63
	do 6 l=2,n
	do 7 m=jm1,jm2
	cs(m)=-(f(m)+f(m-n65))*.5
  7	continue
c
	jm1=jm1+n65
	jm2=jm1+n63
  6	continue
c
	do 8 l=1,nm65
	f(l)=cs(l)
  8	continue
  11	continue
c
c ****   calcul de la fonction en teta=0 et teta=pi
c
	do 12 m=1,n64
	cc(m)=0
  12	continue
c
	jm1=n6565
	do 13 l=3,n,2
	do 14 m=1,n64
	cc(m)=cc(m)+f(m+jm1)
  14	continue
	jm1=jm1+n6565
  13	continue
c
	do 15 m=1,n64
	f1(m)=0
  15	continue
c
	jm1=n65
	do 16 l=2,n ,2
	do 17 m=1,n64
	f1(m)=f1(m)+f(m+jm1)
  17	continue
	jm1=jm1+n6565
  16	continue
	do 18 m=1,n64
	cs(m)=(f(m)+f(m+nm65))*.5
  18	continue
c
	do 19 m=1,n64
	cs(m)=cs(m)+cc(m)
  19	continue
c
	do 20 m=1,n64
	fn(m)=cs(m)+f1(m)
  20	continue
c
	do 21 m=1,n64
	f1(m)=cs(m)-f1(m)
  21	continue
c 
c ****  calcul de coef. de fourier de la foncyion ponderee a partir de c
c ****  de tchebychev.
c
	jm1=n66
	jm2=n65+n64
	do 22 m=jm1,jm2
	cs(m)=0
  22	continue
c
	do 23 l=2,n,2
	do 24 m=jm1,jm2
	cs(m)=cs(m)+f(m)
  24	continue
c
	do 25 m=jm1,jm2
	cs(m+n6565)=-f(m)
  25	continue
	jm1=jm1+n6565
	jm2=jm1+n63
  23	continue
c
	jm1=n66
	jm2=jm1+n63
	do 26 m=jm1,jm2
	f(m)=0
  26	continue
c
	jm1=jm1+n6565
	jm2=jm1+n63
	do 27 l=4,n,2
	do 38 m=jm1,jm2
	f(m)=-cs(m)
  38	continue
	jm1=jm1+n6565
	jm2=jm1+n63
  27	continue
c
	call tfiyms(n,n64,f,cc)
c
c	 ****   restitution de la fonction

c
	jm1=n2651
	jm2=n365
	jm3=1
	jm4=n64
	do 28 l=1,n3
	senn=sen(l)
	do 29 m=jm3,jm4
	f(m)=(cc(m+jm2)-cc(m+n2651))*senn
  29	continue
	jm1=jm1+n65
	jm2=jm2-n6565
	jm3=jm3+n65
	jm4=jm3+n63
  28	continue
c
	jm2=n3*n65
	jm3=1
	jm4=n64
	do 30  l=1,n3
	do 31 m=jm3,jm4
	cs(m+jm2)=cc(m+n2651)+f(m)
  31	continue
	do m=jm3,jm4
	cs(m+n2651)=cc(m+jm2)-f(m)
	enddo
c
	jm2=jm2-n6565
	jm3=jm3+n65
	jm4=jm3+n63
30	continue
c
	do 32 l=1,n64
	cs(l)=f1(l)
32	continue
c
	jm1=n265+1
	jm2=jm1+n63
	do 33 m=jm1,jm2
	cs(m)=cc(m)
  33	continue
c
	jm1=nm65
	jm2=jm1+n63
	do 34 m=1,n64
	cs(m+jm1)=fn(m)
  34	continue
c
	if(index.eq.0) return
c
	jm1=n65+1
	jm2=jm1+n63
	do 35 l=2,n1
	cosen=cose(l)
	do 36 m=jm1,jm2
	cs(m)=cs(m)*cosen
  36	continue
c
	jm1=jm1+n65
	jm2=jm1+n63
  35	continue
	do 37 m=1,n64
	cs(m)=0
  37	continue
c
100	format(10x,10d12.4)
101	format(1x,'in')
	return
	end
