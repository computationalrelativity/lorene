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

	subroutine chezms(n,n64,f,cc,cs)

	implicit double precision(a-h,o-z)

c
C 	version du 04.12.1996
c
c## routine modifiee le 08.10.1993: ajout des save
c				     et suppression des variables declarees 
c				     mais non utilisees
c
c  ***	routine pour le calcul simultane des coeff. de tche-
c  ***	bytchef de n64 fonctions. la routine est completment
c  ***	craytinisee.
c  ***	f doit etre echantillonee dans n+1 points.. dimension minimal
c  *** 	des tableaux: =(n64+1)*(n+3), ou n+1 est le nombre des points
c  ***	d'echantillonage.
c
c		cette routine est legerement plus rapide que la chems.
c		l'output est dans cc au lieu de cs
c
c
c		arguments de la routine:
c
c	n	= nombre des degres de liberte-1.
c	n64	= nombre des fonctions d'ont on veut calculer la trans-
c	 	formation de tchebytchev. 
c	f	=tableau imput. les valeurs echantillonnes des n64
c		 fonctions sont stockees en 'parallel': dans f(1),f(2)..
c		 f(n64) on stocke le valeures de n64 fonctions en x=x0
c		 dans f(1+n64),f(2+n64),...f(n64+n64), les valeures
c		 des n64 fonctions dans le point x1 et ainsi de suite.
c		   si n64 est un multiple de 8, pour des raisons de
c		 craytinisation les donnes sont stockees dan la fa-
c		 con suivante: f(1),f(2),...f(n64), f(n64+1+1),
c		 f(n64+1+2),...f(n64+1+n64).
c	cc	=output. les coeff. de tchebytchev sont sto-
c			 en 'parallel'.
c	cs	=tableau de travail
c
c
c	tous les tests de routine ont ete effectues les 20/11/1986
C
C $Id$
C $Log$
C Revision 1.2  2012/03/30 12:12:42  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.1  1997/10/23  08:16:43  eric
c Initial revision
c
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

	dimension cc(*),cs(*),f(*)
	dimension sen(513),fn11(513),y(513)
	data ndim/0/
	data nfon/0/
c##
	save ndim, nfon, n1, n2, n21, an2, x0, pi, sen, n65
	save n63, n65n65, nm65, nm651, nm652, n6519, jinit
c
	n513=513
	if(n.lt.n513) go to 12
  200	format(10x,'dimensions insuffisantes dans chem64,n=',i4)
	call exit
c
c		preparation des quantites necessaires pour le calcul.
c		ces quantites sont calculees la
c		premiere fois qu'on appelle la routine.
c		elles sont recalculees toutes les fois qu'on change
c		n.
c
  12	continue
c	if(n.eq.ndim) go to 10
	n1=n+1
	n2=n/2
	n21=n2+1
	an2=1.d+00/n2
	x0=0
	pi=2*acos(x0)/n
	do 14 l=1,n1
	sen(l)=sin((l-1)*pi)
  14	continue
  10	continue
c	if(nfon.eq.n64.and.n.eq.ndim) go to 18
c
c		preparation des quantites necessaires pour le calcul.
c		ces quantites sont calculees la premiere fois qu'on
c		appelle la routine. elles sont recalculees toutes
c		les fois qu'on change n64.
c
	n65=n64
	ndim=n
	if((n64/8)*8.eq.n64) n65=n64+1
	n63=n64-1
	n65n65=n65+n65
	nm65=n65*n
	nm651=nm65+n65
	nm652=n65*n21
	n6519=n65*(n21-2)+1
	jinit=3*n65+1
	nfon=n64
  18	continue
c
c
c		pondesration des la partie antisimmetrique des fonc-
c		tions.
c
c
	do 3 l=1,n65
	fn11(l)=   (f(l)-f(l+nm65))*.5
  3	continue
c
	n21l=nm652
	n20l=n6519
	do 1 l=1,n2
	n2lx=n21l-n20l+1
	sen20l=sen(n21-l)+.5
	lx64=n20l+n63
	do 6 lx=n20l,lx64
	cs(lx)=(f(lx)-f(lx+n2lx))*sen20l
  6	continue
	do 20 lx=n20l,lx64
	f(lx)=f(lx)-cs(lx)
	f(lx+n2lx)=f(lx+n2lx)+cs(lx)
  20	continue
	n21l=n21l+n65
	n20l=n20l-n65
1	continue
c
c		calcul de la t.f. fonction precedente.
c
	call tfm3s(n,n64,f,cc)
c
	do 4 m=1,n65
	y(m)=0
	cs(m+n65)=0
  4	continue
c
	jl=n65n65
	do 9 l=jinit,nm65,n65n65
	jl2=l-1
	n65l=l+n63
	do 5 m=l,n65l
	f(m)=cs(m-n65n65)-cc(m)
   5	continue
	do 7 m=l,n65l
	cs(m)=f(m)
   7	continue
	do 8 m=1,n65
	y(m)=y(m)+cs(m+jl2)
  8	continue
  9	continue
c
	do 13 m=1,n65
	f(m)=(fn11(m)-y(m))*an2
  13	continue
c
	jm0=n65
	do 11 l=1,n2
	jm=jm0+1
	jm65=jm0+n65
	do 15 m=jm,jm65
	cc(m)=cs(m)+f(m-jm0)
  15	continue
	jm0=jm0+n65n65
  11	continue
  110	format(1x,20i5)
  100	format(1x,10d12.4)
101	format(1x,'chem64')
	return
	end
c	

