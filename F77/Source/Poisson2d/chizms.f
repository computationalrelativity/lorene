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


	subroutine chizms(n,n64,f,cc,cs)
c
c## routine modifiee le 01.04.1994: ajout des save
c				     et replacement de nmax par n513
c

	implicit double precision(a-h,o-z)

c
c		routine pour la transformation inverse simultanee
c		de n64 fonctions. la routine est completement
c		craytinisee.
c
c	arguments de la routine:
c	
c	n	=nombre des degrees de libertee-1. n doit etre
c		 pair et = a 2**p*3**q*5**m avec p,q,m nombres
c		entiers.
c	n64=	nombre des fonction qui doivent etre transformees.
c	f	= coefficients de tchebytchev des fonctions a transfor-
c		 mer. les coefficients sont stockes 'en parallel'
c		i.e. dans f(1),f(2)....f(n64) il-y-a le premier
c		coefficient des n64 fonctions, dans f(n64+1),f(n64+2),
c		...f(n64+n64) le 2me coefficiet des n64 fonctions,
c		et ainsi de suite. si n64 est un multiple de 8, pour
c		de craytinisation les
c		coefficients sont stockes dans la facon suivante:
c		f(1),f(2),...f(n64) pour le premier coefficient
c		f(n64+1+1),f(n64+1+2),...f(n64+1+n64) pour le 2me
c		coefficient, et ainsi de suite.
c	cc	= tableau de travail
c	cs	= tableau output: dans cs il-y-a le n1 valeures des n64
c		  fonctions ('en parallel')
c		  les dimensions minimum des tableaux
c		  sont (n+3)*(n64+1).
c		 le tableau f est detruit.
c
c	cette routine est specialisee: elle doit etre employee avec fouci2s
c		fuci3s...
c
c		tous les test de routine ont ete executes l 8/10/85.
c
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:32:16  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:35:13  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	dimension cc(*),cs(*),f(*)
	dimension sen(513),f11(513),som1(513),fn21(513),somm(513)
	data ndim/0/
	data nfon/0/
c##
	save ndim, nfon, n1, n2, n21, an2, x0, pi, sen, n65
	save n63, n65n65, nm65
	save n3, n651, nm650, nm623, nm20, nm201, nm2064

	n513=513
	if(n.lt.n513) go to 12
	print 200,n,n513
  200	format(10x,'DIMENSIONS INSUFFISANTES DANS CHIZMS, N=',i4,
     ,	' DIMENSIONS MAX.=',i4)
	call exit
c
  12	continue
c
c		preparation des quantites necessaires pour le calcul.
c		ces quantites sont calculees la premiere fois qu'on
c		appelle la routine et toutes les fois q'on chenge
c		la valeur de n.
c
	if(n.eq.ndim) go to 10
	n1=n+1
	n2=n/2
	n21=n2+1
	n3=n2-1
	x0=0
	pi=2.*acos(x0)/n
	do 11 l=2,n2
	sen(n21-l)=.25/sin(pi*(l-1))+.5
  11	continue
  10	continue
c
c		preparation des quantites necessaires au clacul.
c		ces quantites sont calculees la premiere fois 
c		qu'on appelle la routine et toutes les fois 
c		qu'on change n64.
c
	if(nfon.eq.n64.and.ndim.eq.n) go to 20
	nfon=n64
	ndim=n
	n65=n64
	if((n64/8)*8.eq.n64) n65=n64+1
	nm65=n65*n
	n651=n65+1
	nm650=nm65-n65
	nm655=nm65+n65
	nm623=(n2-1)*n65
	n63=n64-1
	n65n65=n65+n65
	nm20=n2*n65
	nm201=nm20+1
	nm2064=nm20+n64
  20	continue
c
c
c ****   calcul de la fonction en teta=0 et teta=pi
c
	do 7 m=1,n64
	som1(m)=0
	somm(m)=0
  7	continue
c
	do 8 l=n65n65,nm650,n65n65
	do 30 m=1,n64
	somm(m)=somm(m)+f(m+l)
30	continue
   8	continue
	do 5 l=n65,nm65,n65n65
	do 31 m=1,n64
	som1(m)=som1(m)+f(m+l)
31	continue
  5	continue
c	
	do m=1,n64
	cc(m)=(f(m)+f(m+nm65))*.5
	enddo
c
	do 32 m=1,n64
	f11(m)=somm(m)+som1(m)+cc(m)
	fn21(m)=somm(m)-som1(m)+cc(m)
32	continue
c
	do 33 l=1,nm655
	cs(l)=0
33	continue
c
c		la boucle suivante est equivalente a:
c
	jl1=n651
	jl2=jl1+n63
	jl3=jl1+n65n65
	jl4=jl3+n63
	do l=2,n,2
	do m=jl1,jl2
	cs(m)=(cs(m)-f(m))
	enddo
c
	do m=jl3,jl4
	cs(m)=cs(m)+f(m-n65n65)
	enddo
	jl1=jl1+n65n65
	jl2=jl1+n63
	jl3=jl3+n65n65
	jl4=jl3+n63
	enddo
c		
	jl1=1
	jl2=jl1+n63
	do l=1,n1,2
	do m=jl1,jl2
	cs(m)=f(m)
	enddo
	jl1=jl1+n65n65
	jl2=jl1+n63
	enddo
c
	call tfizms(n,n64,cs,cc)
c
c		la boucle suivante est equivalente a:
c
c			lsen=1
c		do 3 l=1,n2-1
c			n21l=(n21+l-1)*n65
c			n20l=(n21-l-1)*n65
c			senn=sen(lsen)
c			do m=1,n64
c			f12=(cc(m+n20l)-cc(m+n21l))*senn
c			cs(m+n21l)=cc(m+n21l)+f12
c			cs(m+n20l)=cc(m+n20l)-f12
c			enddo
c			lsen=lsen+1
c 3			continue
c
c
	lsen=1
	do 41 l=n65,nm623,n65
	n21l=nm201+l
	n20l=nm20-l
	l2=l+l
	senn=sen(lsen)
	do 40 m=n21l,n63+n21l
	f12=(cc(m-l2)-cc(m))*senn
	cs(m)=cc(m)+f12
	cs(m-l2)=cc(m-l2)-f12
40	continue
	lsen=lsen+1
41	continue
c
	do 42 m=1,n64
	cs(m)=f11(m)
	cs(m+nm65)=fn21(m)
42	continue
	do 43 m=nm201,nm2064
	cs(m)=cc(m)
43	continue
c
100	format(1x,'CHIZMS',10d12.4)
101	format(1x,'CHIZMS')
110	format(1x,'CHIZMS',5d24.16)
120	format(1x,'CHIZMS',20i5)
300	format(10x,'DANS CHIZMS LMAX=',i5)
	return
	end
