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


	subroutine cd2yms(n,n64,cc,cs)
c
c## routine modifiee le 08.10.1993: ajout des save
c

	implicit double precision(a-h,o-z)

c
c***** calcul de la derivee 2me de n64 fonctions simultanement.
c****	subroutine completement craytinisee et specialisee pour etre
c	appelle par fuci2s ou fuci3s.
c
c		arguments de la subroutine:
c
c	n	=nombre des degres de libertee-1
c	n64	=nombre des fonctions dont on veut calculer la derivee.
c	cc	=tableau imput contenant les coefficients de tchebytchev
c		 des fonctions a deriver.
c	cs	=output contenant le derivees des fonctions.
c
c		toutes les donnes sont stockees dans le mode
c		'parallel' (voir les subroutines tfm ou,tfinm,chem64,
c		 chinm64 pour une explication detaillee')
c
c		subroutine specializee pour etreappelle par fuci2s
c		et fuci3s
c
c		subroutine ayant aubit tous les test de routine le 
c		16/12/86.
c
C
C $Id$
C $Log$
C Revision 1.2  2012/03/30 12:12:42  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.2  1997/05/23  11:34:25  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:21:05  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

	dimension cc(*),cs(*),somm(520),som1(520)
	data nfon/0/
	data ndim/0/
c##
	save ndim, nfon, n1, n11, nn, n65, nm65, nm651, n65n65
	save nm065, nm965, n10

	n520=520
	if(n64.gt.n520) then
	print 200,n520,n64
  200	format(10x,'dimension insuffisantes dans la sub.ced2m64'
     ,	,' dimensions max. admise=',i4,' n64=',i5)
	call exit
	endif
c
	if(n.eq.ndim.and.nfon.eq.n64) go to 10
	ndim=n
	n1=n+1
	n11=n+2
	n10=n-1
	nn=n+n
c
	nfon=n64
	n65=n64
	if((n64/8)*8.eq.n64)n65=n64+1
	nm65=n*n65
	nm651=nm65+n65
	n65n65=n65+n65
	nm065=(n-1)*n65
	nm965=nm065-n65
  10	continue
c
c
	do 4 l=1,nm651
	cs(l)=0
4	continue
c
	do 11 m=1,n64
	somm(m)=0
	som1(m)=0
  11	continue
c
	do 12 m=nm65+1,nm651
	cc(m)=cc(m)*.5
12	continue
c
	al1=n
	lm1=nm965
	lm=nm65
	do 1 lp=1,n10,2
	al2=al1**2
	al32=(al1-2)**2
	do 13 m=1,n64
	cpm=cc(m+lm)*al1
	somm(m)=somm(m)+cpm
	som1(m)=som1(m)+cpm*al2
13	continue
c
	do 14 m=1,n64
	cs(m+lm1)=som1(m)-somm(m)*al32
14	continue
c
	al1=al1-2.d+00
	lm1=lm1-n65n65
	lm=lm-n65n65
1	continue
c
	do 15 m=1,n64
	som1(m)=0
	somm(m)=0
15	continue
c
	al1=n-1
	lm=nm065
	lm1=nm965
	do 2 lp=2,n10,2
	al2=al1**2
	al32=(al1-2)**2
	do 16 m=1,n64
	cpm=cc(m+lm)*al1
	somm(m)=somm(m)+cpm
	som1(m)=som1(m)+cpm*al2
16	continue
c
	do 17 m=1,n64
	cs(m+lm1-n65)=som1(m)-somm(m)*al32
17	continue
c
	al1=al1-2
	lm=lm-n65n65
	lm1=lm1-n65n65
c
2	continue
c
  100	format(1x,'ced2m64',10d12.4)
  101	format(1x,'ced2m64')
  110	format(1x,'ced2m64',20i5)
c
	return
	end
