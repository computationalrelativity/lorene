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

	subroutine tfm1s(n,n64,y,cc)

	implicit double precision(a-h,o-z)

c
c## routine modifiee le 04.04.1994: ajout des save
c				    suppression de la variable jjj non utilisee
c
c		routine pour le calcul de les tf de n64 fonctions
c		a la fois.
c
c		arguments de la routine:
c
c		n	= nombre des degres de libertee. n doit
c			  etre= a 2**p*3**q*5**r avec p,q,r nombrs
c			  entiers (p.ne.0).
c
c		n64	= nombres des fonctions dont on veut calculer
c			  les tf.
c		y	= tableau contenant en imput les echantillons des fon-
c			  ctions. y doit avoir (n+4)*(n64+1)  dimensions.
c			  les n echantillons des n64 fonctions sont
c			  stockes en parallel ou en serie. en paralel
c			  on a:
c			  y(1),y(2),...y(n64) contient les valeures
c			  des n64 fonctions dans le point x=0,
c			  y(1+n64),y(2+n64),...y(n64+n64)
c		 	les valeures des fonctions pour x=2*pi/n
c		  	y(1+2*n64),y(2+2*n64),...y(n64+2*n64)
c		 	les valeures pour x=(2*pi/n)*2 et ainsi de 
c			 jusqu' a x=(2*pi/n)*(n-1).
c
c		cc	= output: les coeeficients de cc impaires contienent
c			 les cosinus et les coefficients paires les sinus.
c			 (cc(2)=0 par definition. l'argument de cc est donc
c			 compris entre 1 et n+1. 
c
c		routine testee le 14/9/1985 
c
C
C $Id$
C $Log$
C Revision 1.2  2012/03/30 12:12:44  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.1  1997/10/23  08:19:46  eric
c Initial revision
c
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/


	dimension cc(*),y(*),ifax(64),trigs(1600)
	data ndim/0/
	data nfon/0/
c##
	save ndim, nfon, ifax, trigs, n65, n63, n6565, nm65, n66

	n1040=1040
	if(n.gt.n1040) then
	print 800,n
  800	format(10x,'dimension insuffisantes dans la routine'
     ,	,' tfm, n=',i5,' > a la dimension de trigs*2/3=',i5)
		endif
c
	if(n.eq.ndim) go to 1
	call fax(ifax,n,3)
	call fftrig(trigs,n,3)
  1	continue
c
	if(nfon.eq.n64.and.ndim.eq.n) go to 4
	nfon=n64
	ndim=n
	n65=n64
	if((n64/8)*8.eq.n64)n65=n64+1
	n63=n64-1
	n6565=n65+n65
	nm65=n*n65
	n66=n65+1
  4	continue
c
c
c
	do 10 l=1,nm65
	cc(l)=y(l)*2
10	continue
c
c		calcul de la tf en parallel
c
	call fft991(cc,y,trigs,ifax,n65,1,n,n64,-1)
c
  101	format(1x,'tfm')
  100	format(1x,'tfm',10d12.4)
c
	return
	end

