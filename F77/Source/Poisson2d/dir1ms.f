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

	subroutine dir1ms(n,n64,ind,cs,cc)

	implicit double precision(a-h,o-z)

c
c		routine pour la multiplication et la division
c		de n64 fonctions par r. (div (par) r (avec resolution
c		f(in) a l'origine
c
c		arguments de la routine:
c
c		n	=nombrede degres de liberte-1
c			 cas contraire.
c		ind	=parametre: si mult=1 la multiplication par x est
c			 effectuee si mult=-1 la division est effectuee.
c		n64	=n64 est le nombre de fonctions qui doivent etre
c			 divisees (ou multipliees) par x
c		cs	=tableau contenant le n+1 coefficients de tcheby-
c			 tchev de la fonction qui doit etre divisee ou
c			 multipliee par x.
c		cc	=tableau contenant les coeff. de la fonction divisee.
c			 les dimensions des tableaux cc et cs doibent
c			 etre > (n64+1)*n1
c
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:29:46  hyc
c *** empty log message ***
c
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	dimension cs(*),cc(*),y(513)
	n1=n+1
	if(n64.gt.513) then
	print*,'dimensions insuff. dans la routine dir1ms, n64=',n64
	endif
c
	n0=n-1
	n65=n64
	if((n64/8)*8.eq.n64) n65=n64+1
	n63=n64-1
	n064=(n-2)*n65+1
	nn64=n064+n65
	n164=nn64+n65
c
	if(ind.eq.1) then
c
c		on effectue la multiplication
c
	do 13 j=1,n64
	cc(j)=cs(j)-cs(j+n65)
  13	continue
c
	jm1=n65+1
	jm2=jm1+n63
	do 1 l=2,n0
	do 2 m=jm1,jm2
	cc(m)=cs(m)-(cs(m+n65)+cs(m-n65))*.5
  2	continue
	jm1=jm1+n65
	jm2=jm1+n63
  1	continue
	do 5 m=nn64,nn64+n63
	cc(m)=cs(m)-cs(m-n65)*.5
  5	continue
c
	do 6  m=nn64,nn64+n63
	cc(m+n65)=-cs(m)
  6	continue
	return
	endif
c
c		on effectue la division
c
		if(ind.eq.-1) then
c
	jm1=n65+1
	jm2=jm1+n63
	do 3 l=1,n0
	do 7 m=jm1,jm2
	cc(m)=cs(m-n65)+2*cs(m)
  7	continue
	do 8 m=jm1,jm2
	cs(m)=cc(m)
  8	continue
	jm1=jm1+n65
	jm2=jm1+n63
  3	continue
c	
	do 9 m=n164,n164+n63
	cc(m)=0
  9	continue
c
  	do 10 m=nn64,nn64+n63
	cc(m)=cs(m)
  10	continue
c
	jm1=n064-1
	jm2=nn64-1
	do 4 ll=2,n
	do 11 m=1,n64
	y(m)=cs(m+jm1)+cc(m+jm2)
  11	continue
c
	do 12 m=1,n64
	cc(m+jm1)=y(m)
  12	continue
c	
	jm2=jm1
	jm1=jm1-n65
  4	continue
	endif
  100	format(1x,10e11.3)
  101	format(1x,' ')
	return
	end
