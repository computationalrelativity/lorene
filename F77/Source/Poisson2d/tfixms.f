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


	subroutine tfixms(n,n64,cc,y)
c
c## routine modifiee le 08.10.1993: ajout des save
c

	implicit double precision(a-h,o-z)

c
c		routine pour la tf inverse multiple.
c
c		arguments de la suroutine:
c
c		n	=nombre des dgres de libertee. (n nombre pair=
c			 a 2**p*3**q*5**r p,q,r nombres entiers.
c
c		cc	=coefficients cosinus et sinus de n64 fonctions.
c		y	= output des fonctions dans l'espace des x.
c			  les valeures des ces fonctions sont stockees
c			  comme dans la sub tf.
c		cette routine est specialisee et ne doit etre employee
c		que avec les routines chixms,fuci2s,fuci3s...
c
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:36:59  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:22:08  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	dimension y(*),cc(*),ifax(64),trigs(1600)
	data ndim/0/
	data nfon/0/
	data jjj/0/
c##
	save ndim, nfon, jjj, ifax, trigs, n65, nm65

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
c
  1	continue
	if(nfon.eq.n64.and.ndim.eq.n) go to 4
c
c		si n64 est un multiple de 8 le temps calcul du cray
c		etre multiplie par un  facteur 4. pour eviter
c		cela on socke les donnees tous les n64+1 intervalles
c		au lieu de n64.(dans le cas de l input en parallel)
c
	nfon=n64
	ndim=n
	n65=n64
	if((n64/8)*8.eq.n64) n65=n64+1
c
	nm65=n*n65
c
  4	continue
c
c		tf inverse en parallel
c
	call fft991(cc,y,trigs,ifax,n65,1,n,n64,1)
c
	do 12 l=1,nm65
	y(l)=cc(l)*.5
12	continue
  101	format(1x,'tfm')
  100	format(1x,'tfm',10d12.4)
	return
c
	end
