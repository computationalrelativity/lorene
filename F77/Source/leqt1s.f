C
C   Copyright (c) 1997 Silvano Bonazzola
C   Copyright (c) 1998 Jerome Novak
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

	subroutine leqt1s(a,n,nlc,nuc,ia,b,m,ib,ijob,xl,ier)

C-------------------------------------------------------------------------
C Version tampon de LEQT1B qui utilise lapack
C
C Relation entre la matrice M(i,j) et son ecriture a bande MB :
C 1/ pour IMSL
C      M(i,j) = MB(l,b) <-> i=l et j=l+b-nlc-1
C
C 2/ pour LAPACK
C	M(i,j) = MB(nlc+nuc+1+i-j,j),	1 =< j =< n
C					max(1,j-nuc) =< i =< min(n,j+nlc)
C-------------------------------------------------------------------------

			implicit none

C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.3  1998/09/30  12:56:23  novak
c parametre nloc mis a 2049 (suffisant??)
c
c Revision 1.2  1997/05/23  11:39:32  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/19 19:24:51  eric
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

C Variables d'appel
C------------------
	integer	n,nlc,nuc,ia,m,ib,ijob,ier
	real*8	a(ia,*),xl(ia,*),b(ia,*)

C Variables locales
C------------------
		integer nloc,nbloc
		parameter (nloc=2049,nbloc=20)

	real*8	ab(nbloc,nloc)
	integer	ipiv(nloc)
	integer	nb,i,j

C Protection
C-----------
	nb = 2*nlc + nuc + 1
	if (n.gt.nloc.or.nb.gt.nbloc) stop 'Dim LeqLap'

C Traduction de la matrice
C-------------------------

		do j=1,n
	    do i=max(1,j-nuc),min(n,j+nlc)
	ab(nlc+nuc+1+i-j,j) = a(i,j-i+nlc+1)
	enddo
	enddo

	call dgbsv(n,nlc,nuc,m,ab,nbloc,ipiv,b,ib,ier)
	if (ier.ne.0) write (*,*) 'Peut-etre ib : ier',ier,ib

	return
	end
