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
	SUBROUTINE CONRS(N,F,FC)

	IMPLICIT double PRECISION (A-H,O-Z)

C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/10/30  13:42:02  novak
c *** empty log message ***
c
C Revision 1.1  1997/10/21 14:03:23  eric
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	dimension F(*),FC(*)
	N1=N+1
	DO 1 L=1,N 1
	FC(L)=F(L)
1	CONTINUE
	RETURN
	END
