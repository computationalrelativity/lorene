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
       SUBROUTINE CONR3S(NDEG,NDIMX,NDIMY,CC,YY)

		implicit double precision (a-h,o-z)

C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.1  1997/10/23  08:25:19  eric
c Initial revision
c
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/


       DIMENSION CC(NDIMX,NDIMY,*),YY(NDIMX,NDIMY,*),NDEG(3)
C
       NX=NDEG(1)
       NY=NDEG(2)
       NZ=NDEG(3)
       
       DO 3 LZ=1,NZ
       DO 2 LY=1,NY
       DO 1 LX=1,NX
       YY(LX,LY,LZ)=CC(LX,LY,LZ)
1      CONTINUE
2      CONTINUE
3      CONTINUE
       RETURN
       END
