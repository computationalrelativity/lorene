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

	SUBROUTINE PRIMS(N,CC,CS,INDEX,Y)

	IMPLICIT double PRECISION (A-H,O-Z)

C		
C		ROUTINE POUR LE CALCUL DE LA PRIMITIVE D'UNE FONCTION
C	EN ENTREE:
C		N	=NOMBRE DES COEFFICIENTS DE TCHEBYTCHEV-1.
C		INDEX	=PARAMETRE: INDEX=1 ON A EN CS LES COEFFICIENTS
C			 DE TCHEBYTCHEV DE LA PRIMITIVE.
C			 INDEX=2 ON A EN Y LA PRIMITIVE DE LA FONCTION
C		CC	=COEFFICIENTS DE TCHEBYTCHEV DE LA FONCTION DONT
C			 ON VEUT CALCULER LA PRIMITIVE.
C
C	EN SORTIE:
C		CS	=TABLEAU CONTENANT LES COEFFICIENT DE TCHEBYTCHEV
C			 DE LA PRIMITIVE (SEULEMENT SI INDEX=1)
C		Y	=PRIMITIVE DE LA FONCTION.(SEULEMENT SI INDEX.NE.1)
C
C	POUR INDEX.NE.1 L'IMPUT EST DETRUIT
C
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:38:21  hyc
c *** empty log message ***
c
C Revision 1.1  1997/05/08 07:35:32  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	DIMENSION CC(*),CS(*),Y(*)
	N1=N+1
	CS(2)=(CC(3)-CC(1))*.5D+00
	DO 3 L=3,N
	L1=L-1
	CS(L)=.5D+00*(CC(L+1)-CC(L1))/(L1)
   3	CONTINUE
	CS(N1)=-CC(N)/N
  200	FORMAT(10X,'INDEX=',I4)
	IF(INDEX.NE.1) GO TO 5
	SOMM=0
	DO 2 L=2,N
	SOMM=SOMM+CS(L)
  2	CONTINUE
	CS(1)=-SOMM-SOMM-CS(N1)
	IF(INDEX.EQ.1) RETURN
  5	CONTINUE
	CALL CHINS(N,CS,CC,Y)
	Y1=Y(1)
	DO 4 L=1,N1
	Y(L)=Y(L)-Y1
   4	CONTINUE
  100   FORMAT(1X,'PRIM',10D12.4)
  101	FORMAT(1X,'PRIM')
	RETURN
	END
