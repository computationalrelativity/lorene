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

	SUBROUTINE EXTM1S(NR,INR,IDR,DEN,SOM)

	IMPLICIT double PRECISION (A-H,O-Z)

C
C##	version du 23/11/93: derniere dimension des tableaux avec *
C			     ajout des SAVE
C
C		ROUTINE POUR LE CALCUL DES VALEURS ET DE DE LA DERIVE'E
C		PREMIERE EN r=0 ET r =2 D'UNE FONCTION SCALAIRE OU 
C		VECTORIELLE A 2-D AVEC ECHANTILLONNAGE DENSE A L'ORI-
C		GINE.
C
C
C		ARGUMENTS DE LA ROUTINE:
C
C		NR	=NOMBRE DES DEGRES DE LIBERTE-1
C		INR	=PARAMETRE: SI INR=0 LA FONCTION EST CALCULEE EN
C			 EN r=0, SI INR=1 EN r=2
C		IDR	=PARAMETRE: SI IDR=0 LA FONCTION EST CALCULEE
C			 SI IDR=1, LA DERIVEE 1ere EST CALCULEE
C		DEN	=TABLEAU IMPUT (A 2-D) CONTENANT LES COEFFICIENTS 
C			 DE TCHEBYTCHEV SELON LES 2 AXES.
C		SOM	=OUTPUT
C			 DES COEFFICIENTS SELON theta DE LA FONCTION
C
C		L'IMPUT DEN N'EST PAS DETRUIT
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:28:29  hyc
c *** empty log message ***
c
C Revision 1.1  1997/05/07 16:38:56  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	DIMENSION DEN(*)
	DIMENSION ID1(514)
	DATA NCONT/0/
C##
	SAVE NCONT, ID1
C
	N514=514
	IF(NR.GT.N514) THEN
C
	PRINT*,'DIMENSIONS INSUFF. DANS LA ROUTINE EXTM1S, NR=',NR
	CALL EXIT
	ENDIF
C
	NR1=NR+1
C
C		PREPARATION DES TABLEAUX
C
	IF(NCONT.NE.NR1) THEN
	NCONT=NR1
C
	DO 1 LR=1,NR1
	ID1(LR)=-(LR-1)**2
   1	CONTINUE
C
	ID1(NR1)=-NR**2/2
	ENDIF
C
C		CALCUL DE LA VALEUR DE LA FONCTION EN r=0
C
	IF(INR.EQ.0) THEN
C
	IF(IDR.EQ.0) THEN
C
	SOM=(DEN(1)+DEN(NR1))*.5
C
	DO 4 LR=2,NR
	SOM=SOM+DEN(LR)
   4	CONTINUE
C
	RETURN
	ENDIF
C
C		CALCUL DE LA DERIVEE EN r=0
C
	IF(IDR.EQ.1) THEN
C
	SOM=-DEN(2)
C
	DO 9 LR=3,NR1
	SOM=SOM+DEN(LR)*ID1(LR)
   9	CONTINUE
C
	RETURN
	ENDIF
C
		ENDIF
C
	IF(INR.EQ.1) THEN
C
C		CALCUL VALEUR DE LA FONCTION EN r=2
C
	IF(IDR.EQ.0) THEN
C
	SOM=(DEN(1)+DEN(NR1))*.5
C
	DO 14 LR=2,NR,2
	SOM=SOM-DEN(LR)
  14	CONTINUE
C
	DO 17 LR=3,NR,2
	SOM=SOM+DEN(LR)
  17	CONTINUE
C
	RETURN
	ENDIF
C
C		CALCUL DE LA FONCTION EN r=2
c
	IF(IDR.EQ.1) THEN
C
	SOM=-DEN(2)
C
	DO 22 LR=3,NR1,2
	SOM=SOM-DEN(LR)*ID1(LR)
  22	CONTINUE
C
	DO 25 LR=4,NR,2
	SOM=SOM+DEN(LR)*ID1(LR)
  25	CONTINUE
C
	RETURN
	ENDIF
C
	ENDIF
C
	RETURN
	END
