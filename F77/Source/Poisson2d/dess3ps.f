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
	SUBROUTINE DESS3PS(NDL,NDR,NDT,NDF,JCOEFF,IPAR,INDD,SCAL)
C
C		ROUTINE POUR DESIALISER UN SCALAIRE.
C		LES DERNIERS JCOEF COEFFICIENTS SONT MIS A ZERO,
C		 LE RACCORDEMENT A TRAVERS LES DIFFERENTES
C		COQUILLES EST PRESERVE.
C
C	ARGUMENTS DE LA ROUTINE:
C
C			NDL	=TABLEAU CONTENENT LES DEGRES DE LIBERTE
C				 DE DANS LES DIFFERENTES DIMENSIONS:
C		                 NDL(1)	=NZON: NOMBRE DE COQUILLES
C				 NDL(2),NDL(3,),...NDL(NZON+1)= NOMBRE
C				 DES DEGRES DE LIBERTE EN r DANS LE DIFFE-
C				 RENTES COQUILLES.
C				 NDL(NZON+2),NDL(NZON+3)=
C				 NOMBRE DES DEGRES EN theta ET phi RESPECTI-
C				 VEMENT
C
C			NDR,NDT,NDF=DIMENSIONS DES TABLEAUX COMME DECLARE
C				    DAN LE PROGRAMME APPELLANT LA ROUTINE..
C			JCOEFF	= NOMBRE DES COEFFICIENTS QUI DOIVENT ETRE
C			          MIS A ZERO
C			IPAR	PARMETRE INDIQUANT LA PARITE:
C
C			IPAR=0 	AUCUNE PARITE
C			IPAR=2 	LE SCALAIRE EST SYMETRIQUE
C				 PAR RAPPORT L'INVERSION z -> -z
C			IPAR=3 	LE SCALAIRE EST ANTI-SYMETRIQUE
C				 PAIR PAR RAPPORT L'INVERSION z -> -z
C			IPAR=4 	LE SCALAIRE EST SYMETRIQUE
C				PAIR PAR RAPPORT L'INVERSION z -> -z ET PAR
C				L'INVERSION x,y -> -x,-y
C			IPAR=5 	LE SCALAIRE EST ANTI-SYMETRIQUE
C				 PAR RAPPORT L'INVERSION z -> -z ET PAR
C				L'INVERSION x,y -> -x,-y
C
C			INDD	= TABLEAU CONTENENT LE TYPE D'ECHANTILLONAGE
C				  DANS CHAQUE COQUILLE:
C			          INDD(LZON)=0 ECHANTILLONAGE RAREFIE
C				  INDD(LZON)=1,2,3 ECHANTILLONAGE DENSE
C			SCAL	= TABLEAU IMPUT A DESIALISER ET EN OUTPUT
C				  LE SCALAIRE DESIALISE
C
C
C $Id$
C $Log$
C Revision 1.2  2012/03/30 12:12:43  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.1  1997/11/28  13:52:40  eric
c Initial revision
c
C
C $Header$
C
C

	IMPLICIT NONE
	character*120 header
	data header/'$Header$'/

	INTEGER NDL,NDR,NDT,NDF,IPAR,JCOEF,INDD,NZON,NR1,NR,NT1,NF,NF1,
     1	LZON,N164,IND,LR,LT,LF,IPA,I314,IPAP,JCOEFF
C
	DOUBLE PRECISION SCAL,CC,VA0,VA1,ZA0,ZA1,COER,COEF
C
	PARAMETER(N164=164)
C
	DIMENSION NDL(*),INDD(*),SCAL(NDR,NDT,NDF,*),CC(N164)
	DIMENSION JCOEFF(*),COER(4,2,2),COEF(4,2,3)
C
	SAVE I314,COER,COEF
	DATA I314/0/
	IF(I314.NE.314) THEN
	I314=314
	DO LF=1,2
	DO LT=1,2
	DO LR=1,4
	COER(LR,LT,LF)=0
	ENDDO
	ENDDO
	ENDDO
C
	COER(1,1,1)=1
	COER(2,1,1)=-.5
	COER(1,2,1)=1
	COER(2,2,1)=.5
C
	DO LT=1,3
	DO LF=1,2
	DO LR=1,4
	COEF(LR,LF,LT)=0
	ENDDO
	ENDDO
	ENDDO
C
	COEF(1,1,1)=1
	COEF(2,1,1)=.5
	COEF(1,2,1)=1
	COEF(2,2,1)=-.5
	DO LT=2,3
	DO LF=1,2
	DO LR=1,4
	COEF(LR,LF,LT)=COEF(LR,LF,1)
	ENDDO
	ENDDO
	ENDDO
C
	COER(1,2,2)=-.75
	COER(2,2,2)=-.25
	COER(1,1,2)=-.25
	COER(2,1,2)=.25
	ENDIF
C
	NZON=NDL(1)
	NT1=NDL(NZON+2)
	NF=NDL(NZON+3)
	NF1=NF+1
C
	IF(IPAR.EQ.4.OR.IPAR.EQ.5) THEN
	IPAP=1
	IPA=0
	IF(IPAR.EQ.5) THEN
	IPA=1
	IPAP=2
	ENDIF
C
		IF(IPAR.EQ.5) THEN
		PRINT 101
		PRINT*,'LE CAS IPAR=',IPAR
		PRINT*,'N A PAS ETE TESTE, Sorry !'
		PRINT*,'(Mais il devrait marcher)'
		STOP
		ENDIF
C
	DO LZON=1,NZON
	NR1=NDL(LZON+1)
	JCOEF=NR1-JCOEFF(LZON)
	NR=NR1-1
	IND=INDD(LZON)
C
C		SI IPAR=2,4 ON IMPOSE QUE LA VALEUR DE LA FONCTION
C		EN r=0 SOIT CONSERVEE
C		SI IPAR=1,3,5 ON IMPOSE QUE LA VALEUR DE LA DERIVE
C		SE LA FONCTION SOIT CONSERVEE
C

	IF(JCOEF.LE.0) GO TO 777
	DO LF=1,NF
	DO LT=1,NT1
	DO LR=1,NR1
	CC(LR)=SCAL(LR,LT,LF,LZON)
	ENDDO
	IF(IND.EQ.0) THEN
C
	CALL EXTR1S(NR,0,IPA,IPA,CC,VA0)
	CALL EXTR1S(NR,1,0,IPA,CC,VA1)
C
	DO LR=JCOEF,NR1
	CC(LR)=0
	ENDDO
	CALL EXTR1S(JCOEF,0,IPA,IPA,CC,ZA0)
	CALL EXTR1S(JCOEF,1,0,IPA,CC,ZA1)
	ZA0=VA0-ZA0
	ZA1=VA1-ZA1
	DO LR=1,2
	CC(LR)=CC(LR)+ZA1*COER(LR,2,IPAP)+ZA0*COER(LR,1,IPAP)
	ENDDO
	ELSE
C
	CALL EXTM1S(NR,0,0,CC,VA0)
	CALL EXTM1S(NR,1,0,CC,VA1)
	DO LR=JCOEF,NR1
	CC(LR)=0
	ENDDO
	CALL EXTM1S(JCOEF,0,0,CC,ZA0)
	CALL EXTM1S(JCOEF,1,0,CC,ZA1)
	ZA0=VA0-ZA0
	ZA1=VA1-ZA1
	DO LR=1,2
	CC(LR)=CC(LR)+ZA0*COEF(LR,1,IND)+ZA1*COEF(LR,2,IND)
	ENDDO
	ENDIF
	DO LR=1,2
	SCAL(LR,LT,LF,LZON)=CC(LR)
	ENDDO
	DO LR=JCOEF,NR1
	SCAL(LR,LT,LF,LZON)=0
	ENDDO
	ENDDO
	ENDDO
  777	CONTINUE
	ENDDO
	ENDIF
C
  101	FORMAT(1X,' ')
	RETURN
	END
