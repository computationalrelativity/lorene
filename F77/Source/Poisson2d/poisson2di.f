C
C   Copyright (c) 1998 Eric Gourgoulhon
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
	SUBROUTINE POISSON2DI(NDL1,NDR,NDT,NDF,INDD,ERRE,SOU,POT)

C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1998/07/20  12:44:07  eric
c Augmentation NDR0, NDT0
c
C Revision 1.1  1998/07/01 13:24:24  eric
C Initial revision
C
C
C $Header$
C
C

		IMPLICIT NONE

	INTEGER NDL1(*), NDR, NDT, NDF, INDD(*)
	REAL*8 ERRE(NDR,*)
	REAL*8 SOU(NDR,NDT,NDF,*), POT(NDR,NDT,NDF,*)

	character*100 header
	data header/'$Header$'/

	INTEGER NDR0, NDT0, NDF0, NDZ0, N64
	INTEGER ND64Q, ND2Z, NDEQ
	PARAMETER (NDR0=130, NDT0=70, NDF0=4, NDZ0=4, N64=20)
C##	PARAMETER (ND2Z=MAX(NDZ0,NDF0,8), NDEQ=NDZ0+8)
	PARAMETER (ND2Z=8, NDEQ=NDZ0+8)
	PARAMETER (ND64Q=(NDR0+2)*(NDT0+2)*NDF0)
	REAL*8 CC(ND64Q), CS(ND64Q), C64(ND64Q)

	INTEGER NDL(NDEQ)

	REAL*8 TRA0(NDR0,NDT0),TRA1(NDR0,NDT0),TRA2(NDR0,NDT0)
	REAL*8 TRA3(NDR0,NDT0)
	REAL*8 TRAB0(NDR0,NDT0,ND2Z), TRAB1(NDR0,NDT0,ND2Z)
	REAL*8 DEN1(NDR0,NDT0,ND2Z), DEN2(NDR0,NDT0,ND2Z)
	REAL*8 BB(NDR0,12), ERRE0(NDR0,ND2Z)
 	REAL*8 SOLHH(NDR0,NDT0,8,NDZ0)

	INTEGER LZON, NR1, LR, LY, NZOE, NZON, NY1, NF, LZ, LT, LF

C******************************************************************************

C 
C... Test de dimension
C
	IF (NDR.GT.NDR0) THEN
		WRITE(*,*) 'POISSON_2DI: NDR trop grand !'
		WRITE(*,*) '  NDR, NDR0 : ', NDR, NDR0
		STOP
	ENDIF

	IF (NDT.GT.NDT0) THEN
		WRITE(*,*) 'POISSON_2DI: NDT trop grand !'
		WRITE(*,*) '  NDT, NDT0 : ', NDT, NDT0
		STOP
	ENDIF

	IF (NDF.GT.NDF0) THEN
		WRITE(*,*) 'POISSON_2DI: NDF trop grand !'
		WRITE(*,*) '  NDF, NDF0 : ', NDF, NDF0
		STOP
	ENDIF



C
C... Recuperation information nombre de points
C
	NZOE = NDL1(1)
	IF (NZOE.GT.NDZ0) THEN
		WRITE(*,*) 'POISSON_2DI: NZOE trop grand !'
		WRITE(*,*) '  NZOE, NDZ0 : ', NZOE, NDZ0
		STOP
	ENDIF

	NZON = NZOE - 1
	NY1 = NDL1(NZOE+2)
	NF = NDL1(NZOE+3)
	
C
C... Tableau NDL decrivant les zones non compactifies:
C
	NDL(1) = NZON 
		DO LZON = 1, NZON
	NDL(1+LZON) = NDL1(1+LZON)
		ENDDO
	NDL(NZON+2)=NY1
	NDL(NZON+3)=NF	
	NDL(NZON+4)=NDL1(NZOE+4)

C
C... Tableau ERRE0 decrivant les points de collocation en r dimensionne comme
C		(NDR0,*), contrairement a ERRE qui est dimensionne comme (NDR,*)
C
		DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
			DO LR=1,NR1
	ERRE0(LR,LZON)=ERRE(LR,LZON)
			ENDDO
	ENDDO

C
C Mise a zero des tableaux de travail
C
	DO LR = 1, ND64Q
		CC(LR) = 0
		CS(LR) = 0 
		C64(LR) = 0
	ENDDO

	DO LT = 1, NDT0
		DO LR = 1, NDR0
			TRA0(LR,LT) = 0
			TRA1(LR,LT) = 0
			TRA2(LR,LT) = 0
			TRA3(LR,LT) = 0
		ENDDO
	ENDDO
	
	DO LZ = 1, ND2Z
		DO LT = 1, NDT0
			DO LR = 1, NDR0
				TRAB0(LR,LT,LZ) = 0
				TRAB1(LR,LT,LZ) = 0
				DEN1(LR,LT,LZ) = 0
				DEN2(LR,LT,LZ) = 0
			ENDDO
		ENDDO
	ENDDO

	DO LZ = 1, 12
		DO LR = 1, NDR0
			BB(LR,LZ) = 0 
		ENDDO
	ENDDO

	DO LZ = 1, NDZ0
		DO LF = 1, 8
			DO LT = 1, NDT0
				DO LR = 1, NDR0
					SOLHH(LR,LT,LF,LZ) = 0
				ENDDO
			ENDDO
		ENDDO
	ENDDO


C
C Calcul des coefficients de la source :
C
C
	CALL FCEZ3S(NDL,NDR,NDT,NDF,N64,2,1,C64,CC,CS,DEN2,DEN1,SOU)
	CALL FCEZ3S(NDL,NDR,NDT,NDF,N64,1,0,C64,CC,CS,DEN2,DEN1,SOU)


C
C Preparation de la source pour GR2P3S
C
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LY=1,NY1
	DO LR=1,NR1
		POT(LR,LY,1,LZON) = SOU(LR,LY,1,LZON)
	ENDDO
	ENDDO
	ENDDO
C
C Appel GR2P3S
C
	CALL GR2P3S(NDL1,NDR,NDT,NDF,INDD,0,CC,C64,BB,DEN1,
     1		    DEN2,TRAB0,TRAB1,ERRE,SOLHH,POT)

C
C Retour dans l'espace des configurations
C
	CALL FCIQ3S(NDL1,NDR,NDT,NDF,INDD,1,0,0,2,C64,CC,CS,DEN1, 
     1	DEN2,ERRE,POT)
	CALL FCIQ3S(NDL1,NDR,NDT,NDF,INDD,2,1,0,2,C64,CC,CS,DEN1,
     1	DEN2,ERRE,POT)


	RETURN 

	END

