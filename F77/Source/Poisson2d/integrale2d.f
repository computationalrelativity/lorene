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
	SUBROUTINE INTEGRALE2D(NDL1,NDR,NDT,NDF,ERRE,SOU,RESU)

C
C $Id$
C $Log$
C Revision 1.3  2003/09/11 16:04:49  j_novak
C Changed the maximal size in r to 134 in integral2d.f
C
C Revision 1.2  2002/03/25 09:16:58  m_bejger
C Increased the number of domains (NZOE) from 4 to 5
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.2  1998/07/20  12:44:19  eric
c Augmentation NDR0, NDT0
c
C Revision 1.1  1998/07/06 15:54:24  eric
C Initial revision
C
C
C $Header$
C
C
		IMPLICIT NONE

	INTEGER NDL1(*), NDR, NDT, NDF
	REAL*8 ERRE(NDR,*), SOU(NDR,NDT,NDF,*), RESU

	character*100 header
	data header/'$Header$'/

	INTEGER NDR0, NDT0, NDF0, NDZ0, N64
	INTEGER ND64Q, ND2Z
	PARAMETER (NDR0=135, NDT0=70, NDF0=4, NDZ0=5, N64=20)
C##	PARAMETER (ND2Z=MAX(NDZ0,NDF0,8), NDEQ=NDZ0+8)
	PARAMETER (ND2Z=8)
	PARAMETER (ND64Q=(NDR0+2)*(NDT0+2)*NDF0)
	REAL*8 CC(ND64Q), CS(ND64Q), C64(ND64Q)

	REAL*8 TRA0(NDR0,NDT0),TRA1(NDR0,NDT0)
	REAL*8 DEN1(NDR0,NDT0,ND2Z), DEN2(NDR0,NDT0,ND2Z)
	REAL*8 ERRE0(NDR0,ND2Z)

	INTEGER IND, LZON, NR1, LR, NZOE, NY1, NF, LZ, LT

C******************************************************************************

C 
C... Test de dimension
C
	IF (NDR.GT.NDR0) THEN
		WRITE(*,*) 'INTEGRALE2D: NDR trop grand !'
		WRITE(*,*) '  NDR, NDR0 : ', NDR, NDR0
		STOP
	ENDIF

	IF (NDT.GT.NDT0) THEN
		WRITE(*,*) 'INTEGRALE2D: NDT trop grand !'
		WRITE(*,*) '  NDT, NDT0 : ', NDT, NDT0
		STOP
	ENDIF

	IF (NDF.GT.NDF0) THEN
		WRITE(*,*) 'INTEGRALE2D: NDF trop grand !'
		WRITE(*,*) '  NDF, NDF0 : ', NDF, NDF0
		STOP
	ENDIF

C
C... Recuperation information nombre de points
C
	NZOE = NDL1(1)
	IF (NZOE.GT.NDZ0) THEN
		WRITE(*,*) 'INTEGRALE2D: NZOE trop grand !'
		WRITE(*,*) '  NZOE, NDZ0 : ', NZOE, NDZ0
		STOP
	ENDIF

	NY1 = NDL1(NZOE+2)
	NF = NDL1(NZOE+3)
	
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
		ENDDO
	ENDDO
	
	DO LZ = 1, ND2Z
		DO LT = 1, NDT0
			DO LR = 1, NDR0
				DEN1(LR,LT,LZ) = 0
				DEN2(LR,LT,LZ) = 0
			ENDDO
		ENDDO
	ENDDO


C
C Passage dans l'espace des coefficients en theta et en r:
C
	CALL FCEZ3S(NDL1,NDR,NDT,NDF,N64,2,0,C64,CC,CS,DEN2,DEN1,SOU)
	CALL FCEZ3S(NDL1,NDR,NDT,NDF,N64,1,0,C64,CC,CS,DEN2,DEN1,SOU)
C
C Terme "monopolaire" l=0 du developpement de SOU en cos(l*theta) ---> TRA0:
C 
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	TRA0(LR,LZON)=SOU(LR,1,1,LZON)
	ENDDO
	ENDDO

C
C Calcul de    int_0^{+infini} SOU(r)(l=0) r dr   
C
	IND=2
	CALL GPAR2S(NDL1,NDR0,IND,C64,ERRE0,TRA0,TRA1)

	RESU=TRA1(1,1)+TRA1(2,1)

	WRITE (*,*) 'INTEGRALE2D: Integ(l=0): ', REAL(TRA1(1,1)),
     1		'+',REAL(TRA1(2,1)),' = ', REAL(RESU)

	RETURN 

	END

