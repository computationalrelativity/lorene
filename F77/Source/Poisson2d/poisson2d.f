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
	SUBROUTINE POISSON2D(NDL1,NDR,NDT,NDF,INDD,ERRE,SOUMAT,SOUQUAD,
     1			     ALAMB,POT)

C
C $Id$
C $Log$
C Revision 1.2  2002/03/25 09:16:59  m_bejger
C Increased the number of domains (NZOE) from 4 to 5
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.2  1998/07/20  12:43:52  eric
c Augmentation NDR0, NDT0
c
C Revision 1.1  1998/07/01 13:27:06  eric
C Initial revision
C
C
C $Header$
C
C

		IMPLICIT NONE

	INTEGER NDL1(*), NDR, NDT, NDF, INDD(*)
	REAL*8 ERRE(NDR,*), SOUMAT(NDR,NDT,NDF,*)
	REAL*8 SOUQUAD(NDR,NDT,NDF,*), ALAMB, POT(NDR,NDT,NDF,*)

	character*100 header
	data header/'$Header$'/

	INTEGER NDR0, NDT0, NDF0, NDZ0, N64
	INTEGER ND64Q, ND2Z, NDEQ
	PARAMETER (NDR0=130, NDT0=70, NDF0=4, NDZ0=5, N64=20)
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

	INTEGER IND, LZON, NR1, LR, LY, NZOE, NZON, NY1, NF, LZ, LT, LF
	REAL*8 AMAQ, AMAT, X0, X1, C1

C******************************************************************************

C 
C... Test de dimension
C
	IF (NDR.GT.NDR0) THEN
		WRITE(*,*) 'POISSON_2D: NDR trop grand !'
		WRITE(*,*) '  NDR, NDR0 : ', NDR, NDR0
		STOP
	ENDIF

	IF (NDT.GT.NDT0) THEN
		WRITE(*,*) 'POISSON_2D: NDT trop grand !'
		WRITE(*,*) '  NDT, NDT0 : ', NDT, NDT0
		STOP
	ENDIF

	IF (NDF.GT.NDF0) THEN
		WRITE(*,*) 'POISSON_2D: NDF trop grand !'
		WRITE(*,*) '  NDF, NDF0 : ', NDF, NDF0
		STOP
	ENDIF



C
C... Recuperation information nombre de points
C
	NZOE = NDL1(1)
	IF (NZOE.GT.NDZ0) THEN
		WRITE(*,*) 'POISSON_2D: NZOE trop grand !'
		WRITE(*,*) '  NZOE, NDZ0 : ', NZOE, NDZ0
		STOP
	ENDIF

	NZON = NZOE - 1
	NY1 = NDL1(NZOE+2)
	NF = NDL1(NZOE+3)

	WRITE (*,*) '------------- POISSON2D -------------------'
C	WRITE (*,*) 'NZOE, NZON : ', NZOE, NZON	
C	WRITE (*,*) 'NY1, NF : ', NY1, NF
C		DO LZON = 1, NZOE
C	WRITE (*,*) 'LZ : ', LZON, '   NR : ', NDL1(1+LZON)
C		ENDDO
	
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
C Les termes quadratiques sont contenus dans SOUQUAD; on les passe dans l'espace
C  des coefficients en theta et en r:
C
C
	CALL FCEZ3S(NDL1,NDR,NDT,NDF,N64,2,0,C64,CC,CS,DEN2,DEN1,SOUQUAD)
	CALL FCEZ3S(NDL1,NDR,NDT,NDF,N64,1,0,C64,CC,CS,DEN2,DEN1,SOUQUAD)
C
C Terme "monopolaire" l=0 du developpement de SOUQUAD en cos(l*theta) ---> TRA0:
C 
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	TRA0(LR,LZON)=SOUQUAD(LR,1,1,LZON)
	ENDDO
	ENDDO

C  Calcul de la masse generee par le termes quadratiques
C
C Calcul de   AMAQ = int_0^{+infini} SOUQUAD(r)(l=0) r dr   [cf. Eq.(4.6)] 
C
	IND=2
	CALL GPAR2S(NDL1,NDR0,IND,C64,ERRE0,TRA0,TRA1)

C	PRINT*,'TRA1=',TRA1(1,1),TRA1(2,1)
	AMAQ=TRA1(1,1)+TRA1(2,1)

	WRITE (*,*) 'Integrale(l=0) de SOUQUAD : ', REAL(TRA1(1,1)), 
     1		REAL(TRA1(2,1)), REAL(AMAQ)

C
C
C Les termes matiere sont contenus dans SOUMAT
C  Passage en Tchebyshev en theta et en r:
C
	CALL FCEZ3S(NDL,NDR,NDT,NDF,N64,2,0,C64,CC,CS,DEN2,DEN1,SOUMAT)
	CALL FCEZ3S(NDL,NDR,NDT,NDF,N64,1,0,C64,CC,CS,DEN2,DEN1,SOUMAT)
C
C Terme "monopolaire" l=0 du developpement de SOUMAT en cos(l*theta) ---> TRA2:
C 
	DO LZON=1,NZON
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	TRA2(LR,LZON)=SOUMAT(LR,1,1,LZON)
	ENDDO
	ENDDO
C
C Calcul de   AMAT = int_0^{+infini} SOUMAT(r)(l=0) r dr   [cf. Eq.(4.6)] 

	IND=1
	CALL GPAR2S(NDL,NDR0,IND,C64,ERRE0,TRA2,TRA3)
	AMAT=TRA3(1,1)

	WRITE (*,*) 'Integrale(l=0) de SOUMAT : ', REAL(AMAT)

C
C----------------------------------------
C Identite du viriel GRV2 [ Eq.(4.6), Eq.(6.36) ] 
C----------------------------------------
C   ALAMB = Abs(lambda), lambda etant defini par l'Eq.(6.38)
C    ALAMB doit etre egal a 1 pour une solution exacte
C
	X0=0
	X1=1 
	ALAMB=0

C##	IF(AMAQ.NE.X0) ALAMB=ABS(AMAT/AMAQ)
	IF(AMAQ.NE.X0) ALAMB = - AMAT / AMAQ 
C
C	WRITE (*,*) 'ALAMB, 1-ALAMB : ', ALAMB, REAL(X1-ALAMB)
C
C Preparation de la source totale pour dzeta suivant l'Eq.(6.37):
C
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LY=1,NY1
	DO LR=1,NR1
		POT(LR,LY,1,LZON) = ALAMB * SOUQUAD(LR,LY,1,LZON)
	ENDDO
	ENDDO
	ENDDO
C
C ... On ajoute la matiere dans les zones non compactifiees : 
	DO LZON=1,NZON
	NR1=NDL1(LZON+1)
	DO LY=1,NY1
	DO LR=1,NR1
		POT(LR,LY,1,LZON) = POT(LR,LY,1,LZON) + SOUMAT(LR,LY,1,LZON)
	ENDDO
	ENDDO
	ENDDO

	IND=4
C
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	TRA0(LR,LZON)=POT(LR,1,1,LZON)
	ENDDO
	ENDDO
C
	CALL GR2P1S(NDL1,NDR0,2,ERRE0,TRA2,TRA0)
C
	CALL GR2P3S(NDL1,NDR,NDT,NDF,INDD,0,CC,C64,BB,DEN1,
     1		    DEN2,TRAB0,TRAB1,ERRE,SOLHH,POT)
C
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	POT(LR,1,1,LZON)=TRA0(LR,LZON)
	ENDDO
	ENDDO
C
C
	CALL FCIQ3S(NDL1,NDR,NDT,NDF,INDD,1,0,0,2,C64,CC,CS,DEN1, 
     1	DEN2,ERRE,POT)
C
C
C Constante d'integration t.q. valeur nulle a l'infini
C
	NR1=NDL1(NZOE+1)
	C1=POT(NR1,1,1,NZOE) ! valeur du coef l=0 de dzeta a l'infini 
C
	DO LZON=1,NZOE
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	POT(LR,1,1,LZON)=POT(LR,1,1,LZON)-C1
	ENDDO
	ENDDO

CC
C Retour dans l'espace des configurations en theta (on y etait deja en r):
C
	CALL FCIQ3S(NDL1,NDR,NDT,NDF,INDD,2,0,0,2,C64,CC,CS,DEN1,
     1	DEN2,ERRE,POT)

	WRITE (*,*) '-------------------------------------------'

	RETURN 

	END

