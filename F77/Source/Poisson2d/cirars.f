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
       SUBROUTINE CIRARS(N,INDEX,F,CC,CS)

		implicit double precision (a-h,o-z)

C
C **** INVERSION DE  LA TF DE TCHEBYTCHEV D'UNE FONCTION ECHANTILLONEE
C      RAREMENT A L'ORIGINE. (INTERVALLE D'ECHANTILLONAGE [0-1]
C **** DIMENSION MINIMALE DE TABLEUX=N+1
C **** ROUTINE DISTRUPTIVE9
C
C           ARGUMENTS DE LA ROUTINE:
C
C      N    = NOMBRE DES DEGRES DE LIBERTE'
C      INDEX = PARAMETRE: IL DOIT ETRE =0 SI LA FONCTION D'ENTREE
C             EST PAIRE, =1 SI IMPAIRE.
C      F    = IMPUT: TABLEAU CONTENANT LES COEFFICIENT DE TCHE-
C           BITCHEV DE LA FONCTION A INVERSER.
C      CC   = TABLEAU DE TRAVAIL.
C      CS   = OUTPUT.
C      
C           ROUtine testee le 20/41985.
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
c Revision 1.1  1997/10/21  14:08:13  eric
c Initial revision
c
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/


       DIMENSION CC(*),CS(*),F(*),SEN(1025),COSE(1025)

	data	ncontr/0/

	save	NCONTR,N1,N11,N2,N21,N3,X0,PI,PI2,X5,SEN,COSE

       IF(NCONTR.EQ.N) GO TO 999
       IF(N.LT.1025) GO TO 9
       PRINT 200,N
  200  FORMAT(10X,'DIMENSIONS INSUFFISANTES DAN CHINPR:N=',I5)
       CALL EXIT
  9    CONTINUE
       NCONTR=N
       N1=N+1
       N11=N+2
       N2=N/2
       N21=N2+1
       N3=N2-1
       X0=0
       PI=2*ACOS(X0)/N
       PI2=PI/2
       X5=.5
       SEN(1)=.25/SIN((N2-1)*PI)+X5
       DO 10 L=2,N3
       SEN(L)=.25/SIN((N2-L)*PI)+X5
       COSE(L)=1/SIN((L-1)*PI2)
  10   CONTINUE
       DO 1000 L=N21,N1
       SEN(L)=.25/SIN((N2-L)*PI)+X5
       COSE(L)=1/SIN((L-1)*PI2)
1000	CONTINUE
       COSE(N2)=1/SIN(N3*PI2)
  999  CONTINUE
C
C
       IF(INDEX.EQ.0) GO TO 11
C
       CS(1)=-F(1)
       F(N1)=-F(N)
C
       DO 12 L=2,N
       CS(L)=-(F(L)+F(L-1))*.5D+00
  12   CONTINUE
       DO 13 L=1,N
       F(L)=CS(L)
  13   CONTINUE
  11   CONTINUE
C
C ****   CALCUL DE LA FONCTION EN TETA=0 ET TETA=PI
C
       SOMM=0
       DO 4 L=3,N,2
       SOMM=SOMM+F(L)
4      CONTINUE
       SOM1=0
       DO 5 L=2,N ,2
       SOM1=SOM1+F(L)
5      CONTINUE
C
       FC=   (F(1)+F(N1))*.5
       F11=SOMM+SOM1    +FC
       FN21=SOMM-SOM1    +FC
C
C ****  CALCUL DE COEF. DE FOURIER DE LA FONCYION PONDEREE A PARTIR DE C
C ****  DE TCHEBYCHEV.
C
       CS(1)=0
       DO 2 L=1,N2
       FL=F(L+L)
       CS(L)=CS(L)+FL
       CS(L+1)=-FL
       CC(L)=F(L+L-1)
2      CONTINUE
       CC(N21)=F(N1)
       CS(1)=0
       CS(N21)=0
C
C ***    INVERSION POUR OBTENIR LA FONCTION PONDEREE.
C
       CALL TFINS(N,CC,CS,F)
C
C ****   RESTITUTION DE LA FONCTION
C
       DO 3 L=1,N3
       N21L=N21+L
       N20L=N21-L
       F12=(F(N20L)-F(N21L))*SEN(L)
       CS(N20L)=F(N21L)+F12
       CS(N21L)=F(N20L)-F12
  3    CONTINUE
C
       CS(1)=FN21
       CS(N21)=F(N21)
       CS(N1)= F11
C
       IF(INDEX.EQ.0) RETURN
       DO 14 L=2,N1
       CS(L)=CS(L)*COSE(L)
  14   CONTINUE
       CS(1)=0
100    FORMAT(10X,10D12.4)
101    FORMAT(1X,'IN')
       RETURN
       END
