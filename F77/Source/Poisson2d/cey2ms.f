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
       SUBROUTINE CEY2MS(N,N64,F,CS,CC)

		implicit double precision (a-h,o-z)

C
C  ***     ROUTINE POUR LE CALCUL SIMULTANE DES COEFF. DE TCHE-
C  ***     BYTCHEF DU 2ME GENRE DE N64 FONCTIONS.
C           LA ROUTINE EST COMPLETEMENT CRAYTINISEE.
C  ****    F DOIT ETRE ECHANTILLONEE DANS N+1 POINTS.. DIMENSION MINIMAL
C ***  DES TABLEAUX: =(N64+1)*(N+3), OU N+1 EST LE NOMBRE DES POINTS
C ***  D'ECHANTILLONAGE.
C
C           ARGUMENTS DE LA ROUTINE:
C
C      N    = NOMBRE DES DEGRES DE LIBERTE-1.
C      N64  = NOMBRE DES FONCTIONS D'ONT ON VEUT CALCULER LA TRANS-
C           FORMATION DE TCHEBYTCHEV. 
C      F    =TABLEAU IMPUT. LES VALEURS ECHANTILLONNES DES N64
C            FONCTIONS SONT STOCKEES EN 'PARALLEL': DANS F(1),F(2)..
C            F(N64) ON STOCKE LE VALEURES DE N64 FONCTIONS EN X=X0
C            DANS F(1+N64),F(2+N64),...F(N64+N64), LES VALEURES
C            DES N64 FONCTIONS DANS LE POINT X1 ET AINSI DE SUITE.
C              SI N64 EST UN MULTIPLE DE 8, POUR DES RAISONS DE
C            CRAYTINISATION LES DONNES SONT STOCKEES DAN LA FA-
C            CON SUIVANTE: F(1),F(2),...F(N64), F(N64+1+1),
C            F(N64+1+2),...F(N64+1+N64).
C      CS   =TABLEAU DE TRAVAIL
C      CC   =OUTPUT. LES COEFF. DE TCHEBYTCHEV SONT STO-
C                EN 'PARALLEL'.
C
C      Tous les tests de routine ont ete effectues les 20/3/1986.
C
C
C $Id$
C $Log$
C Revision 1.2  2012/03/30 12:12:42  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.1  1997/10/23  08:34:54  eric
c Initial revision
c
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

       DIMENSION CC(*),CS(*),F(*)
       DIMENSION SEN(513)
       DATA NDIM/0/
       DATA NFON/0/

	save	N1,N2,N21,X0,PI,SEN,N65,NDIM,N63,N6565,NM65,NM650,NM651
	save	N62,NM652,N6519,NM6,NM653,NM62,NM20,NFON
C
       N513=513
       IF(N.LT.N513) GO TO 12
  200  FORMAT(10X,'DIMENSIONS INSUFFISANTES DANS CEY2MS,N=',I4)
       PRINT 200,N,N513
       CALL EXIT
C
C           PREPARATION DES QUANTITES NECESSAIRES POUR LE CALCUL.
C           CES QUANTITES SONT CALCULEES LA
C           PREMIERE FOIS QU'ON APPELLE LA ROUTINE.
C           ELLES SONT RECALCULEES TOUTES LES FOIS QU'ON CHANGE
C           N.
C
  12   CONTINUE
       IF(N.EQ.NDIM) GO TO 10
       N1=N+1
       N2=N/2
       N21=N2+1
       X0=0
       PI=2*ACOS(X0)/N
       DO 11 L=1,N1
       SEN(L)=SIN((L-1)*PI)-.5
  11   CONTINUE
  10   CONTINUE
       IF(NFON.EQ.N64.AND.N.EQ.NDIM) GO TO 18
C
C           PREPARATION DES QUANTITES NECESSAIRES POUR LE CALCUL.
C           CES QUANTITES SONT CALCULEES LA PREMIERE FOIS QU'ON
C           APPELLE LA ROUTINE. ELLES SONT RECALCULEES TOUTES
C           LES FOIS QU'ON CHANGE N64.
C
       N65=N64
       NDIM=N
       IF((N64/8)*8.EQ.N64) N65=N64+1
       N63=N65-1
       N6565=N65+N65
       NM65=N65*N
       NM650=NM65+1
       NM651=NM65+N65
       N62=N64-1
       NM652=N65*N21
       N6519=N65*(N21-2)+1
       NM6=NM650-N65
       NM653=NM6-N6565
       NM62=NM6+N62
       NM20=N2*N65+1
       NFON=N64
  18   CONTINUE
C
C
C           PONDESRATION DES LA PARTIE ANTISIMMETRIQUE DES FONC-
C           TIONS.
C
       N21L=NM652
       N20L=N6519
       DO 1 L=1,N2
       N2LX=N21L-N20L+1
       SEN20L=SEN(N21-L)
       LX64=N20L+N62
       N22L=N20L+N2LX
       N23L=N22L+N62
C
       DO 6 LX=N20L,LX64
       CS(LX)=(F(LX)+F(LX+N2LX))*SEN20L
  6    CONTINUE
C
       DO 20 LX=N20L,LX64
       CC(LX)=F(LX)+CS(LX)
  20   CONTINUE
C
       DO 21 LX=N22L,N23L
       CC(LX)=F(LX)+CS(LX-N2LX)
  21   CONTINUE   
       N21L=N21L+N65
       N20L=N20L-N65
1      CONTINUE
C
       JM1=NM20
       JM2=JM1+N62
       DO 22 M=JM1,JM2
       CC(M)=F(M)*2
  22   CONTINUE   
C
       DO 3 M=NM20,NM652
C      CC(M)=F(M)*2
  3    CONTINUE
C
C           CALCUL DE LA T.F. FONCTION PRECEDENTE.
C
       CALL TFMYS(N,N64,CC,F)
C
C           RERAREGEMENT DES COEFFICIENTS
C
       DO 4 M=NM6,NM62
       CC(M)=-F(M+N65)*.5
  4    CONTINUE
C      
       JM1=NM653
       JM2=JM1+N62
       DO 24 L=1,N-2,2
       DO 23 M=JM1,JM2
       CC(M)=CC(M+N6565)-F(M+N65)
  23   CONTINUE   
       JM1=JM1-N6565
       JM2=JM1+N62
  24   CONTINUE   
C
       DO 25 M=1,N65
       CC(M+N65)=F(M)*.5
  25   CONTINUE   
C
       JM1=N6565+1
       JM2=JM1+N64-1
       DO 26 L=3,N,2
       DO 27 M=JM1,JM2
       CC(M)=F(M+N65)
  27   CONTINUE   
       JM1=JM1+N6565
       JM2=JM1+N64-1
  26   CONTINUE   
C
       DO 28 M=1,N65
       CC(M)=0
  28   CONTINUE   
C
       DO 29 M=NM650,NM650+N62
       CC(M)=0
  29   CONTINUE   
C
  110  FORMAT(1X,3I5,5E12.4)
  100  FORMAT(1X,10E12.4)
101    FORMAT(1X,'CHEM64')
       RETURN
       END
C
