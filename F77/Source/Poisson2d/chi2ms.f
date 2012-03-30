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

       SUBROUTINE CHI2MS(N,N64,F,CS,CC)

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
C      Tous les tests de routine ont ete effectues les 20/3/1986
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
c Revision 1.2  1997/05/23  11:33:12  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:40:23  hyc
C Initial revision
C
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

	save	N1,N2,N3,N21,X0,PI,SEN,N65,NDIM
	save	N63,N66,N6565,NM65,NM650,NM651,N62,NM652
	save	N6519,NM6,NM7,NM62,NM20,NM19,NFON

C
       N513=513
       IF(N.LT.N513) GO TO 12
  200  FORMAT(10X,'DIMENSIONS INSUFFISANTES DANS CHI2MS,N=',I4)
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
       N3=N2-1
       N21=N2+1
       X0=0
       PI=2*ACOS(X0)/N
       DO 11 L=2,N
       SEN(L)=.5+.25/SIN(PI*(L-1))
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
       N66=N65+1
       N6565=N65+N65
       NM65=N65*N
       NM650=NM65+1
       NM651=NM65+N65
       N62=N64-1
       NM652=N65*N21
       N6519=N65*(N21-2)+1
       NM6=NM650-N65
       NM7=NM65-N65
       NM62=NM6+N62
       NM20=N2*N65+1
       NM19=NM20-1
       NFON=N64
  18   CONTINUE
C
C
C           COMBINAISON LINEAIRE DES COEFFICIENTS.
C
       DO 1 M=1,N64
       CC(M)=-F(M+N65)*2
  1    CONTINUE
C
       N20L=N66
       N19L=N65
       N18L=0
       N16L=N6565
       N21L=N65+N64
       DO 2 L=2,N2
       L2=L+L
       DO 3 M=N20L,N21L
       CS(M)=F(M+N19L)
       CC(M)=F(M+N18L)-F(M+N16L)
  3    CONTINUE
       N20L=N20L+N65
       N21L=N21L+N65
       N19L=N19L+N65
       N18L=N18L+N65
       N16L=N16L+N65
  2    CONTINUE
       DO 5 M=1,N64
       CC(M+NM19)=F(M+NM7)*2
  5    CONTINUE
       CALL TFINMS(N,2,N64,CC,CS,F)
       N21L=NM652
       N20L=N6519
       DO 8 L=1,N3
       N2LX=N21L-N20L+1
       SEN20=SEN(N21-L)
       LX64=N20L+N63
       DO 9 LX=N20L,LX64
       CS(LX)=(F(LX+N2LX)+F(LX))*SEN20
   9   CONTINUE
       DO 13 LX=N20L,LX64
       CC(LX)=F(LX)-CS(LX)
       CC(LX+N2LX)=F(LX+N2LX)-CS(LX)
  13   CONTINUE
       N21L=N21L+N65
       N20L=N20L-N65
  8    CONTINUE
C
C
       DO 14 M=NM20,NM652
       CC(M)=-F(M)*.5
  14   CONTINUE
C
       DO 15 M=NM650,NM651
       CC(M)=0
  15   CONTINUE
       DO 16 M=1,N65
       CC(M)=0
  16   CONTINUE
  110  FORMAT(1X,2I5,5E16.4)
  100  FORMAT(1X,10D12.4)
101    FORMAT(1X,'CHEM64')
       RETURN
       END
C
