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

       SUBROUTINE CHIYMS(N,N64,F,CC,CS)

		implicit double precision (a-h,o-z)

C
C           ROUTINE POUR LA TRANSFORMATION INVERSE SIMULTANEE
C           DE N64 FONCTIONS. LA ROUTINE EST COMPLETEMENT
C           CRAYTINISEE.
C
C      ARGUMENTS DE LA ROUTINE:
C      
C      N    =NOMBRE DES DEGREES DE LIBERTEE-1. N DOIT ETRE
C            PAIR ET = A 2**p*3**q*5**m AVEC p,q,m NOMBRES
C           ENTIERS.
C      N64= NOMBRE DES FONCTION QUI DOIVENT ETRE TRANSFORMEES.
C      F    = COEFFICIENTS DE TCHEBYTCHEV DES FONCTIONS A TRANSFOR-
C            MER. LES COEFFICIENTS SONT STOCKES 'EN PARALLEL'
C           I.E. DANS F(1),F(2)....F(N64) IL-Y-A LE PREMIER
C           COEFFICIENT DES N64 FONCTIONS, DANS F(N64+1),F(N64+2),
C           ...F(N64+N64) LE 2ME COEFFICIET DES N64 FONCTIONS,
C           ET AINSI DE SUITE. SI N64 EST UN MULTIPLE DE 8, POUR
C           DE CRAYTINISATION LES
C           COEFFICIENTS SONT STOCKES DANS LA FACON SUIVANTE:
C           F(1),F(2),...F(N64) POUR LE PREMIER COEFFICIENT
C           F(N64+1+1),F(N64+1+2),...F(N64+1+N64) POUR LE 2ME
C           COEFFICIENT, ET AINSI DE SUITE.
C      CC   = TABLEAU DE TRAVAIL
C      CS   = TABLEAU OUTPUT: DANS CS IL-Y-A LE N1 VALEURES DES N64
C             FONCTIONS ('EN PARALLEL')
C             LES DIMENSIONS MINIMUM DES TABLEAUX
C             SONT (N+3)*(N64+1).
C            LE TABLEAU F EST DETRUIT.
C
C      CETTE ROUTINE EST SPECIALISEE: ELLE DOIT ETRE EMPLOYEE AVEC FOUCI
C  2S
C           FUCI3S...
C
C           Tous les test de routine ont ete executes l 8/10/85.
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:32:31  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:35:16  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

       DIMENSION CC(*),CS(*),F(*)
       DIMENSION SEN(513),F11(513),SOM1(513),FN21(513),SOMM(513)
       DATA NDIM/0/
       DATA NFON/0/

	save	N1,N2,N21,N3,X0,PI,SEN,NFON,NDIM,N65,NM65
	save	N651,NM650,NM655,NM623,N63,N65N65,NM20,NM201,NM2064

       N513=513
       IF(N.LT.N513) GO TO 12
       PRINT 200,N,N513
  200  FORMAT(10X,'DIMENSION INSUFFISANTES DAS CHINM64,N=',I4,
     , ' DIMENSIONS MAX.=',I4)
       CALL EXIT
C
  12   CONTINUE
C
C           PREPARATION DES QUANTITES NECESSAIRES POUR LE CALCUL.
C           CES QUANTITES SONT CALCULEES LA PREMIERE FOIS QU'ON
C           APPELLE LA ROUTINE ET TOUTES LES FOIS Q'ON CHENGE
C           LA VALEUR DE N.
C
       IF(N.EQ.NDIM) GO TO 10
       N1=N+1
       N2=N/2
       N21=N2+1
       N3=N2-1
       X0=0
       PI=2.*ACOS(X0)/N
       DO 11 L=2,N2
       SEN(N21-L)=.25/SIN(PI*(L-1))+.5
  11   CONTINUE
  10   CONTINUE
C
C           PREPARATION DES QUANTITES NECESSAIRES AU CLACUL.
C           CES QUANTITES SONT CALCULEES LA PREMIERE FOIS 
C           QU'ON APPELLE LA ROUTINE ET TOUTES LES FOIS 
C           QU'ON CHANGE N64.
C
       IF(NFON.EQ.N64.AND.NDIM.EQ.N) GO TO 20
       NFON=N64
       NDIM=N
       N65=N64
       IF((N64/8)*8.EQ.N64) N65=N64+1
       NM65=N65*N
       N651=N65+1
       NM650=NM65-N65
       NM655=NM65+N65
       NM623=(N2-1)*N65
       N63=N64-1
       N65N65=N65+N65
       NM20=N2*N65
       NM201=NM20+1
       NM2064=NM20+N64
  20   CONTINUE
C
C
C ****   CALCUL DE LA FONCTION EN TETA=0 ET TETA=PI
C
       DO 7 M=1,N64
       SOM1(M)=0
       SOMM(M)=0
  7    CONTINUE
C
       DO 8 L=N65N65,NM650,N65N65
       DO 30 M=1,N64
       SOMM(M)=SOMM(M)+F(M+L)
30     CONTINUE
   8   CONTINUE
       DO 5 L=N65,NM65,N65N65
       DO 31 M=1,N64
       SOM1(M)=SOM1(M)+F(M+L)
31     CONTINUE
  5    CONTINUE
C      
       DO 1000 M=1,N64
       CC(M)=(F(M)+F(M+NM65))*.5
1000	CONTINUE
C
       DO 32 M=1,N64
       F11(M)=SOMM(M)+SOM1(M)+CC(M)
       FN21(M)=SOMM(M)-SOM1(M)+CC(M)
32     CONTINUE
C
       DO 33 L=1,NM655
       CS(L)=0
33     CONTINUE
C
C           LA BOUCLE SUIVANTE EST EQUIVALENTE A:
C
       JL1=N651
       JL2=JL1+N63
       JL3=JL1+N65N65
       JL4=JL3+N63
       DO 1100 L=2,N,2
       DO 1101 M=JL1,JL2
       CS(M)=(CS(M)-F(M))
1101	CONTINUE
C
       DO 1102 M=JL3,JL4
       CS(M)=CS(M)+F(M-N65N65)
1102	CONTINUE
       JL1=JL1+N65N65
       JL2=JL1+N63
       JL3=JL3+N65N65
       JL4=JL3+N63
1100	CONTINUE
C           
       JL1=1
       JL2=JL1+N63
       DO 1200 L=1,N1,2
       DO 1201 M=JL1,JL2
       CS(M)=F(M)
1201	CONTINUE
       JL1=JL1+N65N65
       JL2=JL1+N63
1200	CONTINUE
C
       CALL TFIYMS(N,N64,CS,CC)
C
       LSEN=1
       DO 41 L=N65,NM623,N65
       N21L=NM201+L
       N20L=NM20-L
       L2=L+L
       SENN=SEN(LSEN)
       DO 40 M=N21L,N63+N21L
       F12=(CC(M-L2)-CC(M))*SENN
       CS(M)=CC(M)+F12
       CS(M-L2)=CC(M-L2)-F12
40     CONTINUE
       LSEN=LSEN+1
41     CONTINUE
C
       DO 42 M=1,N64
       CS(M)=F11(M)
       CS(M+NM65)=FN21(M)
42     CONTINUE
       DO 43 M=NM201,NM2064
       CS(M)=CC(M)
43     CONTINUE
C
100    FORMAT(1X,'CHINM64',10D12.4)
101    FORMAT(1X,'CHINM64')
 110   FORMAT(1X,'CHIM64',5D24.16)
 300   FORMAT(10X,'DANS CHINM64 LMAX=',I5)
  120  FORMAT(1X,'CHIM64',20I5)
       RETURN
       END
