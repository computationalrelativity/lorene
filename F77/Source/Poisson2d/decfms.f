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

       SUBROUTINE DECFMS(N,N64,CC,CS,Y)

		implicit double precision (a-h,o-z)

       DIMENSION CC(*),CS(*),Y(*)
C
C
C           ROUTINE POUR LE CALCUL DE L'OPERATEUR
C
C           COS(TETA)/SIN(TETA)*D/DTETA 
C
C           POUR LES FONCTIONS DEVELOPEES EN POLYNOMES DE TCHEBYTCHEV
C           DU 1ME TYPE.
C
C      ARGUMENTS DE LA ROUTINE:
C
C           N   =NOMBRE DE DEGRES DE LIBERTE'-1
C           N64 =PARAMETRE DEFINISSANT LE NOMBRE DES FONCTIONS QUI
C                DOIVENT ETRE TRAITEES SIMULTANEMENT
C           CC  =TABLEAU IMPUT CONTENANT LES COEFFICIENTS DE TCHBYTCHEV
C               DES N64 FONCTION A TRAITER STOCKES EN PARALEL,
C               (CC(1),CC(2),...CC(N64) SONT LE PREMIERS COEFFICIENTS
C               DES N64 FONCTIONS, CC(N64+1),CC(N64+2),...CC(N64+N64)
C               LES 2EMES COEFFICIENTS DE N64 FONCTIONS ET AINSI DE 
C               SUITE.
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:30:32  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:34:41  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

       DATA NDIM,NEQ/0,0/

	save	N0,N3,N65,N63,N66,N366,N266,N166,NM66
	save	N6565,N6566,NDIM,NEQ

C
C           INITIALISATION
C
       IF(NDIM.NE.N) THEN
       N0=N-1
       N3=N-3
       ENDIF
C
       IF(NDIM.EQ.N.AND.NEQ.EQ.N64) GO TO 800
       N65=N64
       IF((N64/8)*8.EQ.N64) N65=N64+1
       N63=N64-1
       N66=N65+1
       N366=N3*N65+1
       N266=N366+N65
       N166=N266+N65
       NM66=N166+N65
       N6565=N65+N65
       N6566=N6565+1
       NDIM=N
       NEQ=N64
  800  CONTINUE
C
C           MULTIPLICATION DES COEFFICIENTE DES FONCTION A TRAITER 
C           PAR LA MATRICE DE L'OPERATEUR EN QUESTION (MATRICE 
C           A 2 DIAGONALES APRES UN'ASTUCIEUSE COMBINAISON LINEAIRE
C           DES LIGNES)
C
C           MULTIPLICATION PARTIE PAIRE DE LA FONCTION
C
       DO 1000 M=1,N64
       CS(M)=-CC(M+N6565)*2
1000	CONTINUE
C
       JM1=N6566
       JM2=JM1+N63
       AL1=-2
       AL2=-4
       DO 1 L=3,N-3,2
       DO 2 M=JM1,JM2
       CS(M)=CC(M)*AL1+CC(M+N6565)*AL2
  2    CONTINUE 
       JM1=JM1+N6565
       JM2=JM1+N63
       AL1=AL2
       AL2=AL2-2
  1    CONTINUE 
C
C      MULTIPLICATION PARTIE IMPAIRE
C
       JM1=N66
       JM2=JM1+N63
       DO 1100 M=JM1,JM2
       CS(M)=-(CC(M)+CC(M+N6565)*3)
1100	CONTINUE
C
       JM1=JM1+N6565
       JM2=JM1+N63
       AL1=-3
       AL2=-5
       DO 1200 L=4,N-2,2
       DO 1201 M=JM1,JM2
       CS(M)=CC(M)*AL1+CC(M+N6565)*AL2
1201	CONTINUE
       JM1=JM1+N6565
       JM2=JM1+N63
       AL1=AL2
       AL2=AL2-2
1200	CONTINUE
C
C           INVERSION DE L'ASTUCIEUSE COMBINAISON LINEAIRE
C
       AL1=-N
       DO 1300 M=NM66,NM66+N63
       Y(M)=CC(M)*AL1
1300	CONTINUE
C
       AL1=-N0
       DO 1400 M=N166,N166+N63
       Y(M)=CC(M)*AL1
1400	CONTINUE
C
       AL1=-N+2
       AL2=1
       DO 1500 M=N266,N266+N63
       Y(M)=CC(M)*AL1+Y(M+N6565)
1500	CONTINUE
C      
       JM1=N366
       JM2=JM1+N63
       DO 6 L=1,N-2
       DO 7 M=JM1,JM2
       CC(M)=Y(M+N6565)+CS(M)
  7    CONTINUE 
       DO 8 M=JM1,JM2
       Y(M)=CC(M)
  8    CONTINUE 
       JM1=JM1-N65
       JM2=JM1+N63
  6    CONTINUE 
C
  100  FORMAT(1X,10E12.4)
  101  FORMAT(1X,' ')
       RETURN
       END
