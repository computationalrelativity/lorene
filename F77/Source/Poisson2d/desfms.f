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

       SUBROUTINE DESFMS(N,N64,CC,CS,Y)

		implicit double precision (a-h,o-z)

       DIMENSION CC(*),CS(*),Y(*)
C
C           ROUTINE POUR LE CALCUL DE L'OPERATEUR
C
C           COS(TETA)/SIN(TETA)*D/DTETA -1/SIN(TETA)**2
C
C           POUR LES FONCTIONS DEVELOPEES EN SERIE DE FOURIER
C           DU 2ME TYPE.
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
C           CS  =TABLEAU DE TYRAVAIL
C           Y   = TABLEAU OUTPUT CONTENENT LES COEFFICIENTS DE L'IMPUT
C                 APRES APPLICATION DE L'OPERATEUR
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:30:13  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:34:29  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

       DATA NDIM,NEQ/0,0/

	save	N0,N2,N3,N4,N65,N63,N66,N466,N366,N266,N166
	save	NM66,N6565,AN2,AN1,NDIM,NEQ

C
C           INITIALISATION
C
       IF(NDIM.NE.N) THEN
       N0=N-1
       N2=N-2
       N3=N-3
       N4=N-4
       ENDIF
C
       IF(NDIM.EQ.N.AND.NEQ.EQ.N64) GO TO 800
       N65=N64
       IF((N64/8)*8.EQ.N64) N65=N64+1
       N63=N64-1
       N66=N65+1
       N466=N4*N65+1
       N366=N466+N65
       N266=N366+N65
       N166=N266+N65
       NM66=N166+N65
       N6565=N65+N65
       AN2=-N2
       AN1=-N0
       NDIM=N
       NEQ=N64
  800  CONTINUE
C
C           MULTIPLICATION DES COEFFICIENTE DES FONCTION A TRAITER 
C           PAR LA MATRICE DE L'OPERATEUR EN QUESTION (MATRICE 
C           A 2 DIAGONALES APRES UN'ASTUCIEUSE COMBINAISON LINEAIRE
C           DES LIGNES)
C
       JM1=N66
       JM2=N65+N64
       DO 1 L=1,N3
       AL=-L
       DO 2 M=JM1,JM2
       CS(M)=(CC(M)+CC(M+N6565))*AL
  2    CONTINUE 
       JM1=JM1+N65
       JM2=JM1+N63
  1    CONTINUE 
C
C           INVERSION DE L'ASTUCIEUSE COMBINAISON LINEAIRE
C
       DO 4 M=N266,N266+N63
       Y(M)=CC(M)*AN2
  4    CONTINUE 
       DO 5 M=N166,N166+N63
       Y(M)=CC(M)*AN1
  5    CONTINUE 
C
       JM1=N366
       JM2=JM1+N63
       DO 6 L=2,N-2,2
       JM=N0-L
       C1=JM/DFLOAT(JM+2)
C
       DO 7 M=JM1,JM2
       CC(M)=C1*Y(M+N6565)+CS(M)
  7    CONTINUE 
       DO 8 M=JM1,JM2
       Y(M)=CC(M)
  8    CONTINUE 
       JM1=JM1-N6565
       JM2=JM1+N63
  6    CONTINUE 
C
       JM1=N466
       JM2=JM1+N63
       DO 9 L=3,N-3,2
       JM=N0-L
       C1=JM/DFLOAT(JM+2)
       DO 10 M=JM1,JM2
       CC(M)=Y(M+N6565)*C1+CS(M)
  10   CONTINUE 
       DO 11 M=JM1,JM2
       Y(M)=CC(M)
  11   CONTINUE 
       JM1=JM1-N6565
       JM2=JM1+N63
   9   CONTINUE 
C
       DO 12 M=1,N64
       Y(M)=0
  12   CONTINUE 
       DO 13 M=NM66,NM66+N63
       Y(M)=0
  13   CONTINUE 
       RETURN
       END
