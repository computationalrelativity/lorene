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

       SUBROUTINE DISEMS(N,N64,ITCH,IND,CC,CS,Y)

		implicit double precision (a-h,o-z)

       DIMENSION CC(*),CS(*),Y(*)
C
C           ROUTINE POUR LA MULTIPLICATION (DIVISION) PAR SIN(TETA)**2
C           DES FONCTIONS DEVELOPEES EN SERIE DE TCHEBYTCHEV DU 1ER OU
C           2ME TYPE.
C
C      ARGUMENTS DE LA ROUTINE:
C
C           N   =NOMBRE DE DEGRES DE LIBERTE'-1
C           N64 =PARAMETRE DEFINISSANT LE NOMBRE DES FONCTIONS QUI
C                DOIVENT ETRE TRAITEES SIMULTANEMENT
C           ITCH =PARAMETRE: SI ITCH=1 ON MULTIPLIE LA FONCTION PAR
C                SIN(TETA)**2,SI =2 ON EFFECTUE LA DIVISION PAR
C                SIN(TETA)**2
C           IND =PARAAMETRE SI IND=1,ON EFFECTUE LA DIVSION (OU MULTI-
C                PLICATION) PAR SIN(TETA)**2 DES FONCTIONS DEVELOPPEES
C                EN POLYNOMES DE TCHEBYTCHEV DU 1ER TYPE, SI IND=2, ON 
C                EFFECTUE LES MEMES OPERATIONS SUR LES POLYNOMES DU 
C                2ME TYPE.
C           CC  =TABLEAU IMPUT CONTENANT LES COEFFICIENTS DE TCHBYTCHEV
C               DES N64 FONCTION A TRAITER STOCKES EN PARALEL,
C               (CC(1),CC(2),...CC(N64) SONT LE PREMIERS COEFFICIENTS
C               DES N64 FONCTIONS, CC(N64+1),CC(N64+2),...CC(N64+N64)
C               LES 2EMES COEFFICIENTS DE N64 FONCTIONS ET AINSI DE 
C               SUITE.
C           CS  =TABLEAU DE TRAVAIL
C           Y   =OUTPUT: LES RESULTATS DES OPERATIONS DE DIVISION
C                OU MULTIPLICATION SONT STOCKES EN PARALLEL DANS Y.
C
C           Routine teste' le 3 mars 1987
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:29:03  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:34:47  hyc
C Initial revision
C
C
C $Header$
C
C
	character* 100 header
	data header/'$Header$'/

       DATA NDIM,NEQ/0,0/

	save	N1,N0,N2,N3,N4,N5,N65,N63,N66,N466,N366,N266,N166
	save	NM66,N6565,N6566,N1M65,N1301,A75,A25,ANM,ANM1,AB1
	save	NDIM,NEQ

C
C           INITIALISATION
C
       IF(NDIM.NE.N) THEN
       N1=N+1
       N0=N-1
       N2=N-2
       N3=N-3
       N4=N-4
       N5=N-5
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
       N6566=N6565+1
       N1M65=N1*N65
       N1301=N6565+N6566
       A75=.75
       A25=-.25
       ANM=2./N
       ANM1=N*.125
       AB1=N*.25-.25 
       NDIM=N
       NEQ=N64
  800  CONTINUE
C
C11111111111111111111111111111111111111111111111111111111111111111111111
C  111111
C
C           OPERATIONS SUR LES FONCTIONS DEVELOPEES EN POLYNOMES
C           DE TCHEBYTCHEV DU 1ER TYPE.
C
C
C           MULTIPLICATION PAR SIN(TETA)**2 DES POLYNOMES DU 1ER TYPE
C
       IF(IND.EQ.1) THEN
       IF(ITCH.EQ.1)THEN
C
C           PARTIE PAIRE
C
       DO 1 M=1,N64
       Y(M)=(CC(M)-CC(M+N6565))*.5
   1   CONTINUE
C
       DO 2 M=N6566,N6566+N63
       Y(M)=CC(M)*.5-(CC(M+N6565)+CC(M-N6565))*.25
   2   CONTINUE
C
       JM1=N1301
       JM2=JM1+N63
       DO 3 L=5,N0,2
       DO 4 M=JM1,JM2
       Y(M)=CC(M)*.5-(CC(M-N6565)+CC(M+N6565))*.25
   4   CONTINUE
       JM1=JM1+N6565
       JM2=JM1+N63
   3   CONTINUE
C
       DO 5 M=NM66,NM66+N63
       Y(M)=(CC(M)-CC(M-N6565))*.5
   5   CONTINUE
C
C           MULTIPLICATION PARTIE IMPAIRE:
C
       DO 6 M=N66,N66+N63
       Y(M)=(CC(M)-CC(M+N6565))*.25
   6   CONTINUE
C
       JM1=N66+N6565
       JM2=JM1+N63
       DO 7 L=4,N2,2
       DO 8 M=JM1,JM2
       Y(M)=CC(M)*.5-(CC(M-N6565)+CC(M+N6565))*.25
   8   CONTINUE
       JM1=JM1+N6565
       JM2=JM1+N63
   7   CONTINUE
C
       DO 9 M=N166,N166+N63
       Y(M)=(CC(M)-CC(M-N6565))*.25
   9   CONTINUE
       RETURN
       ENDIF
C
       IF(ITCH.EQ.2) THEN
C      
C           DIVISION PAR SIN(TETA)**2 DES POLYNOMES DU PREMIER TYPE
C
C           COMBINAISON LINEAIRE DES COEFFICIENTS DU 2ME TERME DU SYSTEM
C  E.
C           PARTIE PAIRE.
C
       JM1=N6566
       JM2=JM1+N63
       DO 10 M=JM1,JM2
       CS(M)=CC(M)+CC(M-N6565)*.5
  10   CONTINUE
       DO 11 M=JM1,JM2
       CC(M)=CS(M)
  11   CONTINUE
C
       JM1=JM1+N6565
       JM2=JM1+N63
       DO 12 L=5,N0,2
       DO 13 M=JM1,JM2
       CS(M)=CC(M)+CC(M-N6565)
  13   CONTINUE
       DO 14 M=JM1,JM2
       CC(M)=CS(M)
  14   CONTINUE
       JM1=JM1+N6565
       JM2=JM1+N63
  12   CONTINUE
C
       DO 15 M=JM1,JM2
       CS(M)=CC(M)+2*CC(M-N6565)
  15   CONTINUE
       DO 16 M=JM1,JM2
       CC(M)=CS(M)
  16   CONTINUE
C
C           INVERSION MATRICE POLYNOMES DU 1ER TYPE (PARTIE PAIRE)
C
       JM1=N266
       JM2=JM1+N63
       DO 17 M=JM1,JM2
       Y(M)=CC(M)*4
  17   CONTINUE
C
       JM1=JM1-N6565
       JM2=JM1+N63
       DO 18 L=1,N5,2
       DO 20 M=JM1,JM2
       CS(M)=CC(M)*4+Y(M+N6565)
  20   CONTINUE
       DO 19 M=JM1,JM2
       Y(M)=CS(M)
  19   CONTINUE
       JM1=JM1-N6565
       JM2=JM1+N63
  18   CONTINUE
C
       DO 21 M=1,N64
       CS(M)=CC(M)+.5*Y(M+N6565)
  21   CONTINUE
       DO 22 M=1,N64
       Y(M)=CS(M)*2.
  22   CONTINUE
C
C           COMBINAISON LINEAIRE COOEFFICIENT DU 2ME MEMBRE: PARTIE PAIR
C  E
C           (POLYNOMES DU 1ER TYPE)
C
       JM1=N66+N6565
       JM2=JM1+N63
       DO 23 L=2,N2,2
       DO 24 M=JM1,JM2
       CS(M)=CC(M)+CC(M-N6565)
  24   CONTINUE
       DO 25 M=JM1,JM2
       CC(M)=CS(M)
  25   CONTINUE
       JM1=JM1+N6565
       JM2=JM1+N63
  23   CONTINUE
C
C           INVERSION MATRICE PARTIE IMPAIRE (POLYNOMES 1ER TYPE)
C
       JM1=N366
       JM2=JM1+N63
       DO 26 M=JM1,JM2
       Y(M)=CC(M)*4
  26   CONTINUE
       JM1=JM1-N6565
       JM2=JM1+N63
C
       DO 27 L=2,N4,2
       DO 28 M=JM1,JM2
       CS(M)=CC(M)*4+Y(M+N6565)
  28   CONTINUE
       DO 29 M=JM1,JM2
       Y(M)=CS(M)
  29   CONTINUE
       JM1=JM1-N6565
       JM2=JM1+N63
  27   CONTINUE
C
       DO 30 M=N166,N1M65
       Y(M)=0
  30   CONTINUE
C
       RETURN
       ENDIF  
             ENDIF 
C
C22222222222222222222222222222222222222222222222222222222222222222222222
C  222222 
C
C           MULTIPLICATION PAR SIN(TETA)**2 DES POLYNOMES DU 2ME TYPE
C
       IF(IND.EQ.2) THEN
       IF(ITCH.EQ.1) THEN
C
C           MULTIPLICATION DES COEFFICIENTE DES FONCTION A TRAITER 
C           PAR LA MATRICE DE L'OPERATEUR EN QUESTION (MATRICE 
C           A 2 DIAGONALES APRES UN'ASTUCIEUSE COMBINAISON LINEAIRE
C           DES LIGNES)
C
C
       JM1=N66
       JM2=N65+N64
       DO 41 M=JM1,JM2
       Y(M)=CC(M)*A75+CC(M+N6565)*A25
  41   CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
       DO 42 M=JM1,JM2
       Y(M)=CC(M)*.5+CC(M+N6565)*A25
  42   CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
       DO 43 L=4,N2
       DO 44 M=JM1,JM2
       Y(M)=(CC(M+N6565)+CC(M-N6565))*A25+CC(M)*.5 
  44   CONTINUE
       JM1=JM1+N65
       JM2=JM1+N63
  43   CONTINUE
C
       DO 45 M=N266,N266+N63
       Y(M)=CC(M-N6565)*A25+CC(M)*.5
  45   CONTINUE
C
       DO 46 M=N166,N166+N63
       Y(M)=CC(M-N6565)*A25+CC(M)*A75
  46   CONTINUE
C
       DO 40 M=1,N64
       Y(M)=0
  40   CONTINUE
C
       DO 62 M=NM66,N1M65
       Y(M)=0
  62   CONTINUE
       RETURN
       ENDIF
C
       IF(ITCH.EQ.2) THEN
C
C           DIVISION PAR SIN(TETA)**2
C
C           SOLUTION DU SYSTEME: PARTIE PAIRE:(TCHEBITCHEV TYPE 2)
C
C           COMBINAISON LINEAIRE DES COEFFICIENTS DU 2ME MEMBRE
C           TERMES PAIRES.
C
       JM1=N1301
       JM2=JM1+N63
       DO 47 L=3,N2,2
       AL=(L+1)/2
       DO 48 M=JM1,JM2
       CS(M)=CC(M)*AL+CC(M-N6565)
  48   CONTINUE
       DO 49 M=JM1,JM2
       CC(M)=CS(M)
  49   CONTINUE
       JM1=JM1+N6565
       JM2=JM1+N63
  47   CONTINUE
C
C           INVERSION PARTIE PAIRE.
C
       B1=ANM1
       C1=1/B1
       DO 50 M=N266,N266+N63
       Y(M)=CC(M)*C1
  50   CONTINUE
C
       B1=ANM1-.25
       B2=B1-.25
       JM1=N466
       JM2=JM1+N63
       DO 51 L=1,N5,2
       C1=1./B1
       DO 52 M=JM1,JM2
       CS(M)=(Y(M+N6565)*B2+CC(M))*C1   
  52   CONTINUE
       DO 53 M=JM1,JM2
       Y(M)=CS(M)
  53   CONTINUE
       JM1=JM1-N6565
       JM2=JM1+N63
       B1=B2
       B2=B1-.25
  51   CONTINUE
C
C           SOLUTION PARTIE IMPAIRE (TCHEBYTCHEV TYPE 2)
C
C           COMBINAISON LINEAIRE DES COEFFICIENTS DU 2ME MEMBRE
C           TERMES IMPAIRES.
C           
       JM1=N66+N6565
       JM2=JM1+N63
       DO 54 L=2,N2,2
       AL=L+1
       DO 55 M=JM1,JM2
       CS(M)=CC(M)*AL+CC(M-N6565)
  55   CONTINUE
       DO 56 M=JM1,JM2
       CC(M)=CS(M)
  56   CONTINUE
       JM1=JM1+N6565
       JM2=JM1+N63
  54   CONTINUE
C
C           INVERSION PARTIE IMPAIRE
C
       DO 57 M=N166,N166+N63
       Y(M)=CC(M)*ANM
  57   CONTINUE
C
       B1=AB1
       B2=B1-.5 
       JM1=N366
       JM2=N366+N63
       DO 59 L=2,N2,2
       C1=1./B1
       DO 60 M=JM1,JM2
       CS(M)=(CC(M)+Y(M+N6565)*B2)*C1
  60   CONTINUE
       DO 61 M=JM1,JM2
       Y(M)=CS(M)
  61   CONTINUE
       JM1=JM1-N6565
       JM2=JM1+N63
       B1=B2
       B2=B2-.5
  59   CONTINUE
C
       DO 63 M=N266,N1M65
       Y(M)=0
  63   CONTINUE
C
       DO 64 M=1,N64
       Y(M)=0
  64   CONTINUE
C
 100   FORMAT(1X,10E12.4)
 101   FORMAT(1X,' ')
       RETURN
       ENDIF
       ENDIF
       END
