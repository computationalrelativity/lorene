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

       SUBROUTINE DERFMS(N,N64,CC,CS,Y)

		implicit double precision (a-h,o-z)

       DIMENSION CC(*),CS(*),Y(*)
C
C           ROUTINE POUR LE CALCUL DE L'OPERATEUR
C
C           1/X*D/DX-1/X**2
C
C           APPLIQUE' AUX FONCTIONS IMPAIRES AVEC ECHANTILLONAGE RAREFIE
C   A
C           L'ORIGINE.
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
C           CS  =TABLEAU DE TRAVAIL
C           Y   =TABLEAU OUTPUT
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
c Revision 1.2  1997/05/23  11:30:23  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:34:38  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

       DATA NDIM,NEQ/0,0/

	save	N0,N2,N3,N4,N65,N63,N66,N466,N366,N266
	save	N166,NM66,N6565,N1M65,NDIM,NEQ

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
       N1M65=(N+1)*N65
       NDIM=N
       NEQ=N64
  800  CONTINUE
C
C           MULTIPLICATION DES COEFFICIENTE DES FONCTION A TRAITER 
C           PAR LA MATRICE DE L'OPERATEUR EN QUESTION (MATRICE 
C           A 1 DIAGONALES APRES UN'ASTUCIEUSE COMBINAISON LINEAIRE
C           DES LIGNES)
C
       JM1=1
       JM2=JM1+N63
       AL=8
       DO 1 L=1,N0
       DO 2 M=JM1,JM2
       CS(M)=CC(M+N65)*AL
  2    CONTINUE
       JM1=JM1+N65
       JM2=JM1+N63
       AL=AL+8
  1    CONTINUE
C
C           INVERSION DE L'ASTUCIEUSE COMBINAISON LINEAIRE
C
       DO 3 M=N166,N1M65
       Y(M)=0
   3   CONTINUE
C               
       JM1=N266
       JM2=JM1+N63
       B2=N-2
       DO 4 M=JM1,JM2
       Y(M)=CS(M)
   4   CONTINUE
C
       JM1=JM1-N65
       JM2=JM1+N63
       B2=N-1
       DO 5 L=1,N-2
       C2=1/B2
       C3=1-C2
       DO 6 M=JM1,JM2
       CC(M)=CS(M)+Y(M+N65)*C2+Y(M+N6565)*C3
  6    CONTINUE
       DO 7 M=JM1,JM2
       Y(M)=CC(M)
   7   CONTINUE
       JM1=JM1-N65
       JM2=JM1+N63
       B2=B2-1
  5    CONTINUE
  101  FORMAT(1X,' ')
  100  FORMAT(1X,10E12.4)
       RETURN
       END
