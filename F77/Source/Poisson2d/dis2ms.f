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
       SUBROUTINE DIS2MS(N,N64,CC,CS,Y)

		implicit double precision (a-h,o-z)

C
C           ROUTINE EFFECTUANT LA DIVISION PAR SIN(TETA) 
C           D'UNE FONCTION DEVELOPEE EN SERIE DE POLYNOMES DE
C            TCHEBITCHEV DU 2ME GENRE.
C
C           DIVISION PAR SIN(TETA) D' UNE FONCTION DEVELOPEE EN DERIE
C           DE POLYNOMES DE TCHEBYTCHEV DU 2ME TYPE.
C
C           ARGUMENTS DE LAROUTINE:
C           N   =NOMPRE DE DEGRES DE LIBERTE-1
C           N64 =NOMBRE DE FONCTION A CRAYTINISER.
C           CC  =TABLEAU CONTENANT LES (N+1)*N64 COEFFICIENTS
C                    DU DEVELLOPEMENT EN POLYNOMES DE TCHEBYTCHEV
C                DU 2ME TYPE DES FONCTIONS A DIVISER PAR SIN(TETA).
C           CS  = TABLEAU DE TRAVAIL
C           Y   =OUTPUT DES FONCTIONS DIVISEES.
C                LE STOCKAGE EST EN PARALLEL (VOIR PAR EXEMPLE LA
C                ROUTINE TFNMS. DIMENSIONS MINIMES DES TABLEAUX
C                =(N+1)*(N64+1)
C
C           ROUTINE COMPLETEMENT CRAYTINISEE ET AYANT TESTEE AVEC
C           LE PROTOCOL ABITUEL LE 8/2/87.
C
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
c Revision 1.1  1997/10/23  08:41:03  eric
c Initial revision
c
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

       DIMENSION CC(*),Y(*),CS(*)
       DATA NDIM,NEQ/0,0/

	save	N1,N0,N2,N4,N65,N63,N653,N654,N655,N465,N4655
	save	N365,N3655,N265,N2655,N065,N0656,NM65,N6565,NM66
	save	NEQ,NDIM

C
       IF(N.EQ.NDIM) GO TO 80
       N1=N+1
       N0=N-1
       N2=N-2
       N4=N-4
  80   CONTINUE
C
       IF(NEQ.EQ.N64.AND.NDIM.EQ.N) GO TO 81
       N65=N64
       IF((N64/8)*8.EQ.N64) N65=N64+1
       N63=N64-1
       N653=N65*3
       N654=N653+N65
       N655=N654+N65
       N465=N4*N65
       N4655=N465+N64
       N365=N465+N65
       N3655=N365+N64
       N265=N365+N65
       N2655=N265+N64
       N065=N265+N65
       N0656=N065+1
       NM65=N065+N65
       N6565=N65+N65
       NM66=NM65+N65
       NEQ=N64
       NDIM=N
  81   CONTINUE
C
C           CALCUL DE LA FONCTION A DIVISER POUR TETA=PI/2
C
       DO 1 M=1,N64
       CS(M)=CC(M+N65)
  1    CONTINUE  
C
       JM1=N655
       JM2=JM1+N63
       DO 2 L=6,N,4
       DO 3 M=1,N64
       CS(M)=CS(M)+CC(M+JM1)
  3    CONTINUE  
       JM1=JM1+N654
       JM2=JM1+N63
  2    CONTINUE  
C
       JM1=N653
       JM2=JM1+N63
       DO 4 L=4,N,4
       DO 5 M=1,N64
       CS(M)=CS(M)-CC(M+JM1)
  5    CONTINUE  
       JM1=JM1+N654
       JM2=JM1+N63
  4    CONTINUE  
C
C           INVERSION DE LA MATRICE DIVISION PAR SIN(TETA),(L 1ER COEFF.
C           N'EST PAS DETREMINE')
C
       DO 6 M=N365+1,N3655
       Y(M)=CC(M+N65)*2
  6    CONTINUE  
C
       DO 7 M=N265+1,N2655
       Y(M)=CC(M+N65)*2
  7    CONTINUE  
C      
       JM1=N465+1
       JM2=N4655
       DO 8 L=1,N4
       DO 9 M=JM1,JM2
       CS(M)=Y(M+N6565)+CC(M+N65)*2
  9    CONTINUE  
C
       DO 10 M=JM1,JM2
       Y(M)=CS(M) 
  10   CONTINUE  
       JM1=JM1-N65
       JM2=JM1+N63
  8    CONTINUE  
C
       DO 11 M=N0656,NM66
       Y(M)=0
  11   CONTINUE  
C
C      CAL CALCUL DU 1ER COEFFICIENT. LA DETERMINATION DU 1ER COEFFICIEN
C  NT
C      EST EEFECTUEE EN IMPOSANT QUE LA FONCTION QUOTIENT SOIT IDENTIQUE
C  ,
C      POUR TETA=PI/2, A LA FONCTION A DIVISER.
C
C           DETERMINATION DE LA FONCTION A DIVISER PUR TETA=PI/2
C
       JM1=N65+N65
       JM2=JM1+N63
       DO 12 L=3,N,4
       DO 13 M=1,N64
       CS(M)=CS(M)+Y(M+JM1)
  13   CONTINUE  
       JM1=JM1+N654
       JM2=JM1+N63
  12   CONTINUE  
C
       JM1=N654
       JM2=JM1+N63
       DO 14 L=5,N,4
       DO 15 M=1,N64
       CS(M)=CS(M)-Y(M+JM1)
  15   CONTINUE  
       JM1=JM1+N654
       JM2=JM1+N63
  14   CONTINUE  
C
       DO 16 M=1,N64
       Y(M)=CS(M)*2
  16   CONTINUE  
  100  FORMAT(1X,10E12.4)
  101  FORMAT(1X,' ')
       RETURN
       END     
C
