C
C   Copyright (c) 1998 Silvano Bonazzola
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
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.1  1998/06/22  10:38:21  eric
c Initial revision
c
C
C $Header$
C
C

	SUBROUTINE DIVQ1S(NDL,NDIM,IND,RR1,CS,CC)

	IMPLICIT double PRECISION(A-H,O-Z)

C
C##	version du 19.11.1993: derniere dimension des tableaux avec *
C				modification des commentaires 
C
C		ROUTINE EFFECTUANT LA DIVISION OU LA MULTIPLICATION
C		PAR r DANS L'ESPACE DES COEFFICIENTS DANS LES PROBLEMES
C		MULTIZONES (COQUILLES SPHERIQUES), r= RR1(LZ)+x*RR2(LZ)
C		0 < LZ < NOMBRE DES ZONES +1, 0 < x < 1 POUR LZ=1,(RR1(1)=0)
C		0 < x < 2 POUR LZ > 1
C		
C		ARGUMENTS DE LA ROUTINE:
C		
C		NDL	= TABLEAU CONTENANT EN NDL(1) LE NOMBRE DES ZONES
C			  ET DANS NDL(2), NDL(3),... LE NOMBRE DES DGRES
C			  DE LIBERTE DE LA 1ERE ZONE, 2ME ZONE...
C		NDIM	= DIMENSIONS DES TABLEAUX COMME DECLARE DANS LE
C			  PROGRAMME APPELLANT.
C
C		IND	= PDRAPEAU
C
C		    IND=0 	ON EFFECTUE LA LA DIVISION
C				PAR r DES FONCTIONS PAIRES
C	            IND=1	ON EFFECTUE LA LA DIVISION
C				DES FONCTIONS IMPAIRES.
C	            IND=2 	ON EFFECTUE LA MULTIPLICATION PAR r
C			        DES FONCTION PAIRES.
C	            IND=3 	ON EFFECTUE LA MULTIPLICATION PAR r
C			        DES FONCTION IMPAIRES.
C
C		RR1	=TABLEAU CONTENANT LE RAYON INTERNE DES DIFFERENTES
C			 COQUILLES
C		CS	=TABLEAU D'ENTREE: IL CONTIENT LES COEFFICIENTS DE
C			 TCHEBYTCHEV DE LA FONCTION QUI DOIT ETRE DIVISEE
C			 MULTIPLIEE PAR r.
C		CC	=TABLEAU SORTIE: COEFFICIENTS DU RESULTAT
C
C
C		ROUTINE TESTEE LE 20/08/91.
C
C

	character*100 header
	data header/'$Header$'/

	PARAMETER (N257=257)
C
	DIMENSION CC(NDIM,*),CS(NDIM,*),RR1(*),TC(N257),TS(N257)
	DIMENSION NDL(*)
C
C
	NZON=NDL(1)
	DO 40 LZ=2,NZON+1
	IF(NDL(LZ).GT.N257) THEN
	PRINT*,'DIMENSIONS INSUFFISANTES DANS LA ROUTINE DIVQ1S'
	CALL EXIT
	ENDIF
   40	CONTINUE
C
	N1=NDL(2)
	N=N1-1
	CALL DIXR1S(N,IND,CS,CC)
	X2=RR1(2)	
	IF(IND.LT.2 )X2=1/RR1(2)
	DO 1 L=1,N1
	CC(L,1)=CC(L,1)*X2
   1	CONTINUE
C
	IF(NZON.EQ.1) RETURN
C
	ID=1
	IF(IND.EQ.2.OR.IND.EQ.3) ID=0
C
	DO 2 LZON=2,NZON
C
	N1=NDL(LZON+1)
	N=N1-1
C
	R1=RR1(LZON)
	R2=(RR1(LZON+1)-RR1(LZON))*.5
	DO 3 L=1,N1
	TS(L)=CS(L,LZON)
   3	CONTINUE
C
	CALL DIRCMS(N,N257,1,ID,R1,R2,TS,TC)
	DO 4 L=1,N1
	CC(L,LZON)=TC(L)
  4   	CONTINUE
  2   	CONTINUE
C
	RETURN
C  100	FORMAT(1X,10E10.3)
C  101	FORMAT(1X,' ')
	END
