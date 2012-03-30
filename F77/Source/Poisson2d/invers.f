C
C   Copyright (c) 2000 Silvano Bonazzola
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
	SUBROUTINE INVERS(M,NM,A,EPS,E)
C
		implicit double precision (a-h,o-z)

	PARAMETER (N128=513)
C INVERSION D'UNE MATRICE CARRE PAR LA METHODE DES PIVOTS
C
C ARGUMENTS :
C
C M		ORDRE DE LA MATRICE
C NM		PREMIERE DIMENSION DU TABLEAU A
C A(NM,1)	TABLEAU CONTENANT LA MATRICE A INVERSER EN ENTREE
C		ET LA MATRICE INVERSE EN SORTIE
C EPS		ERREUR MAXIMALE TOLEREE (ENTRE)
C E		VALEUR DU DETERMINANT
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
c Revision 1.3  2000/11/28  15:08:22  nbib
c *** empty log message ***
c
c Revision 1.2  2000/11/28  15:03:22  nbib
c Correction probleme de FORMAT (label 1000, ligne 114) sous GNU g77.
c
c Revision 1.1  2000/10/23  14:40:59  nbib
c Initial revision
c
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

	INTEGER Q,S,R
	DIMENSION A(NM,*),Q(N128),S(N128)
C
C
	B(I,J)=A(I,J)
C
	IF(M.GT.N128) THEN
	PRINT*,'DIMENSIONS INSUFFISANTES DAN LA ROUTINE INVERS,M=',M
	CALL EXIT
	ENDIF
C
	N=M
C	INITIALIZE
	D=1.
	DO 1 I=1,N
1	Q(I)=I
C	SEARCH FOR OPTIMAL PIVOTAL ELEMENT
	DO 13 I=1,N
	R=I
	L=I
	C=ABS(B(Q(I),Q(I)))
	DO 5 J=I,N
	DO 5 K=I,N
	T=ABS(B(Q(J),Q(K)))
	IF(.NOT.(C.LT.T))GO TO 5
	R=J
	L=K
	C=T
5	CONTINUE
	L1=L-1
	R=Q(R)
	L=Q(L)
	IF(L1.LT.1) GOTO 70
	DO 7 J=L1,I,-1
7	Q(J+1)=Q(J)
70	Q(I)=R
	S(I)=L
C	INTERCHANGE ROWS IN ORDER TO PLACE THE PIVOTAL ELEMENT ON THE
C	MAIN DIAGONAL
	IF(R.EQ.L)GO TO 9
	D=-D
	DO 8 J=1,N
	T=A(R,J)
	A(R,J)=A(L,J)
8	A(L,J)=T
9	T=A(L,L)
	D=D*T
	IF(C.LT.EPS)GO TO 20
	T=1./T
	A(L,L)=1.
C	DIVIDE THE PIVOTAL ROW BY THE PIVOTAL ELEMENT
	DO 10 J=1,N
10	A(L,J)=A(L,J)*T
C	REDUCE NON PIVOTAL ROWS
	L2=L-1
	IF(L2.EQ.0) GOTO 121
	L1=1
11	DO 12 J=L1,L2
	T=A(J,L)
	A(J,L)=0.
	DO 12 K=1,N
12	A(J,K)=A(J,K)-A(L,K)*T
121	L1=L2+2
	L2=N
	IF(.NOT.(L1.GT.N))GO TO 11
13	CONTINUE
C	INTERCHANGE COLUMNS
	DO 16 J=N,1,-1
	R=Q(J)
	L=S(J)
	IF(R.EQ.L) GOTO 16
	DO 15 K=1,N
	T=A(K,R)
	A(K,R)=A(K,L)
15	A(K,L)=T
16	CONTINUE
18	E=D
	RETURN
20	CONTINUE     
	D=0.
	PRINT 1000,L,I,T
	GO TO 18
1000  	FORMAT(10X,'CICVERDE MATRICE SING. A L ADRESSE',I7,
     1   '   PIVOT NR', I3, ' =',G14.8)
	END
