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

       SUBROUTINE TFINS(N,CC,CS,Y)

		implicit double precision (a-h,o-z)

C
C **** ROUTINE TRANSFORMATION INVERSE DE FOURIER.
C ***** LES DIMENSIONS MINIMUN POUR LES TABLEAUX SONT N+1 POUR CC ET CS
C**** F DOIT AVOIR COMME DIMENSION N+3
C *** LES ARGUMENS DE LA ROUINE SONT: N NOMBRE DES POINTS D'ECHAN-
C *** TILLONAGE,(=2**m*3**l*5**j),
C ***  CC COMPOS.COSINUS,CS COMP.SINUS DE LA TRANSFORMEE DE FOURIER
C ***  Y  FONCTION OUTPUT.
C
C           Routine testee le 10/9/1975.
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:37:24  hyc
c *** empty log message ***
c
C Revision 1.1  1997/05/08 07:39:13  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

       DIMENSION TRIGS(1600),IFAX(40),CC(*),CS(*),AC(7,12),AS(7,12)
       DIMENSION Y(*)

	data	ndim/-1/

	save	N2,N1,N21,NDIM,AC,AS,ifax,trigs

        IF(N.LT.1033) GO TO 10
        PRINT 110,N
 110    FORMAT(10X,'N=',I4,'. N DOIT ETRE < 1033')
        CALL EXIT
  10    CONTINUE
       IF(NDIM.EQ.N) GO TO 20
       CALL FAX(IFAX,N,3)
       CALL FFTRIG(TRIGS,N,3)
       N2=N/2
       N1=N+1
       N21=N2+1
C
       NDIM=N
  300  FORMAT(10X,'LA SUB TF EST APPELLEE POUR LA PREMIERE FOIS')
C
C           PREPARATION DE LA MATRICE D'INVERSION DANS LE CAS
C           DE N<13.
C
       IF(N.GT.12) GO TO 20
       X0=0
       PIG=4*ACOS(X0)/N
       DO 31 L=1,N
       PIGL=PIG*(L-1)
       DO 30 J=2,N21
       J1=J-1
       AC(J,L)=COS(PIGL*J1)
       AS(J,L)=SIN(PIGL*J1)
30     CONTINUE
31     CONTINUE
       DO 32 L=1,N
       AC(1,L)=.5
       AC(N21,L)=AC(N21,L)*.5
32     CONTINUE
  20   CONTINUE
       IF(N.LT.13) GO TO 21
C
C           TF INVERSE POUR N>12.
C
       DO 33 L=1,N2
       L2=L+L
       L1=L2-1
       Y(L1)=CC(L)
       Y(L2)=-CS(L)
33     CONTINUE
       Y(N1)=CC(N21)
       CALL FFT991(Y,CC,TRIGS,IFAX,1,1,N,1,1)
       DO 34 L=1,N
       Y(L)=Y(L)*.5D+00
34     CONTINUE
 100     FORMAT(1X,'TF',10D12.4)
       RETURN
C
C           CALCUL DE LA TF INVERSE DANS LE CAS DE N<13.
C
  21   CONTINUE
  301  FORMAT(10X,'CALCUL DANS LE CAS N<13')
       DO 36 L=1,N
       SOM=CC(1)*AC(1,L)+CC(N21)*AC(N21,L)
       DO 35 J=2,N2
       SOM=SOM+CC(J)*AC(J,L)+CS(J)*AS(J,L)
35     CONTINUE
       Y(L)=SOM
36     CONTINUE
       RETURN
       END
