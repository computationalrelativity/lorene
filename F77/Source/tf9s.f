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

       SUBROUTINE TF9S(N,Y,CC,CS)

		implicit double precision (a-h,o-z)

C
C **** ROUTINE TRANSFORMATION DE FOURIER.
C ***** LES DIMENSIONS MINIMUN POUR LES TABLEAUX SONT N+1 POUR CC ET CS
C**** F DOIT AVOIR COMME DIMENSION N+3
C *** LES ARGUMENS DE LA ROUINE SONT: N NOMBRE DES POINTS D'ECHAN-
C *** TILLONAGE,(=2**m*3**l*5**j),
C ***  Y TABLEAU CONTENANT LA FONCTION DONT ON VEUT CALCULER LA T.F.
C ***  CC COMPOS.COSINUS,CS COMP.SINUS DE LA TRANSFORMEE DE FOURIER
C
C           Routine testee le 10/9/1985
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:38:04  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:42:33  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	DIMENSION TRIGS(1600),IFAX(40),CC(*),CS(*),Y(*)

	DATA NDIM/-1/
	save	ndim,n2,n1,n21,ifax,trigs

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
C
C
  20   CONTINUE
       DO 30 L=1,N1
       Y(L)=Y(L)+Y(L)
30     CONTINUE
       CALL FFT991(Y,CC,TRIGS,IFAX,1,1,N,1,-1)
       DO 31 L=1,N21
       CC(L)=Y(L+L-1)
       CS(L)=-Y(L+L)
31     CONTINUE
 100     FORMAT(1X,'TF',10D12.4)
       RETURN
       END
