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

       SUBROUTINE TFIZMS(N,N64,CC,Y)

		implicit double precision (a-h,o-z)

C
C           ROUTINE POUR LA TF INVERSE MULTIPLE.
C
C           ARGUMENTS DE LA SUROUTINE:
C
C           N   =NOMBRE DES DGRES DE LIBERTEE. (N NOMBRE PAIR=
C                A 2**p*3**q*5**r p,q,r NOMBRES ENTIERS.
C
C           CC  =COEFFICIENTS COSINUS ET SINUS DE N64 FONCTIONS.
C           Y   = OUTPUT DES FONCTIONS DANS L'ESPACE DES x.
C                 LES VALEURES DES CES FONCTIONS SONT STOCKEES
C                 COMME DANS LA SUB TF.
C           CETTE ROUTINE EST SPECIALISEE ET NE DOIT ETRE EMPLOYEE
C           QUE AVEC LES ROUTINES CHIXMS,FUCI2S,FUCI3S...
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:36:19  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:22:11  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

       DIMENSION Y(*),CC(*),IFAX(64),TRIGS(1600)
       DATA NDIM/0/
       DATA NFON/0/
       DATA JJJ/0/

	save	NFON,NDIM,N65,NM65,ifax,trigs

       N1040=1040
       IF(N.GT.N1040) THEN
       PRINT 800,N
  800  FORMAT(10X,'DIMENSION INSUFFISANTES DANS LA ROUTINE'
     , ,' TFM, N=',I5,' > A LA DIMENSION DE TRIGS*2/3=',I5)
             ENDIF
C
       IF(N.EQ.NDIM) GO TO 1
       CALL FAX(IFAX,N,3)
       CALL FFTRIG(TRIGS,N,3)
C
  1    CONTINUE
       IF(NFON.EQ.N64.AND.NDIM.EQ.N) GO TO 4
C
C           SI N64 EST UN MULTIPLE DE 8 LE TEMPS CALCUL DU CRAY
C           ETRE MULTIPLIE PAR UN  FACTEUR 4. POUR EVITER
C           CELA ON SOCKE LES DONNEES TOUS LES N64+1 INTERVALLES
C           AU LIEU DE N64.(DANS LE CAS DE L INPUT EN PARALLEL)
C
       NFON=N64
       NDIM=N
       N65=N64
       IF((N64/8)*8.EQ.N64) N65=N64+1
C
       NM65=N*N65
C
  4    CONTINUE
C
C           TF INVERSE EN PARALLEL
C
       CALL FFT991(CC,Y,TRIGS,IFAX,N65,1,N,N64,1)
C
       DO 12 L=1,NM65
       Y(L)=CC(L)*.5
12     CONTINUE
  101  FORMAT(1X,'TFM')
  100  FORMAT(1X,'TFM',10D12.4)
       RETURN
C
       END
