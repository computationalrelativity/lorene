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
       SUBROUTINE CD1MRS(N,N64,IND,CC,CS,INDEX,F)

		implicit double precision (a-h,o-z)

C
C***** CALCUL DE LA DERIVEE PREMIERE DE N64 FONCTION A ECHANTILLONAGE 
C      RAREFIE A L'ORIGINE.
C           
C           ARGUMENTS DE LA SSUBROUTINE:
C           
C           N   = NOMBRE DES DEGRES DE LIBERTEE.
C           IND = PARMETRE: SI IND=0 ON TRAITE LES FONCTIONS
C                 PAIRES, SI IND=1 LES FONCTION IMPAIRES, SI
C                 IND=2 ON CALCULE LES COEFFICIENTS DE TCHE-
C                 BITCHEF DE LA FONCTION 1/r*D/dr DES FONCTIONS 
C                 PAIRES.
C           N64 = NOMBRE DES FONCTION DONT ON VEUT CALCULER LA 
C                 TRANSFORMEE.
C           CC  = TABLEAU D'ENTREE: IL CONTIENT LES COEFF. DE
C                 TCHEBYTCHEV DE LA FONCTION D'ONT ON VEUT CAL-
C                 CULER LA DERIVEE.
C           CS  = TABLEAU DE TRAVAIL ET TABLEAU OUTPUT SI
C                 INDEX=1.
C           INDEX = PARMETRE: SI INDEX=1 ON A EN SORTIE DANS CS
C                 LES COEFFICIENTS DE LA DERIVEE DE LA FONCTION
C               SI INDEX=2 ON EN SORTIE DANS Y LA FONCTION DE-
C                RIVEE.
C           Y   = TABLEAU OUTPUT
C
C           LES DIMENSIONS DES TABLEAUX DOIVENT ETRE GT.EQ. (N64+1)*(N+2
C  )
C****    SUBROUTINE DISTRUPTIVE.
C
C           SUBROUTINE TESTEE LE 08/9/1986.
C
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:04:12  hyc
c *** empty log message ***
c
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

       DIMENSION CC(*),CS(*),F(*)
C
       DATA NFONC/0/
       DATA NDIM/0/

	save	ndim,nfonc,n1,n11,n65,NM65,N6564,N6565
	save	N651,NM651,N63,NM64,NM66

       IF(NFONC.EQ.N64.AND.NDIM.EQ.N) GO TO 999
C
       NDIM=N
       NFONC=N64
       N1=N+1
       N11=N+2
       N65=N64
       IF((N64/8)*8.EQ.N64) N65=N64+1
       NM65=N65*N
       N6564=NM65+N64
       N6565=N65+N65
       N651=NM65+1
       NM651=NM65+N65
       N63=N64-1
       NM64=NM65-N65
       NM66=NM65+N65
  999  CONTINUE
C
C           CALCUL DE 1/R*D/DR D'UNE FONCTION PAIRE.
C
       IF(IND.EQ.2) THEN
C
C
       JM2=N6564
       DO 1 M=N651,JM2
       CC(M)=CC(M)*.5
  1    CONTINUE
C
       DO 2 L=1,N65
       CS(L)=0
   2   CONTINUE
C
       DO 3 M=1,N64
       F(M)=0
   3   CONTINUE
C
       JM1=NM65
       JM2=JM1-N65
       DO 4 LP=1,N,2
       L=N11-LP
       L1=L-1
       FAC=L1+L1
C
       DO 5 M=1,N64
       F(M)=F(M)-FAC*CC(M+JM1)
       CS(M+JM2)=-F(M)*4
  5    CONTINUE
C
       JM1=JM1-N6565
       JM2=JM1-N65
 4     CONTINUE
C
       DO 6 M=1,N64
       F(M)=0
   6   CONTINUE
C
       JM1=NM64
       JM2=JM1-N65
       DO 7 LP=2,N,2
       L=N11-LP
       L1=L-1
       FAC=L1+L1
       DO 8 M=1,N64
       F(M)=F(M)-FAC*CC(M+JM1)
       CS(M+JM2)=-4*F(M)
   8   CONTINUE
       JM1=JM1-N6565
       JM2=JM1-N65
   7   CONTINUE
       DO 9 M=N651,NM651
       CC(M)=CC(M)*2
   9   CONTINUE
C
       JM1=N651
       JM2=NM66
       DO 20 L=JM1,JM2
       CS(L)=0
  20   CONTINUE
C
       IF(INDEX.EQ.1) RETURN
       CALL CIRAMS(N,0,N64,CS,CC,F)
       RETURN
       ENDIF
C
C           CALCUL DE LA DERIVEE D'UNE FONCTION PAIRE.
C      
C
       IF(IND.EQ.0) THEN
C
       DO 10 M=N651,NM651
       CC(M)=CC(M)*.5
   10  CONTINUE
C
       DO 11 M=1,N64
       F(M)=0
   11  CONTINUE
C
       JM1=NM65
       JM2=JM1-N65
       DO 12 LP=1,N
       L=N1-LP
       FAC=4*L
       DO 13 M=1,N64
       F(M)=F(M)-CC(M+JM1)*FAC
       CS(M+JM2)=F(M)
  13   CONTINUE
       JM1=JM1-N65
       JM2=JM1-N65
  12   CONTINUE
C
       DO 14 M=N651,NM651
       CS(M)=0
  14   CONTINUE
C
       DO 15 M=N651,NM651
       CC(M)=CC(M)*2
  15   CONTINUE
C
       IF(INDEX.EQ.1) RETURN
       CALL CIRAMS(N,1,N64,CS,CC,F)
       RETURN
       ENDIF
C
C           CALCUL DE LA DERIVEE D'UNE FONCTION IMPAIRE.
C
C
       IF(IND.EQ.1) THEN
C
       DO 16 M=1,N64
       F(M)=0
  16   CONTINUE
C
       JM1=NM64
       DO 17 LP=1,N
       L=N1-LP
       FAC=4*L-2
       DO 18 M=1,N64
       F(M)=F(M)-CC(M+JM1)*FAC
       CS(M+JM1)=F(M)
  18   CONTINUE
       JM1=JM1-N65
  17   CONTINUE
C
       DO 19 M=N651,NM651
       CS(M)=0
  19   CONTINUE
       IF(INDEX.EQ.1) RETURN
       CALL CIRAMS(N,0,N64,CS,CC,F)
       ENDIF
       RETURN
  100  FORMAT(1X,10E12.4)
  101  FORMAT(1X,' ')
       END
