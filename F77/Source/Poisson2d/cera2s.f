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

       SUBROUTINE CERA2S(N,INDE,F,CS,CC)

		implicit double precision (a-h,o-z)

C
C      ROUTINE POUR LE DEVELOPPEMENT D'UNE FONCTION EN POLYNOMES
C      DE TCHEBYTCHEV DU PREMIER ORDRE OU DU 2me ORDRE AYANT UNE 
C      SYMMETRIE PAR RAPPORT TETA=PI/2.
C
C           ARGUMENTS DE LA ROUTINE:
C
C      N    = DEGRES DE LIBERTE-1
C      INDE =PARAMETRE: SI INDE=0 LA FONCTION A DEVELOPPER
C            DOIT ETRE PAIRE PAR RAPPORT TETA=PI/2,
C            SI INDE=1 ELLE DOIT TRE IMPAIRE.
C
C      F    + TABLEAU CONTENANT L'ECHEANTILLONAGE DE LA FONCTION
C            QU'ON VEUT TRANSFORMER.
C      CS   = TABLEAU DE TRAVAIL 
C      CC   = TABLEAU CONTENANT EN SORTIE LE COEFFICIENTS DE LA
C             TRANSFORMATION
C           DIMENSIONS MINIMUM DES TABLEAUX NDIM > N+3
C
C           L'IMPUT EST DETRUIT
C
C
C
C           Routine testee le 15/3/1986
C
C
C $Id$
C $Log$
C Revision 1.2  2012/03/30 12:12:42  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.2  1997/05/23  11:33:43  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:40:44  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

       DIMENSION CC(*),CS(*), F(*),SEN(520),COSE(1050)
       DATA NDIM/0/

	save	NDIM,N1,N2,N21,N10,X0,PIG
	save	SEN,PI2,COSE

C
       N1024=1024
       IF(N.GT.N1024) THEN
       PRINT 400,N1024,N
  400  FORMAT(10X,'DAN LA SUB. CHE64 LES DIMENSIONS SONT'
     , ,' INSUFFISANTES, DIMENS. MAX.=',I5,' N=',I5)
       CALL EXIT
       ENDIF
C
       IF(NDIM.EQ.N) GO TO 4
       NDIM=N
       N1=N+1
       N2=N/2
       N21=N2+1
       N10=N-1
       X0=0
       PIG=2*ACOS(X0)/N
       DO 10 L=1,N2
       X=(L-1)*PIG
       SEN(N21-L)=(SIN(X)-.5)
10     CONTINUE
C
       PI2=PIG/2
       DO 11 L=1,N1
       COSE(L)=SIN((L-1)*PI2)
  11   CONTINUE
  4    CONTINUE
C
C
C ****   PONDERIZATION DE LA PARIE ASSYMETRIQUE  PAR LA FONCTION SINUS.
C
       F11=F(1)
       DO 12 L=1,N1
       IF(INDE.EQ.1) GO TO 13
       F(L)=F(L)*COSE(L)
  12   CONTINUE
  13   CONTINUE
C
       DO 1 L=1,N2
       N21L=N21+L
       N20L=N21-L
       F1=F(N21 L  )
       F2=F(N20 L)
       F12=(F2+F1)*SEN(L)
       F(N20L)=F(N20L)+F12
       F(N21L)=F(N21L)+F12
  1    CONTINUE
       F(N21)=F(N21)+F(N21)
C
C ***   CALCUL DE COEFF. DE FOURIER DE LA FONTION PONDEREE.
C
       CALL TF9S(N,F,CC,CS)
C
C           RERARENGEMENT DE COEFFICIENTS DE FOURIER POUR
C           OBTENIR LES COEFFICIENTS DE TCHEBYTCHEV
       
C
       F(N)=-CC(N21)*.5
       DO 2 L=1,N2-1
       LN=N21-L
       LN2=LN+LN
       F(LN2-2)=F(LN2)-CC(LN)
   2   CONTINUE
       F(2)=CC(1)*.5
C
       DO 3 L=1,N2
       F(L+L-1)=CS(L)
   3   CONTINUE
C
C
       DO 5 L=2,N,2
       CC(L)=-F(L)
  5    CONTINUE
       DO 6 L=1,N1,2
       CC(L)=F(L)
  6    CONTINUE
       CC(1)=0
       CC(N1)=0
       IF(INDE.EQ.1)RETURN
C
C           "MULTIPLICATION PAR SINUS" DE COEFFICIENTS.
C
C           ON FAIT L'HYPOTHESE QUE F(N)=0
C
       F(N)=0
       DO 1000 L=1,N-1
       LN=N1-L
       F(LN-1)=-2*CC(LN)-F(LN)
1000	CONTINUE
C
C           DETERMINATION DU Neme COEFFICIENT.
C
       SOMM=0
       DO 1100 L=2,N-1,2
       SOMM=SOMM-F(L)
1100	CONTINUE
C
       DO 1200 L=1,N-1,2
       SOMM=SOMM+F(L)
1200	CONTINUE
C
C           SOMM EST = A LA VALEUR DE LA FONCTION EN TETA=PI/2
C           SOUS L'HYPOTHESE QUE F(N)=0
C
       SOMM=(SOMM-F11)/N
C
C           ON AJOUTTE UN DELTA DE DIRAC CORRESPONDANTE
C           A LA DIFFERENCE ENTRE LA VALEUR DE LA FONCTION
C           EN PI/2 ET LA VALEUR TROUVEE AVEC L'HYPOTHESE CC(N)=0
C
       DO 1300 L=1,N,2
       CC(L)=F(L)-SOMM
1300	CONTINUE
       DO 1400 L=2,N,2
       CC(L)=F(L)+SOMM
1400	CONTINUE
C
100    FORMAT(1X,'CHE64',10D12.4)
101    FORMAT(1X,'CHE64')
109    FORMAT(1X,I5,10D12.4)
  110  FORMAT(1X,3I4,5E12.4)
       RETURN
       END
