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

       SUBROUTINE CERARS(N,INDEX,F,CS,CC)

		implicit double precision (a-h,o-z)

C           
C           ROUTINE POUR LE CALCUL DES COEFFICIENTS DE
C           TCHEBYTCHEV D'UNE FONCTION  AVEC ECHANTILLONAGE
C           RAREFIE' A L'ORIGINE. L'ECHANTILLONAGE EST EFFECTUE'
C           DANS LES POINTS XL=-DCOS((L-1)*PI/(N*N)+PI/2) DE L'IN-
C           TERVALLE 0<X<1.
C      
C               ARGUMENTS DE LA ROUTINE:
C
C           N   = NOMBRE DES DEGRES DE LIBERTE'
C           INDEX  = PARAMETRE: SI LA FONCTION EST PAIRE INDEX
C                 DOIT ETERE=0, SI LA FONCTION EST IMPAIRE
C                INDEX=1
C           F   =FONCTION DONT ON VEUT CALCULER LES COEF-
C                DE TCHEBYTCHEV.
C           CS  = TABLEAU DE TRAVAIL
C           CC  = OUTPUT.(COEEF. DE TCHEBITCHEV DE F)
C......................................................................
C
C      Routine testee le 13/4/85.
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
c Revision 1.2  1997/05/23  11:33:32  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:40:40  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

       DIMENSION CC(*),CS(*), F(*),SEN(1025),COSE(1025)

	data	ncontr/0/

	save	NCONTR,X0,PI,PIT,N1,N2,N21,X5,SEN,COSE

       IF(N.EQ.NCONTR) GO TO 999
C
C           PREPARATION DES FONCTIONS AUXILIAIRES.
C
       NCONTR=N
       IF(N.LT.1025) GO TO 8
       PRINT 200,N
 200   FORMAT(10X,'DIMENSION DU TABLEAU SEN INSUFISANTES DANS LA'
     , ,' ROUTINE CHEIM: N=',I4)
       CALL EXIT
  8    CONTINUE
       X0=0
       PI=2*ACOS(X0)/N
       PIT=PI/2
       N1=N+1
       N2=N/2
       N21=N2+1
       X5=.5
       DO 7 L=1,N1
       SEN(L)=SIN((N2-L)*PI)+X5
       COSE(L)=SIN((L-1)*PIT)
  7    CONTINUE
 999   CONTINUE
C
       N21=N2+1
       IF(INDEX.EQ.0) GO TO 10
C
C           LA FONCTION IMPAIRE EST RENDUE PAIRE PAR MULTIPLICA-
C           TION AVEC SIN(TETA).
C
       DO 1100 L=1,N1
       F(L)=F(L)*COSE(L)
1100	CONTINUE
C
  10   CONTINUE
C           TRASFORMATION DE TCHEBYTCHEV DE LA FONCTION PAIRE.
C
       FN11=   (F(1)-F(N 1))*.5D+00
C
C ****   PONDERIZATION DE LA PARIE ASSYMETRIQUE  PAR LA FONCTION SINUS.
C
       DO 1 L=1,N2
       N21L=N21+L
       N20L=N21-L
       F1=F(N21 L  )
       F2=F(N20 L)
       F12=(F2-F1)*SEN(L)
       F(N20L)=F(N20L)-F12
       F(N21L)=F(N21L)+F12
1      CONTINUE
       CS(N21)=0
C
C ***   CALCUL DE COEFF. DE FOURIER DE LA FONTION PONDEREE.
       CALL TF9S(N,F,CS,CC)
C
C ***  REARENGEMENT DE COEFF DE FOURIER POUR OBTENIR LES COEFF. DE TCHEB
C
       C1=CS(1)
       C21=CS(N21)
       F(1)=0
       SOMM=0
       DO 2 L=2,N2
       F(L)=CC(L)+F(L-1)
       SOMM=SOMM+F(L)
2      CONTINUE
       F(1)=( FN11-SOMM)/N2
       F1=F(1)
       DO 4 L=2,N2
       CC(L)=CS(L)
4      CONTINUE
C
C
C ***    FORMATION DU TABLEAU CS CONTENANT LES COEFF. DE TCHEBYCHEC, ET
C ***   DU PREMIER COEFF.
C
       CC(1)=C1
       CS(N1)=C21
       DO 5 L=1,N2
       LL=L+L
       CS(LL)=-(F(L)+F1)
       CS(LL-1)=CC(L)
5      CONTINUE
       CS(2)=-F1
       IF(INDEX.EQ.1) GO TO 11
       DO 12 L=1,N1
       CC(L)=CS(L)
  12   CONTINUE
       RETURN
  11   CONTINUE
C
C           A PARTIR DE LA TRANSFORMME DE LA FONCTION PAIRE ON
C           OBTIENT LES COEFF. DE LA FONCTION IMPAIRE.
C
       CC(1)=-CS(1)
       DO 9 L=2,N
       CC(L)=-CC(L-1)-CS(L)-CS(L)
  9    CONTINUE
       CC(N1)=0  
 100   FORMAT(10X,10E12.4)
101    FORMAT(1X,'CHE64')
       RETURN
       END
C
