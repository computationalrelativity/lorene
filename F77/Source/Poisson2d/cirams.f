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

       SUBROUTINE CIRAMS(N,INDEX,N64,F,CC,CS)

		implicit double precision (a-h,o-z)

C
C **** INVERSION DE LA TF DE TCHEBYTCHEV DE N64 FONCTIONS ECHANTILLONEES
C      RAREMENT A L'ORIGINE. (INTERVALLE D'ECHANTILLONAGE [0-1]
C **** DIMENSION MINIMALE DE TABLEUX=N+1
C **** ROUTINE DISTRUPTIVE9
C
C           ARGUMENTS DE LA ROUTINE:
C
C      N    = NOMBRE DES DEGRES DE LIBERTE'
C      INDEX = PARAMETRE: IL DOIT ETRE =0 SI LA FONCTION D'ENTREE
C             EST PAIRE, =1 SI IMPAIRE.
C      N64  =NOMBRE DES FONCTIONS A INVERSER.
C      F    = IMPUT: TABLEAU CONTENANT LES COEFFICIENTS DE TCHE-
C           BITCHEV DES FONCTIONS A INVERSER.
C            F DOIT CONTENIR LES COEFFICIENTS DE N64 FONCTIONS 'EN PARAL
C  -
C            LEL' (Cfr. LA ROUTINE TFNMS). SI LE NOMBRE DES FONCTIONS
C            A INVERSER EST UN MULTIPLE DE 8, N64 DOIT ETRE, POUR DES RA
C  I-
C            SONS DE CRAYTINISATION = LE NOMBRE DES FONCTIONS=1.
C
C      CC   = TABLEAU DE TRAVAIL.
C      CS   = OUTPUT.
C
C               LES DIMENSIONS DES TABLEAUX DOIVENT ETRE GT.EQ.
C               (N64+1)*(N+2)
C               
C           Routine testee le 8/10/1986.
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:32:03  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:40:35  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

       DIMENSION CC(*),F1(257),FN(257),CS(*),F(*),SEN(1025),COSE(1025)

       DATA NCONTR/0/,NFONC/0/

	save	M1025,N1,N11,N2,N21,N3,X0,PI,PI2,SEN,COSE
	save	N257,NCONTR,NFONC,N65,N63,NM65,NM651,NM66,N2651
	save	N265,N6565,N633,N365

       IF(NCONTR.EQ.N) GO TO 999
       M1025=257
       IF(N.LT.M1025) GO TO 9
       PRINT 200,N
  200  FORMAT(10X,'DIMENSIONS INSUFFISANTES DAN CIRAMS:N=',I5)
       PRINT*,'EXECUTION INTERROMPUE'
       CALL EXIT
  9    CONTINUE
C
       N1=N+1
       N11=N+2
       N2=N/2
       N21=N2+1
       N3=N2-1
       X0=0
       PI=2.*ACOS(X0)/N
       PI2=PI/2
       SEN(1)=.25/SIN((N2-1)*PI)+.5
       DO 1 L=2,N3
       SEN(L)=.25/SIN((N2-L)*PI)+.5
       COSE(L)=1/SIN((L-1)*PI2)
  1    CONTINUE
       DO 2 L=N21,N1
       SEN(L)=.25/SIN((N2-L)*PI)+.5
       COSE(L)=1/SIN((L-1)*PI2)
  2    CONTINUE
       COSE(N2)=1/SIN(N3*PI2)
C
  999  CONTINUE
C
             IF(N.NE.NCONTR.OR.N64.NE.NFONC) THEN
C
       N257=257
       IF(N64.GT.257) THEN
       PRINT*,'DIMENSION INSUFFISANTES DANS LA SUB. CIRAMS, N64=',N64
       PRINT*,'EXECUTION INTERROMPUE'
       CALL EXIT
       ENDIF
C
       NCONTR=N
       NFONC=N64
C
       N65=N64
       IF((N64/8)*8.EQ.N64)N65=N64+1
       N63=N64-1
       NM65=N*N65
       NM651=NM65+1
       NM66=NM65+1
       N2651=N21*N65
       N265=N2651-N65
       N6565=N65+N65
       N633=N65-1
       N365=N2651-N6565
       ENDIF
C
C
       IF(INDEX.EQ.0) GO TO 11
C
       DO 3 M=1,N64
       CS(M)=-F(M)
  3    CONTINUE
C
       DO 4 M=NM65-N633,NM65
       CC(M)=-F(M)
  4    CONTINUE
C
       DO 5 M=NM66,NM65+N65
       F(M)=CC(M-N65)
  5    CONTINUE
C
       JM1=N65+1
       JM2=JM1+N63
       DO 6 L=2,N
       DO 7 M=JM1,JM2
       CS(M)=-(F(M)+F(M-N65))*.5
  7    CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
  6    CONTINUE
C
       DO 8 L=1,NM65
       F(L)=CS(L)
  8    CONTINUE
  11   CONTINUE
C
C ****   CALCUL DE LA FONCTION EN TETA=0 ET TETA=PI
C
       DO 12 M=1,N64
       CC(M)=0
  12   CONTINUE
C
       JM1=N6565
       DO 13 L=3,N,2
       DO 14 M=1,N64
       CC(M)=CC(M)+F(M+JM1)
  14   CONTINUE
       JM1=JM1+N6565
  13   CONTINUE
C
       DO 15 M=1,N64
       F1(M)=0
  15   CONTINUE
C
       JM1=N65
       DO 16 L=2,N ,2
       DO 17 M=1,N64
       F1(M)=F1(M)+F(M+JM1)
  17   CONTINUE
       JM1=JM1+N6565
  16   CONTINUE
       DO 18 M=1,N64
       CS(M)=(F(M)+F(M+NM65))*.5
  18   CONTINUE
C
       DO 19 M=1,N64
       FN(M)=CC(M)+CS(M)+F1(M)
       F1(M)=CC(M)-F1(M)+CS(M)
  19   CONTINUE
C
C ****  CALCUL DE COEF. DE FOURIER DE LA FONCYION PONDEREE A PARTIR DE C
C ****  DE TCHEBYCHEV.
C
       
       DO 20 M=1,N2651
       CS(M)=0
20     CONTINUE
C
       JM1=1
       JM2=N64
       JM5=N65
       DO 21  L=1,N2
       DO 22 M=JM1,JM2
       CS(M)=CS(M)+F(M+JM5)
       CS(M+N65)=-F(M+JM5)
  22   CONTINUE
C
       JM3=JM1+N65
       JM4=JM3+N63
       JM6=JM5-N65
       DO 23 M=JM3,JM4
       CS(M)=-F(M+JM6)
  23   CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
       JM5=JM5+N65
  21   CONTINUE
C
       JM1=1
       JM2=JM1+N63
       JM3=0
       JM4=JM4+N63
       DO 24 L=1,N21
       DO 25 M=JM1,JM2
       CC(M)=F(M+JM3)
  25   CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
       JM3=JM3+N65
  24   CONTINUE
C
       DO 26 M=1,N64
       CS(M)=0
  26   CONTINUE
C
       DO 27 M=N265+1,N265+N64
       CS(M)=0
  27   CONTINUE
C
C ***    INVERSION POUR OBTENIR LA FONCTION PONDEREE.
C
       CALL TFINMS(N,2,N64,CC,CS,F)
C       ****   RESTITUTION DE LA FONCTION
C
C
       JM1=N2651
       JM2=N365
       JM3=1
       JM4=N64
       DO 28 L=1,N3
       SENN=SEN(L)
       DO 29 M=JM3,JM4
       CC(M)=(F(M+JM2)-F(M+N2651))*SENN
  29   CONTINUE
       JM1=JM1+N65
       JM2=JM2-N6565
       JM3=JM3+N65
       JM4=JM3+N63
  28   CONTINUE
C
       JM2=N3*N65
       JM3=1
       JM4=N64
       DO 30  L=1,N3
       DO 31 M=JM3,JM4
       CS(M+JM2)=F(M+N2651)+CC(M)
       CS(M+N2651)=F(M+JM2)-CC(M)
  31   CONTINUE
C
       JM2=JM2-N6565
       JM3=JM3+N65
       JM4=JM3+N63
30     CONTINUE
C
       DO 32 L=1,N64
       CS(L)=F1(L)
32     CONTINUE
C
       JM1=N265+1
       JM2=JM1+N63
       DO 33 M=JM1,JM2
       CS(M)=F(M)
  33   CONTINUE
C
       JM1=NM65
       JM2=JM1+N63
       DO 34 M=1,N64
       CS(M+JM1)=FN(M)
  34   CONTINUE
C
       IF(INDEX.EQ.0) RETURN
C
       JM1=N65+1
       JM2=JM1+N63
       DO 35 L=2,N1
       COSEN=COSE(L)
       DO 36 M=JM1,JM2
       CS(M)=CS(M)*COSEN
  36   CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
  35   CONTINUE
       DO 37 M=1,N64
       CS(M)=0
  37   CONTINUE
C
100    FORMAT(10X,10D12.4)
101    FORMAT(1X,'IN')
       RETURN
       END
