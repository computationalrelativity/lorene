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

       SUBROUTINE CHELES(N,INV,N64,MM,CY,CC,Y)

C
C           ROUTINE POUR LE DEVELOPPEMENT D' UNE FONCTION EN
C           EN FONCTIONS ASSOCIEES DE LEGENDRE.(ET VICE-VERSA)
C
C           LE PRINCIPE EST LE SUIVANT: ON SUPPOSE QUE LA FONCTION
C           A DEVELOPPER SOIT DEVELOPPABLE EN FONCTIONS ASSOCIEES
C           DE LEGENDRE D'ORDRE m ET DE DEGRE' j :
C                    m
C                   P  (TETA)    (P(m,l))
C                    j
C
C           ON CHERCHE LES COEFFICIENTES DU DEVELOPPEMENT
C
C
C           F(TETA)=A(l)*P(m,l ) SOMME' SUR l. , 0 < l < L
C
C           LES FONCTIONS P(m,l) ONT LA FORME:
C           P(m,l)=SIN(TETA)**m*(B(j)*COS(TETA)**j) OU j VA DE
C           ZERO A l. LA FONCTION F EST DEVELOPPEE EN COEFFICIENTS
C           DE TCHEBYTCHEV DU PREMIER GENRE OU DU 2me SELON QUE
C           m EST PAIRE OU IMPAIRE. LA ROUTINE CALCULE LES MA-
C           TRICES DE TRANSFORMATION DIRCTES ET INVERSES DES COEF-
C           FICIENTS DE TCHEBYTCHEV AU COEFFICIENTS DE LEGENDRE.
C           LE STOCKAGE EST EN "PARALLEL", (CFR LA ROUTINE
C           TFMS.)
C
C
C           LA CONVENTION DES INDICES EST LA SUIVANTE:
C           LA FONCTION ASSOCIEE DE LEGENDRE D'ORDRE
C           m ET DE DEGREE j DEFINIE DANS LA LITTERATURE
C           
C                   m
C                  P
C                   j
C
C           AVEC m.LE.j m.GE.0, j.GE.0 EST STOCKEE DANS UN TA-
C           BLEAU AYANT LES INDICES M=m+1, J=j-m+1.
C
C      NOTE:
C      ------
C
C           ON DOIT REMARQUER QUE LA REALATION DE PASSAGE DE LA-
C           GRANGE A TCHEBYTCHEV N'EST PAS UNIVOQUE POUR m>0.
C           IL FAUT EN GENERALE M+L COEFFICIENTS DE TCHEBYTCHEV
C           POUR DEVELOPPER UNE FONTION DE LEGENDRE D'ORDRE M
C           ET DEGREE L. VICE-VERSA ETANT DONNEA N COEFFICIENTS
C           DE TCHEBYTCHEV ON AURA UNE CORRESPONDANCE UNIVOQUE
C           AVEC N-M F COEFFICIENTS DE LEGENDRE.
C
C               ARGUMENTS DE LA ROUTINE:
C                       ...........................
C      
C           N   = NOMBRE DES DEGRES DELIBERTE-1
C           INV = PARAMETRE: SI INV=0 LA ROUTINE EFFECTUE
C                 LE PASSAGE TCHEBYTCHEV LEGENDRE, SI INV=1
C                 ON A LA TRANSFORMATION LEGENDRE-TCHEBYTCHEV
C                 SI INV>1 ON A LA TRANSFORMATION LEGENDRE
C                 ESPACE DES TETA.
C
C           N64 = NOMBRE DES FONCTIONS QUI DOIVENT ETRE TRANS-
C                 FORMEES SIMULTANEUMENT.
C           MM  = ORDRE DE LA FONCTION DE LEGENDRE.
C           CY  = IMPUT. SI INV=0 CY DOIT CONTENIR LE COEF-
C                CIENTS DE TCHEBYTCHEV DE LA FONCTION A TRAN-
C                SFORMER. SI INV > 0 DANS CY DOIVENT ETRE STO-
C                CQUES LES COEFF. DE LEGENDRE.
C           CC  = TABLEAU DE TRAVAIL. LES DIMENSION DE CE TA-
C                 BLEAU DOIVENT ETRE > (N1+3)*(N64+1)
C           Y   = OUTPUT. LES DIMENSIONS DE CY ET Y DOIVENT
C                 ETERE > (N1+3)*(N64+1).
C
C           ROUTINE ayant testee le 26/5/1986,Routine corrigee
C	    le 26/06/90  en ayant trouve une erreur dans le cas
C	    INV=1 dans le calcul de la partie antisymmetrique.
C
C		version IBM.
C
	IMPLICIT NONE
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
c Revision 1.2  1997/05/23  11:33:21  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:35:05  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/


	double PRECISION
     1  X2,SQ2,WX,WPIG,WPI,WA1,WA2,WP,Y,CS,CC,OUT,WY,WPLJ,TETA
     1  ,Y1,WA3,WA5,WA4,WA6,WA,CY

	INTEGER NDIM,N66,M33,M,MM,N,NL,L,IM,N1,N2,N21,NJ1,NJ2
     1	,MP,J,IS,I,I2,J2,JJ,N64,NFON,NM65,NM651,MP1,INV,LMAX
     1  ,LMIN,LF,L2,JN,IJ,JN65,J0,LN64,MR1,MR2,LN65,J1,NJ21
     1	,MM1,MM2,N65
C
	DIMENSION WA1(65,65,65),WA2(65,65,65),WA3(65,65,33),CC(*)
	DIMENSION CY(*),Y1(132),WP(130,130),WY(132),IM(132),CS(132),Y(*)
	DIMENSION WA4(65,65,33),WA5(65,65,33),WA6(65,65,33)
C
	save X2,SQ2,IM,N1,N2,N21,NL,WX,WPIG,WPI
     1  ,NJ1,NJ2,NJ21,MP
     1  ,NFON,N65,NM65,NM651
C
	save JJ,WP,WA1,WA2,WA3,WA4,WA5,WA6
     1  ,MM1,MM2
C
	DATA NFON,NL/0,0/
C
C           PREPARATION DES PARAMETRES NECESSAIRES AU CALCUL.
C
       NDIM=130
       N66=129
       M33=65
       M=MM+1
C
       IF(N.GE.N66.OR.M.GT.M33) THEN
       PRINT 800,N,M
  800  FORMAT(10X,'DIMENSIONS INSUFISANTES DANS LA ROUTINE CHELES',
     , 'N=',I4,'  M=',I4)
             CALL EXIT
             ENDIF
       IF(N.EQ.NL) GO TO 22
       X2=2
       SQ2=SQRT(X2)
       DO 1 L=1,M33
       IM(L)=0
  1    CONTINUE
       N1=N+1
       N2=N/2
       N21=N2+1
       NL=N
       WX=0
       WPIG=ACOS(WX)
       WPI=WPIG/N
  22   CONTINUE
       IF(IM(M).EQ.314) GO TO 21

       DO L=1,N1
       Y1(L)=0
       ENDDO

       IM(M)=314
       NJ1=N-M+2
       NJ2=(NJ1-1)/2
       NJ21=NJ2+1 
       MP=1
       IF((M/2)*2.EQ.M)MP=0
C
C      PRINT*,'PREPARATION DES MATRICES DE TRANSFORMATION'
C
C
       DO 2 J=1,N21
       DO 3 L=1,N21
       WA1(L,J,M)=0
       WA2(L,J,M)=0
  3    CONTINUE
  2    CONTINUE
C
C
C           PREPARATION DES FONCTION DE LEGENDRE POUR m=M ET POUR
C           j.LE.(N1-M). LES FONCTIONS SONT ECHANTILLONNEES DANS
C           L'INTERVAL 0 < TETA < PI.
C
C           NORMALISATION DES FONCTIONS DE LEGENDRE.
C
       CALL LEGE1(N,NDIM,M,WP)
C
       DO 8 J=1,NJ1
       DO 7 L=1,N1
       WP(L,J)=WP(L,J)/(N*M)
  7    CONTINUE
  8    CONTINUE
C
       DO 9 J=1,NJ1
       DO 10 L=1,N1
       Y(L)=WP(L,J)**2
  10   CONTINUE
       CALL CERARS(N,0,Y,CS,CC)
       CALL INTRAS(N,CC,OUT)
       WY(J)=OUT
  9    CONTINUE
C
       DO 11 J=1,NJ1
       DO 12 L=1,N1
       WPLJ=WP(L,J)/SQRT(WY(J))
       WP(L,J)=WPLJ
  12   CONTINUE
  11   CONTINUE
C
C
C           PREPARATION DES MATRICES DE TRANSFORMATION
C           TCHEBYTCHEV-LEGENDRE. LA MATRICE EST DEFINIE
C           PAR:  SOMME DE P(TETA,l,m)*T(TETA,j)*sin(teta) d(TETA)
C
C           CALCUL POUR LA PARTIE PAIRE DE LA FONCTION
C
       IS=1
       DO 13 I2=1,N21
       IS=-IS
       I=I2+I2-1
       IF(MP.EQ.1) THEN
       DO 14 L=1,N1
       TETA=(L-1)*WPI
       Y1(L)=COS(TETA*(I-1))*IS
  14   CONTINUE
       ENDIF
C
       IF(I.EQ.1.OR.I.EQ.N1)THEN
             DO 15 L=1,N1
             Y1(L)=Y1(L)*.5
  15   CONTINUE
             ENDIF
C
       IF(MP.EQ.0) THEN
       Y1(1)=0
       DO 16 L=1,N1
       TETA=(L-1)*WPI+WPIG
       Y1(L)=-SIN(TETA*I)
  16   CONTINUE
       ENDIF
C
       DO 17 J2=1,NJ21
       J=J2+J2-1
       DO 18 L=1,N1
       Y(L)=WP(L,J)
  18   CONTINUE
C
       DO 19 L=1,N1
       Y(L)=Y(L)*Y1(L)
  19   CONTINUE
       CALL CERARS(N,0,Y,CS,CC)
       CALL INTRAS(N,CC,OUT)
       WA1(I2,J2,M)=-OUT*SQ2
  17   CONTINUE
  13   CONTINUE
C
C           CALCUL DE LA MATRICE DE TRANSFORMATION POUR LA PARTIE
C           IMPAIRE DE LA FONCTION.
C
       DO 20 I2=1,N2
       I=I2+I2
       IF(MP.EQ.1) THEN
       DO 24 L=1,N1
       TETA=(L-1)*WPI+WPIG
       Y1(L)=COS(TETA*(I-1))
  24   CONTINUE
       ENDIF
C
       IF(MP.EQ.0) THEN
       DO 25 L=1,N1
       TETA=(L-1)*WPI+WPIG
       Y1(L)=SIN(TETA*I)
  25   CONTINUE
       ENDIF
C
       DO 26 J2=1,NJ2
       J=J2+J2
       DO 27 L=1,N1
       Y(L)=WP(L,J)
  27   CONTINUE
C
       DO 28 L=1,N1
       Y(L)=Y(L)*Y1(L)
  28   CONTINUE
C
       CALL CERARS(N,0,Y,CS,CC)
       CALL INTRAS(N,CC,OUT)
       WA2(I2,J2,M)=OUT*SQ2
  26   CONTINUE
  20   CONTINUE
C
C
C           CALCUL DES MATRICES POUR LE PASAGE DE LEGENDRE A TCHEBY-
C           TCHEV
C
C
C           CES MATRICES SONT DEFINIES PAR :
C           TRANSFORMATION DE TCHEBYTCHEV DES FONCTIONS DE
C           LEGENDRE. IL-Y-A QUATTRE MATRICES DEUX POUR MM PAIRE
C           OU IMPAIRE, ET DEUX POUR LES DEUX PARITES DE LA 
C           FONCTION A TRANSFORMER.
C
C
       DO 29 J=1,NJ21
       JJ=J+J-1
       DO 30 L=1,N1
       Y(L)=WP(L,JJ)
  30   CONTINUE
C
       IF(MP.EQ.1) THEN
       MM1=(M+1)/2
       CALL CERARS(N,0,Y,CS,CC)
       DO 31 L=1,N21
       WA3(J,L,MM1)=CC(L)/SQ2
  31   CONTINUE
       ENDIF
C
       IF(MP.EQ.0) THEN
       MM2=M/2
       CALL CERA2S(N,0,Y,CS,CC)
       DO 32 L=1,N21
       WA5(J,L,MM2)=CC(L)
  32   CONTINUE
       ENDIF
  29   CONTINUE
C
       DO 33 J=1,NJ2
       JJ=J+J
       DO 34 L=1,N1
       Y(L)=WP(L,JJ)
  34   CONTINUE
       IF(MP.EQ.1) THEN
       CALL CERARS(N,1,Y,CS,CC)
       DO 35 L=1,N21
       WA4(J,L,MM1)=CC(L)/SQ2 
  35   CONTINUE
       ENDIF
C
       IF(MP.EQ.0) THEN
       CALL CERA2S(N,1,Y,CS,CC)
       DO 36 L=1,N21
       WA6(J,L,MM2)=CC(L)/SQ2
  36   CONTINUE
       ENDIF
  33   CONTINUE
  21     CONTINUE
C
       IF(N64.EQ.NFON) GO TO 23
       NFON=N64
       N65=N64
       IF((N64/8)*8.EQ.N64) N65=N64+1
       NM65=N1*N65
       NM651=N*N65+1
  23   CONTINUE
C
C.......................................................................
C  ......
C
C           COMMENCEMENT DE LA TRANSFORMATION TCHEBYTCHEV-LEGENDRE.
C
       MP=1
       IF((M/2)*2.EQ.M)MP=0
       NJ1=N-M+2
       NJ2=(NJ1-1)/2   
       NJ21=NJ2+1
       MP1=MP+1
C
             IF(INV.EQ.0) THEN
C
C           TRANSFORMATION POUR LA PARTIE SYMMETRIQUE  DE LA FON-
C           TION.
C
             LMAX=N21
C
       IF(MP.EQ.0) LMAX=N2
       DO 37 J=1,NJ21
       DO 38 LF=1,N64
       WY(LF)=0
  38   CONTINUE
       LMIN=1
       IF(M.LT.3) LMIN=J
       DO 39 L=LMIN,LMAX
       L2=(L+L-MP1)*N65
       WA=WA1(L,J,M)
       DO 40 LF=1,N64
       WY(LF)=WY(LF)+WA*CY(LF+L2)
  40   CONTINUE
  39   CONTINUE
       JN=(J+J-2)*N65
       DO 41 LF=1,N64
       Y(LF+JN)=WY(LF)
  41   CONTINUE
  37   CONTINUE
C
C           TRANSFORMATION POUR LA PARTIE ANTISYMMETRIQUE.
C
       IJ=0
       IF(MP.EQ.0) IJ=1
       DO 42 J=1,NJ2
       DO 43 L=1,N65
       WY(L)=0
  43   CONTINUE
C
       LMIN=1
       IF(M.LT.3)LMIN=J
       DO 44 L=LMIN,N2
       L2=(L+L+IJ-1)*N65
       WA=WA2(L,J,M)
       DO 45 LF=1,N64
       WY(LF)=WY(LF)+WA*CY(LF+L2)
  45   CONTINUE
  44   CONTINUE
       JN65=(J+J-1)*N65
       DO 46 LF=1,N64
       Y(LF+JN65)=WY(LF)
  46   CONTINUE
  42   CONTINUE
C
       IF(NJ1.EQ.N1) RETURN
C
C           SI m.ne.0 IL-Y-A MOINS DES COEFFICIENTES DE LEGENDRE
C           QUE DE COEFFICIENTS DE TCHEBYTCHEV. (VOIR NOTE INTRO-
C           DUCTIVE). LE TABLEAU DE SORTIE EST MIS=0 PUR J> NJ1-M.
C
       J0=NJ1
       IF(MP.EQ.1) J0=NJ1+1
       DO 47 L=J0,N1
       LN64=(L-1)*N65
       DO 48 LF=1,N64
       Y(LF+LN64)=0
  48   CONTINUE
  47   CONTINUE
       RETURN
             ENDIF
C
C
C
C           TRANSFORMATION LEGENDRE-TCHEBYTCHEV ET LEGENDRE ESPACE
C           DES CONFIGUREATIONS
C
       IF(INV.GE.1) THEN
C
       MR1=(M+1)/2
       MR2=MR1-1
       IF(MP.EQ.1) THEN
C
C           TRANSFORMATION POUR m PAIRE ET DE LA PARTIE PAIRE
C           DE LA FONCTION.
C
       MM1=(M+1)/2
       DO 49 J=1,N21
       JJ=J+J-1
       JN65=(JJ-1)*N65
       LMIN=1
       IF(J.GT.MR1) LMIN=J-MR2
       DO 50 L=1,N64
       WY(L)=0
  50   CONTINUE
       DO 51 L=LMIN,NJ21
       WA=WA3(L,J,MM1)
       LN65=(L+L-2)*N65
       DO 52 LF=1,N64
       WY(LF)=WY(LF)+WA*CY(LF+LN65)
  52   CONTINUE
  51   CONTINUE
C
       DO 53 LF=1,N64
       Y(LF+JN65)=WY(LF)
  53   CONTINUE
  49   CONTINUE
C
C
C           TRANSFORMATION DE LA PARTIE IMPAIRE DE LA FONCTION
C
       DO 54 J=1,N21
       JJ=J+J
       JN65=(JJ-1)*N65
       LMIN=1
       IF(J.GT.MR1) LMIN=J-MR2
       DO 55 LF=1,N64
       WY(LF)=0
  55   CONTINUE
C
       DO 56 L=LMIN,NJ21
       LN65=(L+L-1)*N65
       WA=WA4(L,J,MM1)
       DO 57 LF=1,N64
       WY(LF)=WY(LF)+WA*CY(LF+LN65)
  57   CONTINUE
  56   CONTINUE
C
       DO 58 LF=1,N64
       Y(LF+JN65)=WY(LF)
  58   CONTINUE
  54   CONTINUE
C
       DO 59 L=NM651,NM65
       Y(L)=Y(L)*2
  59   CONTINUE
       ENDIF
C
C
C           TRANSFORMATION POUR m PAIRE.
C
             IF(MP.EQ.0) THEN
C
C
C
       MM2=M/2
C           TRANSFORMATION DE LA PARTIE SYMMETRIQUE.
C
C
C
       DO 60 J=1,N21
       JN65=(J+J-1)*N65
       LMIN=1
       IF(J.GT.MR1) LMIN=J-MR2
       DO 61 LF=1,N64
       WY(LF)=0
  61   CONTINUE
       DO 62 L=1,NJ21
       LN65=(L+L-2)*N65
       WA=WA5(L,J,MM2)/SQ2
       DO 63 LF=1,N64
       WY(LF)=WY(LF)+CY(LF+LN65)*WA
  63   CONTINUE
  62   CONTINUE
C
       DO 64 LF=1,N64
       Y(LF+JN65)=WY(LF)
  64   CONTINUE
  60   CONTINUE
C
C
C           TRANSFORMATION DE LA PARTIE ANTISYMMETRIQUE
C
       DO 65 J=2,N21
       J1=J-1
       LMIN=J1-MR2
       IF(LMIN.LT.1) LMIN=1
       JN65=(J+J-2)*N65
       DO 66 LF=1,N64
       WY(LF)=0
  66   CONTINUE
C
       DO 67 L=LMIN,NJ21
       LN65=(L+L-1)*N65
       WA=WA6(L,J,MM2)
       DO 68 LF=1,N64
       WY(LF)=WY(LF)+WA*CY(LF+LN65)
  68   CONTINUE
  67   CONTINUE
C
       DO 69 LF=1,N64
       Y(LF+JN65)=WY(LF)
  69   CONTINUE
  65   CONTINUE
C
	DO 80 L=1,N64
	Y(L)=0
   80	CONTINUE
       ENDIF
C
             ENDIF
C
C
       IF(INV.EQ.1) RETURN
C
C           TRANSFOMATION TCHEBYTCHEV-ESPACE DES CONFIGURATIONS.
C
       IF(MP.EQ.1) CALL CHINMS(N,N64,Y,CC,CY)
       IF(MP.EQ.0) CALL CHI2MS(N,N64,Y,CC,CY)
       DO 70 L=1,NM65
       Y(L)=CY(L)
  70   CONTINUE
  100  FORMAT(1X,10E10.3)
  101  FORMAT(1X,' ')
  111  FORMAT(1X,10E10.2)
  300  FORMAT(1X,4I4,5E10.3)
       RETURN
       END
