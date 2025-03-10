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

	SUBROUTINE FCIQ3S(NDL,NDR,NDT,NDF,INDD,IN1,IMP,IDR,IND
     1  ,C64,CC,CS,DEN,DENT,ERRE,DENN)
C
C	version du 14.11.1996
C
		IMPLICIT NONE
C
C		ROUTINE POUR LE CALCUL DES TRANSFORMEES INVERSES DE FOURIER
C		E DE TCHEBYTCHEV POUR L'ECHANTILLONAGE REAREFIE' A 3 DIM.
C		DANS LE CAS MULTIZONE. CETTE ROUTINE PERMETS DE TRAITER
C		LES CAS AVEC DIFFERENTS TYPES D'ECHANTILLONAGE.
C		FIN OU RARFIE DANS LE NOYAU CENTRALE, EN u=1/r DANS
C		COQUILLES AINSI QUE DANS LA DERMIERE COQUILLE COMPACTIFIEE
C		EXPLOITANT LE SIMMETRIES DES FONCTIONS A TRANSFORMER (SYM-
C		METRIES EXISTANTES PAR EX. EN COORDONNES SPHERIQUES)
C
C		LE 3ME INDICE EST SUPPOSE' ETRE PERIODIQUES (COORDONNE AZI-
C		MUTHALE FI), LE 2ME VARIE ENTRE 0 ET PI. DANS UN DEVELOPPEMENT
C		EN COORDONNES SPHERIQUES, LES COEFFICIENTS DE FOURIER D'UNE
C		FONCTION SCALAIRE SERONT DE FONCTIONS SYMMETRIQUES EN TETA SI
C		LES COEFFICIENTS  DE FOURIER SONT PAIRES ET ANTISYMMETRIQUES
C		SONT IMPAIRS. LE DEVELOPPEMENT EN TETA EST DONC EFFECTUE' EN
C		POLYNOMES DE TCHEBYTCHEV DU 1ER GENRE POUR LES FONCTIONS
C		PAIRES ET DU 2ME GENRE (EN SERIE DE SINUS) POUR LES FONCTIONS
C		IMPAIRES. ANALOGUEMENT LES COEFFICIENTS Cml(r) SONT DES
C		FONCTIONS SYMMETRIQUES EN r SI m+l EST PAIRE ET
C		ANTISYMMETRIQUES
C		DANS LE CAS OPPOSE'. PARCONSEQUANT LE DEVELOPPEMENT EN r
C		DES COEFFICIENTS  Cml(r) EST EFFECTUE SUR L'INTERVALLE
C		0<r<1 EN TENANT COMPTE DE LA PARITE'. LA TRANSFORMATION
C
C	N.B.	DOIT ETRE PARCONSEQUENT ORDONNEE, C'EST A DIRE IL FAUT
C	---	D'ABORD PROCEDER A LA TRANSFORMATION DE TCHEBYTCHEV  SUR
C		LA VARIABLE r (1ERE INDICE), PUIS A LA TRANSFORMATION EN TETA
C		(2ME INDICE) ET ENFIN LA TRANSFORMATION EN FI (3ME INDICE).
C		CELA PEUT ETRE GENANT. LA ROUTINE FGIR3S EVITE CET INCO-
C		VENIENT.
C
C			LE STOCKAGE DES COEFFICIENTS EST LE SUIVANT (Cfr.
C		LA ROUTINE FUCE3S). DANS DEN(LR,LT,1) IL Y A LE COEFFICIENT
C		CORRESPONDANTS A LA FREQUENCE ZERO DU DEVELOPPEMENT EN COSINUS
C		DANS DEN(LR,LT,2),DEN(LR,LT,3) LES COFFICIENTS COSINUS ET
C		SINUS DE LA FREQUENCE 1, DANS DEN(LR,LT,4), DEN(LR,LT,5)
C		LES COEFFICIENTS DE LA FREQUENCE 2 ET AINSI DE SUITE.
C			 DEN(LR,1,LM),DEN(LR,3,LM).... DEN(LR,2*n+1,LM)
C		SONT LES COEFFICIENTS DU DEVELOPPEMENT SUR LES POLYNOME
C		DE TCHEBYTCHEV DU 1ER ORDRE. DEN(LR,2,LM),DEN(LR,4,LM)...
C		DEN(LR,2*n,LM) LES COEFFICIENTS DE TCHEBYTCHEV DU 2M ORDRE.
C		DEN(1,LR,LM),DEN(3,LR,LM).... SONT LES COEFFICIENTS
C		DU DEVELOPPEMENT EN POLYNOMES DE CHEBYTCHEV DES FON-
C		CTIONS SYMMETRIQUES EN r, ET DEN(2,LR,LM),DEN(4,LR,LM)....
C		DES FONCTIONS ANTISYMMETRIQUES.
C
C		ARGUMENTS DE LA ROUTINE:
C
C			NDL 	=TABLEAU CONTENENT LES PARAMETRES DES DIF-
C				 FERENTES COQUILLES: NDL(1) EST LE NOMBRE
C				 DES COQUILLES, NDL(2),NDL(3),...NDL(NZON+1)
C				 LE NOMBRE DES DEGRES DE LIBERTE EN r DE LA
C				 1ERE 2ME,NZON-EME COQUILL, NDL(NZON+2),
C				 NDL(NZON+3) LES NOMBRES DES DEGRES DE LIBERTE
C				 EN THETA ET PHI.
C		NDR,NDT,NDF	=DIMENSIONS DES DIFFERENTS
C				 TABLEAUX COMME DECLARE DANS LE PROGRAMME AP-
C				 PELLANT.
C
C		POUR DES RAISONS DE CRAYTINISATION NDR,NDT,NDF
C		NE DOIVENT PAS ETRE UN MULTIPLE DE 8.
C
C		INDD	= TABLEAU DEFINISSANT LE TYPE D'ECHANTILLONAGE
C			  DANS CHAQUE COQUILLE: INDD(LZON)=1 ECHANTILLOMAGE
C			  FIN, INDD(LZON) =0 ECHANTILLONAGE RAREFIE (SEUELEMENT
C			  POUR LZON=1, NOYEAU CENTRALE)  INDD(LZON)=3
C			  ECHANTILLONAGE EN u=1/r
C			  SONT VECTORIZEE.
C
C
C		IN1	=PARAMETRE, SI IN1=1 LA TRANSFORMEE EST
C			 EFFECTUEE SUR LE PREMIER INDICE, SI IN1=2 SUR LE
C			 DEUXIEME, SI IN1=3 SUR LE 3me,.
C		IMP	=PARAMETRE DEFINISSANT LA TENSORIALITE' DE LA FONCTION
C			 A TRANSFORMER: IMP=0 SCALAIRE, IMP=1 VECTEUR.
C
C		IDR	=PARAMETRE: DEFINISSANT LES DIFFERENTS OPERATEURS AP-
C			 PLIQUES AUX COEFFICIENTS DE LA FONCTION D'ENTRE'
C			 SELON LETABLEAU SUIVANT:
C
C					SORTIE EN OUTPUT (POUR IND=1)
C					----------------
C		         IDR=0  IN1=1	TRANSFORMATION INVERSE DE TCHEBYTCHEV
C					COEFF.-> r. (IND ARBITRAIRE)
C			 IDR=1 IN1=1   CALCUL DES COEFFICIENTS DE TCHEBYTCHEV
C				        DE L'OPERATEUR 1/r*D/DR POUR LES FON-
C					CTIONS PAIRES ET DE 1/r*D/DR-1/r**2
C					pour les fonctions impaires en r.
C					Dans la zone compactifiee, retourne
C					d/du = -r^2 d/dr (u = 1/r)
C
C			 IDR=2  IN1=1	CALCUL DE D2/DR2
C					Dans la zone compactifiee, retourne
C					u^2 d^2/du^2  (u=1/r)
C
C			 IDR=3  IN1=1   COEFFICIENTS DES LA FONCTION APRES
C					MULTIPLICATION PAR r**2
C			 IDR=4 IN1=1    COEFFICIENTS DE LA FONCTION APRES
C					DIVISION PAR r**2.
C
C			 IDR=5 IN1=1    DIVISION PAR r. (SEULEMNT POUR IND=1)
C			 IDR=6 IN1=1    MULTIPLICATION PAR r.(SEUL.POUR IND=1)
C
C			 IDR=7 IN1=1    DERIVATION PAR RAPPORT A r.(IND=1)
C				        EN OUTPUT DANS LA ZONE COMPACTIFIEE
C			                ON A 1/r * D/dr
C			 IDR=8 IN1=1    DERIVATION PAR RAPPORT A r, DANS LA
C				        ZONE COMPACTIFIEE ON D/du OU u=1/r
C
C	ATTENTION !
C       -----------
C			SI IN1=1 ET 2 > IDR < 7 LES OPERATIONS INDIQUEES SONT
C			EXECUTESS SEULEMENT DANS LE NOYEAU SPHERIQUE.
C---------------------------------------------------------------------------
C
C			 IDR=0,IN1=2	TRANSFORMATION INVERSE DE TCHEBYTCHEV
C					DU 1ER OU 2ME TYPE COEFF.-> teta.(IND
C					ARBITRAIRE)
C			 IDR=1 IN1=2	CALCUL DES COEFFICIENTS DE LA FONCTION
C					APRES L'APPLICATION DE L'OPERATEUR
C					COS(TETA)/SIN(TETA)*D/DTETA POUR LES
C					FONCTIONS DEVELLOPPES EN POLYNOMES
C					DE TCHEBYTCHEV DU 1ER TYPE ET DE
C					L'OPERATEUR
C					COS(TET)/SIN(TET)*D/DTET-1/SIN(TET)**2
C					POUR LES FONCTIONS DEVELOPPEES EN POLY-
C					NOMES DU 2ME TYPE.
C			IDR=2 IN1=2     COEFFICIENTS APRES APPLICATION DE L'O-
C					PERATEUR D2/DTET2
C			IDR=3 IN1=2     MULTIPLICATION PAR SIN(TETA)**2
C
C			IDR=4 IN1=2	DIVISION PAR SIN(TETA)**2
C
C				 POUR LES TRANSFORMES DE LEGENDRE :
C
C		LA CONVENTION DES INDICES EST LA SUIVANTE:
C		LA FONCTION ASSOCIEE DE LAGRANGE D'ORDRE
C		m ET DE DEGREE j DEFINIE DANS LA LITTERATURE
C
C				 m
C				P
C				 j
C
C		AVEC m.LE.j m.GE.0, j.GE.0
C		EST STOCKEE DANS UN TABLEAU AYANT LES INDICES
C			M=m+1, J=j-m+1.	POUR m PAIR
C 		ET
C			M=m+1, J=j-m.	POUR m IMPAIR
C
C			IDR=5 IN1=2	TRANSFORMATION TCEBYTCHEV -> FONCTION
C								m
C					ASSOCIEES DE LEGENDRE P  (teta) POUR
C								l
C					m > 0.
C			IDR=6 IN1=2     TRANSFORMATION INVERSE DE LA
C			PRECEDENTE.
C
C			IDR=7 IN1=2     TRASFORMATION TCHEBYTCHEV -> LEGENDRE
C				       	CAS m=0 inclu.
C			IDR=8 IN1=2     TRANSFORMATION INVERSE DE LA
C			PRECEDENTE.
C			IDR=9 IN1=2     DERIVEE PAR RAPPORT A theta
C
C
C			IDR=0 IN1=3     TRANSFORMATION INVERSE DE FOURIER
C					COEFF.-> ESPACE DE FI (IND ARBITRAIRE)
C			IDR=1 IN1=3     CALCUL DES COEFFICIENTS APRES L'AP-
C					PLICATION DE L'OPERATEUR D/DFI
C			IDR=2 IN1=3	CALCUL DES COEFFICIENTS APRES L'AP-
C					LICATION DE L'OPERATEUR D2/DFI2
C
C		IND	= PARAMETRE: SI IND=1 ON A EN SORTIE LES COEFFICIENTS
C			  DE FOURIER OU DE TCHEBYTCHEV APRES L'APPLICATION
C			  DES OPERATEURS DEFINIS PLUS HAUT, SI IND=2 ON A
C			  EN OUTPUT LES FONCTIONS APRES L'APPLICATION
C			  DES OPERATEURS DEFINIS PLUS HAUT A L'IMPUT.
C
C		N.B.		LA TRANSFORMATION DANS L'ESPACE DE CONFI-
C               ---	        GURATIONS EST CORRECTE SEULEMENT POUR LES
C				TRANSFORMATIONS QUI RESPECTENT LA PARITE'
C				LES CAS IDR=5,6,7 ET IN1=1 NE RESPECTENT
C				PAS LA PARITE' ET PAR CONSEQUENT SEULEMENT
C				LES COEFFICIENTS SONT CALCULES. (IND=2 ESCLU)
C
C		C64,CC,CS= TABLEAUX DE TRAVAIL: DIMENSION MINIME=
C			   NDT)*((MAX(NDEG(1),NDEG(2))+3)
C
C		DENT	= TABLEAU DE TRAVAIL 3DIM. DIMENSIONS MIN. NDEG(1)
C			  NDEG(2),NDEG(3)/2+1
C		ERRE	= TABLEAU CONTENENT LA VARIABLE r
C		DENN	=TABLEAU A 4 DIMENSIONS CONTENANT LA FONCTION
C			 A TRANSFORMER EN IMPUT, ET LA TRANSFORMEE EN
C			 OUTPUT.
C
C		Routine ayant testee avec le protocol ordinaire le 28/3/1987
C
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
c Revision 1.2  1997/05/23  11:43:48  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/14 23:14:55  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

	INTEGER NDR,NDT,NDF,NDL,NZON,NR1,NR,NY1,NF,N64,LN64,INDD,IMP,IDR
     1	,IND,NDEG,IN1,LSZ,LF,LY,LR,LZON
C
	PARAMETER (LN64=64)


	double PRECISION DEN,DENT,DENN,ERRE,CC,C64,CS,R1,R2,S2,S1,X1

	DIMENSION NDL(*)
	DIMENSION DEN(NDR,NDT,*),C64(*),CC(*),CS(*),NDEG(3),INDD(*)
	DIMENSION DENN(NDR,NDT,NDF,*),ERRE(NDR,*)
	DIMENSION DENT(NDR,NDT,*)
C
	NZON=NDL(1)
	NY1=NDL(NZON+2)
	NF=NDL(NZON+3)
	NDEG(2)=NY1
	NDEG(3)=NF
C
	IF(IN1.EQ.1) THEN

	IF ( (NZON.GT.1).AND.(IDR.GE.3).AND.(IDR.LE.6) ) THEN
	WRITE (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	WRITE (*,*) 'FCIQ3S: ATTENTION: l''operation demandee ne sera'
	WRITE (*,*) '   effectuee que dans le noyau central'
	WRITE (*,*) ' IN1, IDR : ', IN1, IDR
	WRITE (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	ENDIF

	N64=MIN(NY1*NF,LN64)
	LSZ=1
C
	IF(INDD(1).EQ.0) THEN
C
	LSZ=2
	NR1=NDL(2)
	NR=NR1-1
	NDEG(1)=NR1
	R2=ERRE(NR1,1)
	R1=0
	N64=MIN(NY1*NF,LN64)
	IF(IDR.EQ.1.OR.IDR.EQ.7.OR.IDR.EQ.8) S2=1/R2
	IF(IDR.EQ.2) S2=1/R2**2
C
	DO 1 LF=1,NF
	DO 2 LY=1,NY1
	DO 3 LR=1,NR1
	DEN(LR,LY,LF)=DENN(LR,LY,LF,1)
   3	CONTINUE
   2	CONTINUE
   1	CONTINUE
C
	IF(IDR.EQ.8) THEN
	CALL FCIR3S(NDEG,NDR,NDT,N64,IN1,IMP,7,IND,C64,CC,CS,DENT,DEN)
	ELSE
	CALL FCIR3S(NDEG,NDR,NDT,N64,IN1,IMP,IDR,IND,C64,CC,CS,DENT,DEN)
	ENDIF
C
	IF(IDR.NE.0) THEN
	IF(IDR.EQ.3)  S2=R2**2
	IF(IDR.EQ.2.OR.IDR.EQ.1.OR.IDR.EQ.4)S2=1/R2**2
	IF(IDR.EQ.6.OR.IDR.EQ.7.OR.IDR.EQ.8) S2=1/R2
	IF(IDR.EQ.5) S2=R2
C
	DO LF=1,NF
	DO LY=1,NY1
	DO LR=1,NR1
	DEN(LR,LY,LF)=DEN(LR,LY,LF)*S2
	ENDDO
	ENDDO
	ENDDO
	ENDIF
C
	DO LF=1,NF
	DO LY=1,NY1
	DO LR=1,NR1
	DENN(LR,LY,LF,1)=DEN(LR,LY,LF)
	ENDDO
	ENDDO
	ENDDO
	ENDIF
C
	DO 10 LZON=LSZ,NZON
C
	N64=MIN(NY1*NF,LN64)
	NR1=NDL(LZON+1)
	NDEG(1)=NR1
	NR=NR1-1
C
	R1=ERRE(1,LZON)
	R2=(ERRE(NR1,LZON)-R1)*.5
C
	IF(IDR.EQ.1.OR.IDR.EQ.7.OR.IDR.EQ.8) S2=1/R2
	IF(IDR.EQ.2) S2=1/R2**2
C
	DO 51 LF=1,NF
	DO 41 LY=1,NY1
	DO 31 LR=1,NR1
	DEN(LR,LY,LF)=DENN(LR,LY,LF,LZON)
   31	CONTINUE
   41	CONTINUE
   51	CONTINUE
C
	IF(IDR.LE.2) CALL FUCI3S(NDEG,NDR,NDT,N64,2,1,IDR,IND,C64,CC,CS,DEN)
	IF(IDR.EQ.7.OR.IDR.EQ.8)
     1	CALL FUCI3S(NDEG,NDR,NDT,N64,2,1,1,IND,C64,CC,CS,DEN)
C
	DO 4 LF=1,NF
	DO 5 LY=1,NY1
	DO 6 LR=1,NR1
	DENN(LR,LY,LF,LZON)=DEN(LR,LY,LF)
   6	CONTINUE
   5	CONTINUE
   4	CONTINUE
C
	IF(IDR.EQ.1.OR.IDR.EQ.2.OR.IDR.EQ.7.OR.IDR.EQ.8) THEN
C@@@C

	IF(INDD(LZON).EQ.1.OR.IDR.EQ.1.OR.IDR.EQ.8) THEN
	DO 7 LF=1,NF
	DO 8 LY=1,NY1
	DO 9 LR=1,NR1
	DENN(LR,LY,LF,LZON)=DENN(LR,LY,LF,LZON)*S2
   9	CONTINUE
   8	CONTINUE
   7	CONTINUE
	ELSE
C
	IF(IDR.EQ.7) THEN
	S1=-R1/R2
	X1=-1
	DO LF=1,NF
	DO LY=1,NY1
	DO LR=1,NR1
	DEN(LR,LY,1)=DENN(LR,LY,LF,LZON)
	ENDDO
	ENDDO
C
	CALL DIRCMS(NR,NDR,NY1,0,S1,X1,DEN,DENT)
C
	DO LY=1,NY1
	DO LR=1,NR1
	DENN(LR,LY,LF,LZON)=DENT(LR,LY,1)
	ENDDO
	ENDDO
	ENDDO
	ENDIF
C
	IF(IDR.EQ.2) THEN
	S1=-R1/R2
	X1=-1
	DO LF=1,NF
	DO LY=1,NY1
	DO LR=1,NR1
	DEN(LR,LY,1)=DENN(LR,LY,LF,LZON)
	ENDDO
	ENDDO
C
	CALL DIRCMS(NR,NDR,NY1,0,S1,X1,DEN,DENT)
	CALL DIRCMS(NR,NDR,NY1,0,S1,X1,DENT,DEN)
	DO LY=1,NY1
	DO LR=1,NR1
	DENN(LR,LY,LF,LZON)=DEN(LR,LY,LF)
	ENDDO
	ENDDO
	ENDDO
	ENDIF
	ENDIF
C
	IF(IDR.EQ.1.AND.INDD(LZON).EQ.1) THEN
	DO 25 LF=1,NF
	DO 21 LY=1,NY1
	DO 22 LR=1,NR1
	DEN(LR,LY,1)=DENN(LR,LY,LF,LZON)
   22	CONTINUE
   21	CONTINUE
	CALL DIRCMS(NR,NDR,NY1,1,R1,R2,DEN,DENT)
	DO 23 LY=1,NY1
	DO 24 LR=1,NR1
	DENN(LR,LY,LF,LZON)=DENT(LR,LY,1)
  24	CONTINUE
  23	CONTINUE
  25	CONTINUE
	ENDIF
	ENDIF
  10	CONTINUE
C
	RETURN
	ENDIF
C
	IF(IN1.GT.1) THEN
C
	DO 20 LZON=1,NZON
	NR1=NDL(LZON+1)
	NDEG(1)=NR1
	N64=MIN(NR1*NY1,LN64)
C
	DO 11 LF=1,NF
	DO 12 LY=1,NY1
	DO 13 LR=1,NR1
	DEN(LR,LY,LF)=DENN(LR,LY,LF,LZON)
   13	CONTINUE
   12	CONTINUE
   11	CONTINUE
C
	CALL FCIR3S(NDEG,NDR,NDT,N64,IN1,IMP,IDR,IND,C64,CC,CS,DENT,DEN)
C
	DO 14 LF=1,NF
	DO 15 LY=1,NY1
	DO 16 LR=1,NR1
	DENN(LR,LY,LF,LZON)=DEN(LR,LY,LF)
   16	CONTINUE
   15	CONTINUE
   14	CONTINUE
   20	CONTINUE
C
	ENDIF
C
  100	FORMAT(1X,10E10.3)
  101	FORMAT(1X,' ')
	RETURN
	END

