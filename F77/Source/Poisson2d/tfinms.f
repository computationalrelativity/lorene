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

       SUBROUTINE TFINMS(N,IP,N64,CC,CS,Y)

		implicit double precision (a-h,o-z)

C
C           ROUTINE POUR LA TF INVERSE MULTIPLE.
C
C           ARGUMENTS DE LA SUROUTINE:
C
C           N   =NOMBRE DES DGRES DE LIBERTEE. (N NOMBRE PAIR=
C                A 2**p*3**q*5**r p,q,r NOMBRES ENTIERS.
C
C           IP  =PARAMETRE: SI IP=1 L'INPUT ET L'OUTPUT
C               DOIVENT ETER EN SERIE, SI IP=0 ILS DOIVENT
C               ETRE EN PARALLEL, SI IP=2 COMME POUR IP=0
C               MAIS SENS CONTROL.
C
C           N64 =NOMBRE DES FONCTIONS DONT ON VEUT CALCULER
C                LA TF INVERSE.
C           CC  =COEFFICIENTS COSINUS DES N6 FONCTIONS.
C                LES COEFFICIENTS SONT STOKES DANS LA FACON'
C                SUIVANTE:
C                    DANS CC(1),CC(2),...CC(N64) SONT STOCKES
C                LES COFICIENTS CORRESPONDANTS A k=0 DES TF
C                DES N64 FONCTIONS,
C                DANS CC(1+N64),CC(2+N64),...CC(N64+N64)
C                LES COEFFICIRNTS CORRESPONDANTS A k=1
C                    ..........................................
C                DANS CC(1+(N/2+1)*N64),CC(2+(N/2+1)*N64),...
C                CC(N64+N64*(N/2+1)) LES COEFF. CORRESPONDANTS
C               A k=N/2 POUR LA SORTIE EN PARALLEL. (CON-
C               FRONTER AVEC LA ROUTINE TFM)
C               DAN LE CAS DE LA SORTIE SERIE LES DONNEES
C               IMPUT SONT STOCKEES EN SERIE I.E
C               CC(1),CC(2),...CC(N21) POUR LA PREMIERE
C               FONCTION, CC(N21+1),CC(N21+2),...CC(N21+N21)
C               PUR LA 2ME FONCTION, CC(2*N21+1),CC(2*N21+2),..
C               CC(3*N21) POUR LA 3ME FONCTION ET AINSI DE 
C               SUITE. DE MEME POUR CS.
C               CS COEFFICIENTS SINUS. IDEM COMME POR CC.
C           Y   = OUTPUT DES FONCTIONS DANS L'ESPACE DES x.
C                 LES VALEURES DES CES FONCTIONS SONT STOCKEES
C                 COMME DANS LA SUB TF.
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:37:47  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:42:36  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

       DIMENSION Y(*),CC(*),CS(*),IFAX(64),TRIGS(1600)
       DATA NDIM/0/
       DATA NFON/0/
       DATA JJJ/0/

	save	ndim,n2,n21,jump,NFON,N65,NM643,JJJ,ifax,trigs

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
       NDIM =N
       N2=N/2
       N21=N2+1
       JUMP=N+3
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
       NM643=JUMP*N64
       IF(IP.EQ.1) GO TO 4
       IF((N64/8)*8.EQ.N64) THEN
       N65=N64+1
       IF(JJJ.EQ.1.OR.IP.EQ.2) GO TO 3
       IF(IP.EQ.4) GO TO 3
       JJJ=1
       write(*,*) 'LE NOMBRE DE FONCTIONS A TRANSFORMER '
     1	,'EST UN MULTIPLE DE 8'
       write(*,*) 'POUR DES RAISONS LES DONNES DOIVENT'
       write(*,*) 'ETRE STOCKEES AVEC UN PAS DE N64+1'
       write(*,*) 'SI ON VEUT CONTINUER LE CALCUL TAPER 1.'
C
       read(*,*) ICONT
       IF(ICONT.EQ.1) GO TO 3
       PRINT 404
 404   FORMAT(10X,'CALCUL ARRETE PARCEQUE LES DONNES ONT ETE STOCKEES'
     , ,' D UNE FACON INCORECTE')
       CALL EXIT
C
  3    CONTINUE
       ENDIF
  4    CONTINUE
C
       N63=N65-1
       NM65=N*N65
       NM652=N21*N65
C
       IF(IP.EQ.1) GO TO 5
C
C           TF INVERSE EN PARALLEL
C
C
C           STOCAGE DANS Y DES COEFF. COS ET SINUS DES N64
C           FONCTIONS. LE STOCKAGE EST EFFECTUE DANS LA FACON
C           SUIVANTE:
C           DANS Y(1),Y(2),...Y(N64) SONT STOCKES LES N64
C           VALEURES DES COEFFICIENTS COS CORRESPONDANTES A k=0.
C           DANS Y(1+N64),Y(2+N64),...Y(N64+N64) LES CORFFICIENTS
C           SINUS (POUR k=0)
C           DANS Y(1+2*N64),Y(2+2*N64),...Y(N64+2*N64) IL-Y-A
C           LES COEFF COS POUR k=1, 
C           DANS Y(1+3*N64),Y(2+3*N64),...Y(N64+3*N64) LES COEFF.
C           SINUS POUR k=1 ET AINSI DE SUITE.
C
C           TV INVERSE EN SERIE
C
       IF(IP.EQ.4) GO TO 8
       LJ=1
       DO 11 L=1,N21
       LJ2=LJ-1
       LJ265=LJ2+N65
       LJ63=LJ+N63
       DO 10 M=LJ,LJ63
       Y(LJ2+M)=CC(M)
       Y(LJ265+M)=-CS(M)
10     CONTINUE
       LJ=LJ+N65
11     CONTINUE
C
  8    CONTINUE
       CALL FFT991(Y,CC,TRIGS,IFAX,N65,1,N,N64,1)
C
       DO 12 L=1,NM65
       Y(L)=Y(L)*.5
12     CONTINUE
  101  FORMAT(1X,'TFM')
  100  FORMAT(1X,'TFM',10D12.4)
       RETURN
C
  5    CONTINUE
C      
C      TF INVERSE EN SERIE
C
       JM=0
       JN=1
       DO 14 M=1,N64
       LJN2=JN+N2
       JM2=JM-(JN+JN-2)
       DO 13 L=JN,LJN2
       LJ=L+L+JM2
       Y(LJ-1)=CC(L)
       Y(LJ)=-CS(L)
13     CONTINUE
       JN=JN+N21
       JM=JM+JUMP
14     CONTINUE
C
       CALL FFT991(Y,CC,TRIGS,IFAX,1,JUMP,N,N64,1)
       DO 15 L=1,NM643
       Y(L)=Y(L)*.5
15     CONTINUE
       RETURN
       END
