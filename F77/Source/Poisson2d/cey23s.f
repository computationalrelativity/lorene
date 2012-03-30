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
       SUBROUTINE CEY23S(NDEG,NDIMR,NDIMY,NNN64,ITCH,C64,CC,CS,DEN)

		implicit double precision (a-h,o-z)

C
C           ROUTINE POUR LA TRANSFORMATION DE LA PARTIE EN TETA
C           D'UNE FONCTION EN CORDONNES SPHERIQUES EN POLYNOMES DE
C           TCHEBYTCHEV DU 1ER OU 2ME TYPE.
C           ROUTINE HAUTEMENT SECIALISEE QUI DOIT ETRE APPELEE
C           PAR FCER3S OU FGER3S.
C
C
C           ROUTINE COMPLETEMENT CRAYTINIZEE.
C
C           ARGUMENTS DE LA ROUTINE:
C
C               NDEG = TABLEAU, NDEG(3) CONTENANT LES DE-
C                    GREES DE LIBERTE DES TRANSFORMEES
C                    A EFFECTUER, NDEG(1) CONCERNE LE PRE-
C                    MIER INDICE DE LA MATRICE, NDEG(2)
C                    LE 2ME INDICE, NDEG(3) LE 3me INDICE DE LA
C                    MATRICE.     
C                    NDEG DOIT IMPERATIVEMENT ETRE DE LA
C                    FORME 2**m*3**p*5**q POUR LES TRANS-
C                    FORMEEES DE FOURIER (m,p,q NOMBRES 
C                    INTIERS)
C                    ET 2**p*3**p*5*q+1 POUR LES TRANSFOR-
C                    MEES DE TCHEBYTCHEV.
C
C           NDIMY,NDIMZ  =DIMESION DU TABLEAU YY(LR,LY,LZ).
C           POUR DES RAISONS DE CRAYTINISATION NDIMY ET NDIMZ NE DOIT PA
C  S
C           ETRE UN MULTIPLE DE 8.
C
C           NNN64 = PARAMETRE DE LA VECTORIZATION, PAR FOR TESTEXEMPLE
C                NNN64=64 SIGNIFIE QUE 64 FONCTIONS A TRANSFORMER
C                SONT VECTORIZEE.
C
C           ITCH =PARAMETRE, SII ITCH.EQ.1 LA TRANSFORMEE
C                DE TCHEBYTCHEV DU 1ER TYPE, SI ITCH=2 DU 2ME TYPE.
C
C           C64,CC,CS= TABLEAUX DE TRAVAIL: DIMENSION MINIME=
C                  (NNN64+1)*((MAX(NDEG(1),NDEG(2))+3)
C
C           DEN =TABLEAU A 3 DIMENSIONS CONTENANT LA FONCTION
C                A TRANSFORMER EN IMPUT, ET LA TRANSFORMEE EN
C                OUTPUT. LES COEFFICIENTS DE FOURIER (DANS LE CAS
C                D'UNE TRANSFORMATION DE FOURIER) SONT STOCKES
C                DANS LA FACON SUIVANTE: (EXEMPLE DANS LE CAS
C                DE LA TRANSFORMATION DU 1ER INDICE) DANS DEN(1,L,M)
C                IL-Y-A LE COEFFICIENT DE FOURIER COSINUS DE LA 
C                FREQUENCE ZERO, DANS LES COEFFICIENTS PAIRES 
C                (DEN(2,L,M),DEN(4,L,M)...DEN(2*N,L,M) LES COEFI-
C                CIENTS EN COSINUS, DANS LES COEFFICIENTS IMPAIRES
C                DEN(3,L,M),DEN(5,L,M),...DEN(2*N-1) LES TERMES EN 
C                SINUS. PARCONSEQUENT LE TERMES DEN(2,L,M) ET DEN(3,L,M)
C                SONT LES TERMES DU DEVELOPPEMENT DE FOURIER CORRESPON-
C                DANTS A LA MEME FREQUENCE. EN TOTALE IL-Y-A 2*N
C                DGRES DE LIBERTE.
C
C           Routine ayant testee avec le protocol ordinaire le 10/12/198
C  6
C
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
c Revision 1.1  1997/10/23  08:30:29  eric
c Initial revision
c
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/


       DIMENSION DEN(NDIMR,NDIMY,*),C64(*),CC(*),CS(*),NDEG(3)
       DATA NDY,NDR,NDZ/0,0,0/
       DATA NEQ/0/

	save	NDR,NDY,NDZ,NEQ,NR,NZ,NY,NNN65,NN64Y,NN65Y
	save	NR1Z1,MULTY,IRESTY,N6Y,MULT1,IREEE2,II,I2
	save	LZZ1,LZZ2,LDZ,N64Y,N65Y,N6565Y,N365Y,IRR2Y

       NR1=NDEG(1)
       NY1=NDEG(2)
       NZ1=NDEG(3)
C
C           INITIALISATION.
C
       IF(NDY.EQ.NY1.AND.NDZ.EQ.NZ1.AND.NDR.EQ.NR1.AND.NEQ.EQ.NNN64) GO 
     &TO 800
C           PRINT*,'INITIALISATION'
       NDR=NR1
       NDY=NY1
       NDZ=NZ1
       NEQ=NNN64
       NR=NR1-1
       NZ=NZ1-1
       NY=NY1-1
       NNN65=NNN64
       IF((NNN64/8)*8.EQ.NNN64)NNN65=NNN64+1
C
C               PREPARATION DES QUANTITES NCESSAIRES POUR LA TRANSFOR-
C           MATION DU 2me INDEX DU TABLEAU.
C
       NN64Y=NNN64
       NN65Y=NNN65
       NR1Z1=NR1*NZ1
       MULTY=(NR1Z1)/NN64Y
       IRESTY=NR1Z1-MULTY*NN64Y
C
C           OPTIMISATION DE NN64Y. (VOIR PLUS HAUT)
C
       IF(MULTY.GT.0.AND.NNN64.GT.NR1) THEN
       N6Y=(NNN64/NR1)*NR1
       MULT1=NR1Z1/N6Y
       IREEE2=NR1Z1-MULT1*N6Y
       II=0
       I2=0
       IF(IRESTY.GT.0)II=1
       IF(IREEE2.GT.0)I2=1
       IF(MULT1+I2.LE.MULTY+II) THEN
       MULTY=MULT1
       NN64Y=N6Y
       NN65Y=NN64Y
       IF((NN64Y/8)*8.EQ.NN64Y)NN65Y=NN64Y+1
       IRESTY=IREEE2
       ENDIF
       ENDIF
C
       IF(MULTY.EQ.0) THEN
       NN64Y=NR1Z1
       NN65Y=NN64Y
       IF((NN64Y/8)*8.EQ.NN64Y)NN65Y=NN64Y+1
       MULTY=1
       IRESTY=0
       ENDIF
       LZZ1=1
       LZZ2=NN64Y/NR1
       LDZ=LZZ2-LZZ1
       IF(IRESTY.GT.0) THEN
       N64Y=IRESTY
       N65Y=N64Y
       IF((N64Y/8)*8.EQ.N64Y)N65Y=N64Y+1 
       ENDIF
C
       N6565Y=NN65Y+NN65Y
       N365Y=N6565Y+NN65Y
       IRR2Y=NN64Y-(NN64Y/NR1)*NR1
C
C
 800   CONTINUE
C
C           TRANSFORMATION DU 2ME INDICE
C
C22222222222222222222222222222222222222222222222222222222222222222222222
C  2.
C22222222222222222222222222222222222222222222222222222222222222222222222
C  2.
C
C
       IF(NY.LT.2) RETURN
C
       IF(IRR2Y.EQ.0)THEN
C
C           ON EFFECTUE LA TRANSFORMATION DANS LE CAS NN64Y MULTIPLE DE 
C  NR1.
C      
       LZ1=LZZ1
       LZ2=LZZ2
C
       DO 115 LMU=1,MULTY
       JY=0
       DO  102 LY=1,NY1
       JZ=JY 
       DO 100 LZ=LZ1,LZ2
       DO 101 LR=1,NR1
       C64(LR+JZ)=DEN(LR,LY,LZ)
  101  CONTINUE
       JZ=JZ+NR1
  100  CONTINUE
       JY=JY+NN65Y
  102  CONTINUE
C
       IF(ITCH.EQ.1) CALL CHEYMS(NY,NN64Y,C64,CC,CS)
       IF(ITCH.EQ.2) CALL CEY2MS(NY,NN64Y,C64,CS,CC)
C
C           ON STOCKE LES COEFFICIENTS DE LA TRANSFORMATION DANS 
C           DEN(LR,LY,LZ).
C
C
       JY=0
       DO 114 LY=1,NY1
       JZ=JY 
       DO 112 LZ=LZ1,LZ2
       DO 113 LR=1,NR1
       DEN(LR,LY,LZ)=CC(LR+JZ)
 113   CONTINUE
       JZ=JZ+NR1
 112   CONTINUE
       JY=JY+NN65Y
 114   CONTINUE
C
       LZ1=LZ2+1
       LZ2=LZ1+LDZ
  115  CONTINUE
             IF(IRESTY.EQ.0) RETURN
C
C           LE CALCUL EST CONTINUE' SI NR1*NY1 N'EST PAS UN MULTIPLE
C           DE N64Y.
C
       IF(IRESTY.GT.0) THEN
C
       JY=0
       DO 118 LY=1,NY1
       JZ=JY
       DO 116 LZ=LZ1,NZ1
       DO 117 LR=1,NR1
       C64(LR+JZ)=DEN(LR,LY,LZ)
 117   CONTINUE
       JZ=JZ+NR1
 116   CONTINUE
       JY=JY+N65Y
 118   CONTINUE
C
       IF(ITCH.EQ.1) CALL CHEYMS(NY,N64Y,C64,CC,CS)
       IF(ITCH.EQ.2) CALL CEY2MS(NY,N64Y,C64,CS,CC)
C
C           ON REINTRODUIT LES COEFF. DE FOURIER DANS DEN(LR,LY,LZ)
C
       JY=0
       DO 130  LY=1,NY1
       JZ=JY
       DO 128 LZ=LZ1,NZ1
       DO 129 LR=1,NR1
       DEN(LR,LY,LZ)=CC(LR+JZ)
 129   CONTINUE
       JZ=JZ+NR1
 128   CONTINUE
       JY=JY+N65Y
 130   CONTINUE
C
       RETURN
       ENDIF
       ENDIF
C
C           CALCUL DE LA TF DANS LE CAS OU N64Y N'EST PAS UN MULTIPLE DE
C   
C           NR1.
C
       IF(IRR2Y.GT.0) THEN
C           
       NR64Y=NN64Y
       N63R=NR1-NR64Y+1
       NR65Y=NN65Y
       MU2=NR64Y
       LZ1=1
       LR1=1
       LRR=1
       DO 189 LM=1,MULTY
       IF(LRR.GT.NR1) LRR=1
       LRR1=0
       LZ1=(MU2-NR64Y)/NR1+1
       LZ2=.99999+FLOAT(MU2)/NR1
       LR1=+MU2+N63R-LZ1*NR1
       LR2=MU2-(LZ2-1)*NR1
       JY=0
C
       DO 136 LY=1,NY1
       JJZ=JY
       IF(LZ2.GT.LZ1) THEN
       JZ=JY-LRR+1
C
       DO 131 LR=LRR,NR1
       C64(LR+JZ)=DEN(LR,LY,LZ1)
 131   CONTINUE                               
       JJZ=NR1+JZ
       LRR1=1
       ENDIF
C
       IF(LZ2.GT.LZ1+1) THEN
       JZ=JJZ
       DO 132 LZ=LZ1+1,LZ2-1
       DO 133 LR=1,NR1
       C64(LR+JZ)=DEN(LR,LY,LZ)
 133   CONTINUE                               
       JZ=JZ+NR1
 132   CONTINUE                               
       JJZ=JZ
       ENDIF
C
       IF(LR2.GE.LR1.OR.LRR1.EQ.1) THEN
C
       IF(LRR1.EQ.0) THEN
       JZ=JJZ+1-LR1
       DO 134 LR=LR1,LR2
       C64(LR+JZ)=DEN(LR,LY,LZ2)
 134   CONTINUE                               
       ENDIF
       IF(LRR1.EQ.1) THEN
       JZ=JJZ
       DO 135 LR=1,LR2
       C64(LR+JZ)=DEN(LR,LY,LZ2)
 135   CONTINUE                               
       ENDIF
C
       LSR=LR2+1
       ENDIF
       JY=JY+NR65Y
 136   CONTINUE                               
C
       IF(ITCH.EQ.1) CALL CHEYMS(NY,NR64Y,C64,CC,CS)
       IF(ITCH.EQ.2) CALL CEY2MS(NY,NR64Y,C64,CS,CC)
C
C           LES VALEURES CALCLCULEES PAR TFMYS SONT MEMORISEES DANS 
C           DEN(LR,LY,LZ).
C
       JY=0
       DO 162 LY=1,NY1
       JJZ=JY
       IF(LZ2.GT.LZ1) THEN
C
       JZ=JY-LRR+1
       DO 157 LR=LRR,NR1
       DEN(LR,LY,LZ1)=CC(LR+JZ)
 157   CONTINUE          
       JJZ=NR1+JZ
       LRR1=1
       ENDIF
C
       IF(LZ2.GT.LZ1+1) THEN
       JZ=JJZ
       DO 158 LZ=LZ1+1,LZ2-1
       DO 159 LR=1,NR1
       DEN(LR,LY,LZ)=CC(LR+JZ)
 159   CONTINUE          
       JZ=JZ+NR1
 158   CONTINUE          
       JJZ=JZ
       ENDIF
C
       IF(LR2.GE.LR1.OR.LRR1.EQ.1) THEN
       IF(LRR1.EQ.0) THEN
       JZ=JJZ+1-LR1
       DO 160 LR=LR1,LR2
       DEN(LR,LY,LZ2)=CC(LR+JZ)
 160   CONTINUE          
       ENDIF
C
       IF(LRR1.EQ.1) THEN
       JZ=JJZ
       DO 161 LR=1,LR2
       DEN(LR,LY,LZ2)=CC(LR+JZ)
 161   CONTINUE          
       ENDIF
C
       LSR=LR2+1
       ENDIF
       JY=JY+NR65Y
 162   CONTINUE          
C
       LRR=LSR 
       MU2=MU2+NR64Y
 189   CONTINUE                               
C
       IF(IRESTY.EQ.0) RETURN
C
       JY=0
       LR1=LRR
       LZ1=LZ2
       LZ3=LZ2
       IF(LRR.GT.NR1) THEN
       LR1=1
       LZ3=LZ2+1
       ENDIF
       NR64Y=IRESTY
       NR65Y=NR64Y
       IF((NR64Y/8)*8.EQ.NR64Y)NR65Y=NR64Y+1
       N265=NR65Y+NR65Y
       JY=-LR1+1
C
       DO 166 LY=1,NY1
C
       JZ=JY
       IF(LRR.LE.NR1) THEN
       DO 163 LR=LR1,NR1
       C64(LR+JZ)=DEN(LR,LY,LZ2)
 163   CONTINUE          
       LZ3=1+LZ2
       JZ=JZ+NR1
       ENDIF
C
       IF(LZ3.LE.NZ1)THEN
       DO 164 LZ=LZ3,NZ1
       DO 165 LR=1,NR1
       C64(LR+JZ)=DEN(LR,LY,LZ)
 165   CONTINUE          
       JZ=JZ+NR1
 164   CONTINUE          
       ENDIF  
       JY=JY+NR65Y
 166   CONTINUE          
C      
       IF(ITCH.EQ.1) CALL CHEYMS(NY,NR64Y,C64,CC,CS)
       IF(ITCH.EQ.2) CALL CEY2MS(NY,NR64Y,C64,CS,CC)
C           
C           ON STOCKES LES COEFFICIENTS DE FOURIER DANS DEN(LR,LY,LZ).
C
       JY=-LR1+1
       DO 183 LY=1,NY1
C
       JZ=JY
       IF(LRR.LE.NR1) THEN
       DO 180 LR=LR1,NR1
       DEN(LR,LY,LZ2)=CC(LR+JZ)
 180   CONTINUE          
       LZ3=1+LZ2
       JZ=JZ+NR1
       ENDIF
C
       IF(LZ3.LE.NZ1)THEN
       DO 181 LZ=LZ3,NZ1
       DO 182 LR=1,NR1
       DEN(LR,LY,LZ)=CC(LR+JZ)
 182   CONTINUE          
       JZ=JZ+NR1
 181   CONTINUE          
       ENDIF  
       JY=JY+NR65Y
 183   CONTINUE          
C
       RETURN
       ENDIF
C
       RETURN
       END
C
