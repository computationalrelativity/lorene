C Fast Fourier Transform subroutines
C
C Copyright (c) 1978 Clive Temperton (European Centre for Medium-Range 
C                                     Weather Forecasts, Reading, UK)
C Copyright (c) 1980 Russ Rew (National Center for Atmospheric Research,
C                              Boulder, Colorado, USA) 
C
C References:
C ----------
C    C. Temperton: ``Self-Sorting Mixed-Radix Fast Fourier
C          Transforms.'', Journal of Computational Physics 52, 1 (1983)
C    C. Temperton: ``Fast Mixed-Radix Real Fourier Transforms.'',
C          Journal of Computational Physics 52, 340 (1983)
C               
C
C $Id$
C $Log$
C Revision 1.2  2004/10/04 19:22:00  e_gourgoulhon
C Added copyright and references.
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.3  2000/12/14  15:41:16  eric
c subroutine VPASSM (ligne 152) : les DATA sont forces a la double
c  precision par l'ajout de D00 aux valeurs numeriques
c  (cela generait des erreurs double -> simple precision avec les
c   compilateurs g77 et NAG f95 sous Linux).
c
c Revision 1.2  1997/05/23  11:45:49  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:05:32  hyc
C Initial revision
C
C
C $Header$
C
C

       SUBROUTINE FAX(IFAX,N,MODE)

       IMPLICIT double PRECISION (A-H,O-Z)

	character*100 header
	data header/'$Header$'/

       DIMENSION IFAX(*)
       NN=N
       IF (IABS(MODE).EQ.1) GO TO 10
       IF (IABS(MODE).EQ.8) GO TO 10
       NN=N/2
       IF ((NN+NN).EQ.N) GO TO 10
       IFAX(1)=-99
       RETURN
10     K=1
C       TEST FOR FACTORS OF 4
20     IF (MOD(NN,4).NE.0) GO TO 30
       K=K+1
       IFAX(K)=4
       NN=NN/4
       IF (NN.EQ.1) GO TO 80
       GO TO 20
C       TEST FOR EXTRA FACTOR OF 2
30     IF (MOD(NN,2).NE.0) GO TO 40
       K=K+1
       IFAX(K)=2
       NN=NN/2
       IF (NN.EQ.1) GO TO 80
C       TEST FOR FACTORS OF 3
40     IF (MOD(NN,3).NE.0) GO TO 50
        K=K+1
       IFAX(K)=3
       NN=NN/3
       IF (NN.EQ.1) GO TO 80
       GO TO 40
C       NOW FIND REMAINING FACTORS
50       L=5
       INC=2
C       INC ALTERNATELY TAKES ON VALUES 2 AND 4
60      IF (MOD(NN,L).NE.0) GO TO 70
       K=K+1
       IFAX(K)=L
       NN=NN/L
       IF (NN.EQ.1) GO TO 80
       GO TO 60
70     L=L+INC
       INC=6-INC
       GO TO 60
80      IFAX(1)=K-1
C       IFAX(1) CONTAINS NUMBER OF FACTORS
C       IFAX(1) CONTAINS NUMBER OF FACTORS
       NFAX=IFAX(1)
C       SORT FACTORS INTO ASCENDING ORDER
       IF (NFAX.EQ.1) GO TO 110
       DO 100 II=2,NFAX
       ISTOP=NFAX+2-II
       DO 90 I=2,ISTOP
       IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
       ITEM=IFAX(I)
       IFAX(I)=IFAX(I+1)
       IFAX(I+1)=ITEM
90      CONTINUE
100    CONTINUE
110     CONTINUE
       RETURN
       END
C
       SUBROUTINE FFTRIG(TRIGS,N,MODE)

       IMPLICIT double PRECISION (A-H,O-Z)

       DIMENSION TRIGS(*)
       X1=1
       PI=2.*ASIN(X1)
       IMODE=IABS(MODE)
       NN=N
       IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
       DEL=(PI+PI)/DFLOAT(NN)
       L=NN+NN
       DO 10 I=1,L,2
       ANGLE=.5D00*DFLOAT(I-1)*DEL
       TRIGS(I)=COS(ANGLE)
       TRIGS(I+1)=SIN(ANGLE)
10     CONTINUE
       IF (IMODE.EQ.1) RETURN
       IF (IMODE.EQ.8) RETURN
       DEL=.5D00*DEL
       NH=(NN+1)/2
       L=NH+NH
       LA=NN+NN
       DO 20 I=1,L,2
       ANGLE=.5D00*DFLOAT(I-1)*DEL
       TRIGS(LA+I)=COS(ANGLE)
       TRIGS(LA+I+1)=SIN(ANGLE)
20     CONTINUE
       IF (IMODE.LE.3) RETURN
       DEL=.5D00*DEL
       LA=LA+NN
       IF (MODE.EQ.5) GO TO 40
       DO 30 I=2,NN
       ANGLE=DFLOAT(I-1)*DEL
       TRIGS(LA+I)=2*SIN(ANGLE)
30      CONTINUE
       RETURN
40     CONTINUE
       DEL=.5D00*DEL
       DO 50 I=2,N
       ANGLE=DFLOAT(I-1)*DEL
       TRIGS(LA+I)=SIN(ANGLE)
50     CONTINUE
       RETURN
       END
C
C       SUBROUTINE 'VPASSMD' - MULTIPLE VERSION OF 'VPASSA'
C       PERFORMS ONE PASS THROUGH DATA
C       AS PART OF MULTIPLE COMPLEX FFT ROUTINE
C       A IS FIRST REAL INPUT VECTOR
C       B IS FIRST IMAGINARY INPUT VECTOR
C       C IS FIRST REAL OUTPUT VECTOR
C       D IS FIRST IMAGINARY OUTPUT VECTOR
C       TRIGS IS PRECALCULATED TABLE OF SINES ' COSINES
C       INC1 IS ADDRESSING INCREMENT FOR A AND B
C       INC2 IS ADDRESSING INCREMENT FOR C AND D
C       INC3 IS ADDRESSING INCREMENT BETWEEN A'S & B'S
C       INC4 IS ADDRESSING INCREMENT BETWEEN C'S & D'S
C       LOT IS THE NUMBER OF VECTORS
C       N IS LENGTH OF VECTORS
C       IFAC IS CURRENT FACTOR OF N
C       LA IS PRODUCT OF PREVIOUS FACTORS
C
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,
     * IFAC,LA)

       IMPLICIT double PRECISION (A-H,O-Z)

       DIMENSION A(*),B(*),C(*),D(*),TRIGS(*)
       DATA SIN36/0.587785252292473D00/,COS36/0.809016994374947D00/,
     *       SIN72/0.951056516295154D00/,COS72/0.309016994374947D00/,
     *       SIN60/0.866025403784437D00/
C
       M=N/IFAC
       IINK=M*INC1
       JINK=LA*INC2
       JUMP=(IFAC-1)*JINK
       IBASE=0
       JBASE=0
       IGO=IFAC-1
       IF (IGO.GT.4) RETURN
       GO TO (10,50,90,130),IGO
C
C       CODING FOR FACTOR 2
C
10     IA=1
       JA=1
       IB=IA+IINK
       JB=JA+JINK
       DO 20 L=1,LA
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 15 IJK=1,LOT
       C(JA+J)=A(IA+I)+A(IB+I)
       D(JA+J)=B(IA+I)+B(IB+I)
       C(JB+J)=A(IA+I)-A(IB+I)
       D(JB+J)=B(IA+I)-B(IB+I)
       I=I+INC3
       J=J+INC4
15     CONTINUE
       IBASE=IBASE+INC1
       JBASE=JBASE+INC2
20      CONTINUE
       IF (LA.EQ.M) RETURN
       LA1=LA+1
       JBASE=JBASE+JUMP
       DO 40 K=LA1,M,LA
       KB=K+K-2
       C1=TRIGS(KB+1)
       S1=TRIGS(KB+2)
       DO 30 L=1,LA
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 25 IJK=1,LOT
       C(JA+J)=A(IA+I)+A(IB+I)
       D(JA+J)=B(IA+I)+B(IB+I)
       C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
       D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
       I=I+INC3
       J=J+INC4
25     CONTINUE
       IBASE=IBASE+INC1
       JBASE=JBASE+INC2
30      CONTINUE
       JBASE=JBASE+JUMP
40     CONTINUE
       RETURN
C
C       CODING FOR FACTOR 3
C
50     IA=1
       JA=1
       IB=IA+IINK
       JB=JA+JINK
       IC=IB+IINK
       JC=JB+JINK
       DO 60 L=1,LA
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 55 IJK=1,LOT
       C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
       D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-.5D0*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-.5D0*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-.5D0*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-.5D0*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
       I=I+INC3
       J=J+INC4
55     CONTINUE
       IBASE=IBASE+INC1
       JBASE=JBASE+INC2
60     CONTINUE
       IF (LA.EQ.M) RETURN
       LA1=LA+1
       JBASE=JBASE+JUMP
       DO 80 K=LA1,M,LA
       KB=K+K-2
       KC=KB+KB
       C1=TRIGS(KB+1)
       S1=TRIGS(KB+2)
       C2=TRIGS(KC+1)
       S2=TRIGS(KC+2)
       DO 70 L=1,LA
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 65 IJK=1,LOT
       C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
       D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
       C(JB+J)=C1*((A(IA+I)-.5*(A(IB+I)+
     *A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *-S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
       D(JB+J)=S1*((A(IA+I)-0.5*(A(IB+I)+
     *A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *+C1*((B(IA+I)-.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
       C(JC+J)=C2*((A(IA+I)-.5*(A(IB+I)+
     *A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *-S2*((B(IA+I)-.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
       D(JC+J)=S2*((A(IA+I)-.5*(A(IB+I)+
     *A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *+C2*((B(IA+I)-.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
       I=I+INC3
       J=J+INC4
65     CONTINUE
       IBASE=IBASE+INC1
       JBASE=JBASE+INC2
70     CONTINUE
       JBASE=JBASE+JUMP
80     CONTINUE
       RETURN
C
C       CODING FOR FACTOR 4
C
90     IA=1
       JA=1
       IB=IA+IINK
       JB=JA+JINK
       IC=IB+IINK
       JC=JB+JINK
       ID=IC+IINK
       JD=JC+JINK
       DO 100 L=1,LA
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 95 IJK=1,LOT
       C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
       C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
       D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
       D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
       C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
       C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
       D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
       D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
       I=I+INC3
       J=J+INC4
95      CONTINUE
       IBASE=IBASE+INC1
       JBASE=JBASE+INC2
100    CONTINUE
       IF (LA.EQ.M) RETURN
       LA1=LA+1
       JBASE=JBASE+JUMP
       DO 120 K=LA1,M,LA
       KB=K+K-2
       KC=KB+KB
       KD=KC+KB
       C1=TRIGS(KB+1)
       S1=TRIGS(KB+2)
       C2=TRIGS(KC+1)
       S2=TRIGS(KC+2)
       C3=TRIGS(KD+1)
       S3=TRIGS(KD+2)
       DO 110 L=1,LA
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 105 IJK=1,LOT
       C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
       D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
       C(JC+J)= C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *  -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
       D(JC+J)= S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *  +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
       C(JB+J)=C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *  -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
       D(JB+J)=S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *  +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
       C(JD+J)=C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     * -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
       D(JD+J)=S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     * +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
       I=I+INC3
       J=J+INC4
105     CONTINUE
       IBASE=IBASE+INC1
       JBASE=JBASE+INC2
110     CONTINUE
       JBASE=JBASE+JUMP
120     CONTINUE
       RETURN
C
C       CODING FOR FACTOR 5
C
130    IA=1
       JA=1
       IB=IA+IINK
       JB=JA+JINK
       IC=IB+IINK
       JC=JB+JINK
       ID=IC+IINK
       JD=JC+JINK
       IE=ID+IINK
       JE=JD+JINK
       DO 140 L=1,LA
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 135 IJK=1,LOT
       C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
       D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *     -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *     +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *     +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *     -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *     +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *     +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *     -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
       I=I+INC3
       J=J+INC4
135    CONTINUE
       IBASE=IBASE+INC1
       JBASE=JBASE+INC2
140    CONTINUE
       IF (LA.EQ.M) RETURN
       LA1=LA+1
       JBASE=JBASE+JUMP
       DO 160 K=LA1,M,LA
       KB=K+K-2
       KC=KB+KB
       KD=KC+KB
       KE=KD+KB
       C1=TRIGS(KB+1)
       S1=TRIGS(KB+2)
       C2=TRIGS(KC+1)
       S2=TRIGS(KC+2)
       C3=TRIGS(KD+1)
       S3=TRIGS(KD+2)
       C4=TRIGS(KE+1)
       S4=TRIGS(KE+2)
       DO 150 L=1,LA
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 145 IJK=1,LOT
       C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
       D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
       C(JB+J)=
     * C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *-(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *-S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *+(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
       D(JB+J)=
     *S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *-(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *+C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *+(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
       C(JE+J)=
     * C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *+(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *-S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *-(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
       D(JE+J)=
     *S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *+(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *+C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *-(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
       C(JC+J)=
     *C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *-(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *-S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *+(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
       D(JC+J)=
     *S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     * -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     * +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     * +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
       C(JD+J)=
     * C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     * +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     * -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     * -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
       D(JD+J)=
     * S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *  +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     * -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
       I=I+INC3
       J=J+INC4
145     CONTINUE
       IBASE=IBASE+INC1
       JBASE=JBASE+INC2
150    CONTINUE
       JBASE=JBASE+JUMP
160    CONTINUE
       RETURN
       END
C
C       SUBROUTINE 'FFT991' - MULTIPLE REAL/HALF-COMPLEX PERIODIC
C       FAST FOURIER TRANSFORM
C
C       SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
C       THAT IN MRFFT2
C
C       PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
C       IS GIVEN BY COOLEY, LEWIS ' WELCH (J. SOUND VIB., VOL. 12
C       (1970), 315-337)
C
C       A IS THE ARRAY CONTAINING INPUT ' OUTPUT DATA
C       WORK IS AN AREA OF SIZE (N+1)*LOT
C       TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C       IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
C       INC IS THE INCREMENT WITHIN EACH DATA "VECTOR"
C       (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C       JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C       N IS THE LENGTH OF THE DATA VECTORS
C       LOT IS THE NUMBER OF DATA VECTORS
C       ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C       = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C       ORDERING OF COEFFICIENTS:
C       A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C       WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
C
C       ORDERING OF DATA:
C       X(0),X(1),X(2),...,X(N-1)
C
C       VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
C       PARALLEL
C
C       *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
C
C       DEFINITION OF TRANSFORMS:
C       ----------------------------------------------------------------
C  --
C
C       ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C       WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C
C       ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C       B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C
       SUBROUTINE FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)

       IMPLICIT double PRECISION (A-H,O-Z)

       DIMENSION A(*),WORK(*),TRIGS(*),IFAX(*)
C
       NFAX=IFAX(1)
       NX=N+1
       NH=N/2
       INK=INC+INC
       IF (ISIGN.EQ.+1) GO TO 30
C
C       IF NECESSARY, TRANSFER DATA TO WORK AREA
       IGO=50
       IF (MOD(NFAX,2).EQ.1) GOTO 40
       IBASE=1
       JBASE=1
       DO 20 L=1,LOT
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 10 M=1,N
       WORK(J)=A(I)
       I=I+INC
       J=J+1
10     CONTINUE
       IBASE=IBASE+JUMP
       JBASE=JBASE+NX
20     CONTINUE
C
       IGO=60
       GO TO 40
C
C       PREPROCESSING (ISIGN=+1)
C       ----------------------------------------------------------------
C  --
C
30     CONTINUE
       CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
       IGO=60
C
C       COMPLEX TRANSFORM
C       ----------------------------------------------------------------
C  --
C
40      CONTINUE
       IA=1
       LA=1
       DO 80 K=1,NFAX
       IF (IGO.EQ.60) GO TO 60
50     CONTINUE
       CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,
     *       INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
       IGO=60
       GO TO 70
60     CONTINUE
       CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,
     *       2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
       IGO=50
70     CONTINUE
       LA=LA*IFAX(K+1)
80     CONTINUE
C
       IF (ISIGN.EQ.-1) GO TO 130
C
C       IF NECESSARY, TRANSFER DATA FROM WORK AREA
       IF (MOD(NFAX,2).EQ.1) GO TO 110
       IBASE=1
       JBASE=1
       DO 100 L=1,LOT
       I=IBASE
       J=JBASE
CDIR$ IVDEP
       DO 90 M=1,N
       A(J)=WORK(I)
       I=I+1
       J=J+INC
90     CONTINUE
       IBASE=IBASE+NX
       JBASE=JBASE+JUMP
100    CONTINUE
C
C       FILL IN ZEROS AT END
110    CONTINUE
       IB=N*INC+1
CDIR$ IVDEP
       DO 120 L=1,LOT
       A(IB)=0.D00
       A(IB+INC)=0.D00
       IB=IB+JUMP
120    CONTINUE
       GO TO 140
C
C       POSTPROCESSING (ISIGN=-1):
C       ----------------------------------------------------------------
C  --
C
130     CONTINUE
       CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
C
140    CONTINUE
       RETURN
       END
C
C       SUBROUTINE FFT99AD - PREPROCESSING STEP FOR FFT99, ISIGN=+1
C
C       (SPECTRAL TO GRIDPOINT TRANSFORM)
C
       SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)

       IMPLICIT double PRECISION (A-H,O-Z)

       DIMENSION A(*),WORK(*),TRIGS(*)
       NH=N/2
       NX=N+1
       INK=INC+INC
C
C       A(0) ' A(N/2)
       IA=1
       IB=N*INC+1
       JA=1
       JB=2
CDIR$ IVDEP
       DO 10 L=1,LOT
       WORK(JA)=A(IA)+A(IB)
       WORK(JB)=A(IA)-A(IB)
       IA=IA+JUMP
       IB=IB+JUMP
       JA=JA+NX
       JB=JB+NX
10     CONTINUE
C
C       REMAINING WAVENUMBERS
       IABASE=2*INC+1
       IBBASE=(N-2)*INC+1
       JABASE=3
       JBBASE=N-1
C
       DO 30 K=3,NH,2
       IA=IABASE
       IB=IBBASE
       JA=JABASE
       JB=JBBASE
       C=TRIGS(N+K)
       S=TRIGS(N+K+1)
CDIR$ IVDEP
       DO 20 L=1,LOT
       WORK(JA)=(A(IA)+A(IB))-
     *       (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
       WORK(JB)=(A(IA)+A(IB))+
     *       (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
       WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+
     *       (A(IA+INC)-A(IB+INC))
       WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-
     *       (A(IA+INC)-A(IB+INC))
       IA=IA+JUMP
       IB=IB+JUMP
       JA=JA+NX
       JB=JB+NX
20     CONTINUE
       IABASE=IABASE+INK
       IBBASE=IBBASE-INK
       JABASE=JABASE+2
       JBBASE=JBBASE-2
30     CONTINUE
C
       IF (IABASE.NE.IBBASE) GO TO 50
C       WAVENUMBER N/4 (IF IT EXISTS)
       IA=IABASE
       JA=JABASE
CDIR$ IVDEP
       DO 40 L=1,LOT
       WORK(JA)=2.0*A(IA)
       WORK(JA+1)=-2.D00*A(IA+INC)
       IA=IA+JUMP
       JA=JA+NX
40     CONTINUE
C
50     CONTINUE
       RETURN
       END
C
C       SUBROUTINE FFT99BD - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
C       (GRIDPOINT TO SPECTRAL TRANSFORM)
C
       SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)

		IMPLICIT double PRECISION (A-H,O-Z)

       DIMENSION WORK(*),A(*),TRIGS(*)
C
       NH=N/2
       NX=N+1
       INK=INC+INC
C
C       A(0) ' A(N/2)
       SCALE=1.D00/DFLOAT(N)
       IA=1
       IB=2
       JA=1
       JB=N*INC+1
CDIR$ IVDEP
       DO 10 L=1,LOT
       A(JA)=SCALE*(WORK(IA)+WORK(IB))
       A(JB)=SCALE*(WORK(IA)-WORK(IB))
       A(JA+INC)=0.D00
       A(JB+INC)=0.D00
       IA=IA+NX
       IB=IB+NX
       JA=JA+JUMP
       JB=JB+JUMP
10     CONTINUE
C
C       REMAINING WAVENUMBERS
       SCALE=0.5D00*SCALE
       IABASE=3
       IBBASE=N-1
       JABASE=2*INC+1
       JBBASE=(N-2)*INC+1
C
       DO 30 K=3,NH,2
       IA=IABASE
       IB=IBBASE
       JA=JABASE
       JB=JBBASE
       C=TRIGS(N+K)
       S=TRIGS(N+K+1)
CDIR$ IVDEP
       DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB))
     *       +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB))
     *      -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *       +(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *       -(WORK(IB+1)-WORK(IA+1)))
       IA=IA+NX
       IB=IB+NX
       JA=JA+JUMP
       JB=JB+JUMP
20     CONTINUE
       IABASE=IABASE+2
       IBBASE=IBBASE-2
       JABASE=JABASE+INK
       JBBASE=JBBASE-INK
30     CONTINUE
C
       IF (IABASE.NE.IBBASE) GO TO 50
C       WAVENUMBER N/4 (IF IT EXISTS)
       IA=IABASE
       JA=JABASE
       SCALE=2.D00*SCALE
CDIR$ IVDEP
       DO 40 L=1,LOT
       A(JA)=SCALE*WORK(IA)
       A(JA+INC)=-SCALE*WORK(IA+1)
       IA=IA+NX
       JA=JA+JUMP
40     CONTINUE
C
50     CONTINUE
       RETURN
       END
