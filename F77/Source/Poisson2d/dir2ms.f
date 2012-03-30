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

       SUBROUTINE DIR2MS(N,IPAR,MULT,N64,CC,CS,Y)

		implicit double precision (a-h,o-z)

C
C      
C           ROUTINE POUR LA MULTIPLICATION PAR x**2 D' UNE FONCTION
C           AVEC ECHANTLLONAGE RAREFIE' A L'ORIGINE.(FONCTION PAIRE OU
C           IMPAIRE),(DI(vision) R(echantillonage rarefie')2 (X**2)
C           M(multiple) S(simple precision) 
C
C           ARGUMENTS DE LA ROUTINE:
C
C           N   =NOMBREDE DEGRES DE LIBERTE-1
C           IPAR =PARAMETRE: SI LA FONCTION QUI DOIT ETRE DIVISEE
C                OU MULTPLIEE EST PAIIRE IPAR=0, IPAR=1 DANS LE
C                CAS CONTRAIRE.
C           MULT =PARAMETRE: SI MULT=1 LA MULTIPLICATION PAR x**2 EST
C                EFFECTUEE SI MULT=-1 LA DIVISION EST EFFECTUEE.
C           N64 =N64 EST LE NOMBRE DE FONCTIONS QUI DOIVENT ETRE
C                DIVISEES (OU MULTIPLIEES) PAR x**2
C           CS  =TABLEAU CONTENANT LE N+1 COEFFICIENTS DE TCHEBY-
C                TCHEV DE LA FONCTION QUI DOIT ETRE DIVISEE OU
C                MULTIPLIEE PAR x**2.
C           CC  = TABLEAU DE TRAVAIL
C           Y   =TABLEAU CONTENANT LES COEFF. DE LA FONCTION DIVISEE.
C                LES DIMENSIONS DES TABLEAUX CC ET CS DOIVENT
C                ETRE > (N64+1)*N1
C
C           Routine ayant subi tous le test les 8/3/87
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
c Revision 1.3  1997/05/23  11:51:01  hyc
c *** empty log message ***
c
C Revision 1.2  1997/05/23 11:29:34  hyc
C *** empty log message ***
C
C Revision 1.1  1997/03/17 20:34:53  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

      DIMENSION CS(*),CC(*),Y(*)
       DATA NDIM,NEQ/0,0/

	save	N0,N2,N3,N4,N5,N65,N63,N66,N466,N366,N266,N166
	save	NM66,N6565,N6566,N1301,A75,A25,ANM,ANM1,AB1
	save	NDIM,NEQ

C
C           INITIALISATION
C
       IF(NDIM.NE.N) THEN
       N0=N-1
       N2=N-2
       N3=N-3
       N4=N-4
       N5=N-5
       ENDIF
C
       IF(NDIM.EQ.N.AND.NEQ.EQ.N64) GO TO 800
       N65=N64
       IF((N64/8)*8.EQ.N64) N65=N64+1
       N63=N64-1
       N66=N65+1
       N466=N4*N65+1
       N366=N466+N65
       N266=N366+N65
       N166=N266+N65
       NM66=N166+N65
       N6565=N65+N65
       N6566=N6565+1
       N1301=N6565+N6566
       A75=.75
       A25=-.25
       ANM=2./N
       ANM1=N*.125
       AB1=N*.25-.25 
       NDIM=N
       NEQ=N64
  800  CONTINUE
C
C           CAS FONCTIONS PAIRES.
C
       IF(IPAR.EQ.0) THEN
C
C           MULTIPLICATION FONCTIONS PAIRES
C
       IF(MULT.EQ.1) THEN
       JM1=1
       JM2=N64
       DO 3 M=JM1,JM2
       Y(M)=(CC(M)+CC(M+N65))*.5
   3   CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
       DO 1 L=3,N+1
       DO 2 M=JM1,JM2
       Y(M)=(CC(M-N65)+CC(M+N65))*.25+CC(M)*.5
   2   CONTINUE
       JM1=JM1+N65
       JM2=JM1+N63
   1   CONTINUE
C
       DO 14 M=JM1,JM2
       Y(M)=CC(M)+CC(M-N65)*.5
  14   CONTINUE
       RETURN
       ENDIF
C
C           DIVISION FONCTION PAIRES
C
       IF(MULT.EQ.-1) THEN
C
C           COMBINAISON DES COEFFICIENTS DU 2ME MEMBRE DU SYSTEME (FONCT
C  ION
C           PAIRES)
C
       JM1=N66
       JM2=N66+N63
       DO 18 L=1,N
       DO 19 M=JM1,JM2
       CS(M)=CC(M)*2-CC(M-N65)
  19   CONTINUE
       DO 20 M=JM1,JM2
       CC(M)=CS(M)
  20   CONTINUE
       JM1=JM1+N65
       JM2=JM1+N63
  18   CONTINUE
C
C           SOLUTION DU SYSTEME: PARTIE PAIRE
C
       DO 21 M=NM66,NM66+N63
       Y(M)=0
  21   CONTINUE
C
       JM1=N166
       JM2=N166+N63
       DO 4 M=JM1,JM2
       Y(M)=CC(M)*2
   4   CONTINUE
C
       JM1=JM1-N65
       JM2=JM1+N63
       DO 5 L=1,N0
       DO 6 M=JM1,JM2
       CS(M)=CC(M)*2-Y(M+N65)
   6   CONTINUE
       DO 7 M=JM1,JM2
       Y(M)=CS(M)
   7   CONTINUE
       JM1=JM1-N65
       JM2=JM1+N63
   5   CONTINUE
C
       RETURN
       ENDIF
             ENDIF
C
C           CAS FONCTIONS FONCTIONS IMPAIRES
C
       IF(IPAR.EQ.1) THEN
C
C           MULTIPLICATION FONCTIONS IMPAIRES
C
       IF(MULT.EQ.1) THEN
       DO 9 M=1,N64
       Y(M)=CC(M)*.75+CC(M+N65)*.25
   9   CONTINUE
C
       JM1=N65+1
       JM2=N63+JM1
       DO 10 L=2,N0
       DO 11 M=JM1,JM2
       Y(M)=(CC(M-N65)+CC(M+N65))*.25+CC(M)*.5
  11   CONTINUE
       JM1=JM1+N65
       JM2=JM1+N63
  10   CONTINUE
C
       DO 12 M=JM1,JM2
       Y(M)=(CC(M)+CC(M+N65))*.75+CC(M-N65)*.25
  12   CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
       DO 22 M=JM1,JM2
       Y(M)=0
  22   CONTINUE
       RETURN
       ENDIF
C
C           DIVISION FONCTIONS IMPAIRES
C
       IF(MULT.EQ.-1)THEN
C
C           COMBINAISON LINEAIRE DES COEFFICIENTES DU 2ME MEMBRE DU
C           SYSTEME
C
       JM1=N66
       JM2=JM1+N63
       DO 23 L=1,N0
       AL=L+L+1
       DO 24 M=JM1,JM2
       CS(M)=CC(M)*AL-CC(M-N65)
  24   CONTINUE
       DO 25 M=JM1,JM2
       CC(M)=CS(M)
  25   CONTINUE
       JM1=JM1+N65
       JM2=JM1+N63
  23   CONTINUE
C
       JM1=N166
       JM2=JM1+N6565-1
C
       DO 13 M=JM1,JM2
       Y(M)=0
  13   CONTINUE
C
C           SOL;UTION DU SYSTEME: FONCTION IMPAIRES
C
       AL1=(N-2)*.5+.75
       AL2=AL1-.5
       JM1=JM1-N65
       JM2=JM1+N63
       C1=1./AL1
       DO 26  M=JM1,JM2
       Y(M)=CC(M)*C1
  26   CONTINUE
C
       JM1=JM1-N65
       JM2=JM1+N63
       AL1=AL2
       AL2=AL2-.5
       DO 15 L=1,N2
       C1=1./AL1
       DO 16 M=JM1,JM2
       CS(M)=(CC(M)-Y(M+N65)*AL2)*C1
  16   CONTINUE
       DO 17 M=JM1,JM2
       Y(M)=CS(M)
  17   CONTINUE
       JM1=JM1-N65
       JM2=JM1+N63
       AL1=AL2
       AL2=AL2-.5
  15   CONTINUE
C
       ENDIF
       ENDIF
       RETURN
  100  FORMAT(1X,10E12.4)
  101  FOR MAT(1X,' ')
       END
C
