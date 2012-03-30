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

       SUBROUTINE DIMRAS(N,IPAR,MULT,N64,CS,CC)

		implicit double precision (a-h,o-z)

C
C      
C           ROUTINE POUR LA MULTIPLICATION PAR x D' UNE FONCTION
C           AVEC ECHANTLLONAGE RAREFIE' A L'ORIGINE.(FONCTION PAIRE OU
C           IMPAIRE),(DI(vision) M(ultiple) RA(echantillonage rarefie') 
C           S(imple precision).
C
C           ARGUMENTS DE LA ROUTINE:
C
C           N   =NOMBREDE DEGRES DE LIBERTE-1
C           IPAR =PARAMETRE: SI LA FONCTION QUI DOIT ETRE DIVISEE
C                OU MULTPLIEE EST PAIIRE IPAR=0, IPAR=1 DANS LE
C                CAS CONTRAIRE.
C           MULT =PARAMETRE: SI MULT=1 LA MULTIPLICATION PAR x  EST
C                EFFECTUEE SI MULT=-1 LA DIVISION EST EFFECTUEE.
C           N64 =N64 EST LE NOMBRE DE FONCTIONS QUI DOIVENT ETRE
C                DIVISEES (OU MULTIPLIEES) PAR x
C           CS  =TABLEAU CONTENANT LE N+1 COEFFICIENTS DE TCHEBY-
C                TCHEV DE LA FONCTION QUI DOIT ETRE DIVISEE OU
C                MULTIPLIEE PAR x.
C           CC  =TABLEAU CONTENANT LES COEFF. DE LA FONCTION DIVISEE.
C                LES DIMENSIONS DES TABLEAUX CC ET CS DOIBENT
C                ETRE > (N64+1)*N1
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
c Revision 1.2  1997/05/23  11:29:55  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:34:56  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

       DIMENSION CS(*),CC(*)
       N1=N+1
       N0=N-1
       N65=N64
       IF((N64/8)*8.EQ.N64) N65=N64+1
       N63=N64-1
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
       DO 1 L=1,N
       DO 2 M=JM1,JM2
       CC(M)=-(CS(M)+CS(M+N65))*.5
   2   CONTINUE
       JM1=JM1+N65
       JM2=JM1+N63
   1   CONTINUE
C
C           LE DERNIER CORFFICIENT ESTPOSE'=0.
C
       JM1=N*N65+1
       JM2=JM1+N64
       DO 3 M=JM1,JM2
       CC(M)=0
   3   CONTINUE
       RETURN
       ENDIF
C
C           DIVISION FONCTIONS PAIRES
C
       IF(MULT.EQ.-1) THEN
       DO 4 M=1,N64
       CC(M)=-CS(M)
   4   CONTINUE 
       JM1=N65+1
       JM2=JM1+N63
       DO 5 L=2,N
       DO 6 M=JM1,JM2
       CS(M)=-(CS(M)*2+CC(M-N65))
   6   CONTINUE
       DO 7 M=JM1,JM2
       CC(M)=CS(M)
   7   CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
   5   CONTINUE
C
       JM1=N*N65+1
       JM2=JM1+N63
       DO 8 M=JM1,JM2 
       CC(M)=0
   8   CONTINUE
       ENDIF
       RETURN
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
       CC(M)=-CS(M)
   9   CONTINUE
C
       JM1=N65+1
       JM2=N63+JM1
       DO 10 L=2,N
       DO 11 M=JM1,JM2
       CC(M)=-.5*(CS(M-N65)+CS(M))
  11   CONTINUE
       JM1=JM1+N65
       JM2=JM1+N63
  10   CONTINUE
C
       JM1=N65*N+1
       JM2=JM1+N63
       DO 12 M=JM1,JM2
       CC(M)=-CS(M-N65)
  12   CONTINUE
       RETURN
       ENDIF
C
C           DIVISION FONCTIONS IMPAIRES
C
       IF(MULT.EQ.-1)THEN
       JM1=N65*N+1
       JM2=JM1+N63
       DO 13 L=JM1,JM2
       CC(L)=0
  13   CONTINUE
C
       JM1=N0*N65+1
       JM2=JM1+N63
       DO 14 M=JM1,JM2
       CC(M)=-2*CS(M)
  14   CONTINUE
C
       JM1=(N-2)*N65+1
       JM2=JM1+N63
       DO 15 L=1,N0
       DO 16 M=JM1,JM2
       CS(M)=-(CC(M+N65)+2.*CS(M))
  16   CONTINUE
       DO 17 M=JM1,JM2
       CC(M)=CS(M)
  17   CONTINUE
       JM1=JM1-N65
       JM2=JM1+N63
  15   CONTINUE
C
       ENDIF
       ENDIF
       RETURN
       END
       
C
C
