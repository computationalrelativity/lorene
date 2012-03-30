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

	SUBROUTINE DIRA2S(N,IPAR,MULT,CS,CC)

	IMPLICIT double PRECISION (A-H,O-Z)

C		
C	ROUTINE POUR LA MULTIPLICATION PAR x**2 D' UNE FONCTION AVEC
C	ECHANTILLONNAGE RAREFIE' A L'ORIGINE.(FONCTION PAIRE OU IMPAIRE),
C	(DI(vision) M(ultiple) R2(echantillonnage rarefie') S(imple precision)
C
C	ARGUMENTS DE LA ROUTINE:
C
C		N	=NOMBREDE DEGRES DE LIBERTE-1
C		IPAR	=PARAMETRE: SI LA FONCTION QUI DOIT ETRE DIVISEE
C			 OU MULTPLIEE EST PAIRE IPAR=0, IPAR=1 DANS LE
C			 CAS CONTRAIRE.
C		MULT	=PARAMETRE: SI MULT=1 LA MULTIPLICATION PAR x**2 EST
C			 EFFECTUEE, SI MULT=-1 LA DIVISION EST EFFECTUEE.
C		CS	=TABLEAU CONTENANT LES N+1 COEFFICIENTS DE TCHEBY-
C			 TCHEV DE LA FONCTION QUI DOIT ETRE DIVISEE OU
C			 MULTIPLIEE PAR x**2.
C		CC	=TABLEAU CONTENANT LES COEFFICIENTS DE LA FONCTION
C			 DIVISEE.
C
C	(Routine testee le 28/9/1986)
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
c Revision 1.2  1997/05/23  11:29:24  hyc
c *** empty log message ***
c
C Revision 1.1  1997/05/08 07:33:03  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

	DIMENSION CS(*),CC(*)

	N1 =N+1
	N01=N-1
C
C	CAS FONCTIONS PAIRES
C
	IF (IPAR.EQ.0) THEN
C
C	MULTIPLICATION FONCTIONS PAIRES
C
	IF (MULT.EQ.1) THEN
	CC(1)=(CS(1)+CS(2))*.5
	DO L=2,N
	CC(L)=(CS(L-1)+CS(L+1))*.25+CS(L)*.5
	ENDDO
	CC(N1)=(CS(N)+CS(N1))*.5
	RETURN
	ENDIF
C
C
C	DIVISION FONCTIONS PAIRES
C
	IF (MULT.EQ.-1) THEN
C
C		COMBINAISON LINEARE DES COEFFICIENTS DU 2EME MEMBRE
C		DU SYSTEME
C
	DO L=2,N
	CS(L)=CS(L)*2-CS(L-1)
	ENDDO
	CS(N1)=CS(N1)-CS(N)
C
C		INVERSION DE LA MATRICE
C
	CC(N1)=0
	CC(N)=CS(N)*2
	DO L=1,N01
	LL=N-L
	CC(LL)=2*CS(LL)-CC(LL+1)
	ENDDO
	RETURN
	ENDIF
	ENDIF
C
C		CAS FONCTIONS FONCTIONS IMPAIRES
C
	IF (IPAR.EQ.1) THEN
C
C		MULTIPLICATION FONCTIONS IMPAIRES
C
	IF (MULT.EQ.1) THEN
C
	CC(1)=.75*CS(1)+.25*CS(2)
	DO L=2,N01
	CC(L)=(CS(L-1)+CS(L+1))*.25+CS(L)*.5
	ENDDO
	CC(N)=CS(N01)*.25+CS(N)*.5
	CC(N1)=0
	RETURN
	ENDIF
C
C
C		DIVISION FONCTIONS IMPAIRES
C
	IF (MULT.EQ.-1)THEN
C
C		COMBINAISON LINEAIRE DES COEFFICIENTS DU 2EME MEMBRE
C		DU SYSTEME
C
	DO L=2,N
	L1=L-1
	CS(L)=CS(L)*(L+L1)-CS(L1)
	ENDDO
C
C		SOLUTION DU SYSTEME
C
	CC(N)=CS(N)/N
	ALN=N*.5-.25
	DO L=1,N01
	AL1=ALN-.5
	LL=N-L
	CC(LL)=(CS(LL)-CC(LL+1)*AL1)/ALN
	ALN=ALN-.5
	ENDDO
C
	ENDIF
	ENDIF
	RETURN
	END
