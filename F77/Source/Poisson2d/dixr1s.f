C
C   Copyright (c) 1998 Silvano Bonazzola
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
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.1  1998/06/22  10:40:59  eric
c Initial revision
c
C
C $Header$
C
C

	SUBROUTINE DIXR1S(N,IND,CS,CC)

	IMPLICIT double PRECISION (A-H,O-Z)

C
C	    ROUTINE MULTIPLIANT OU DIVISANT UNE FONCTION AVEC ECHAN-
C	    TILLONNAGE RAREFIE' A L'ORIGINE PAR x. LA FONCTION AINSI
C	    TRAITEE CHANGE DE PARITE'. (DI(vision) (par) X (d'une
C	    fonction) R(arefiee) S(imple precision).
C
C		ARGUMENTS DE LA ROUTINE:
C
C	    N	=NOMBRE DE DEGRES DE LIBERTE'
C
C	    IND	=PARAMETRE:
C
C		    IND=0  LA FONCTION PAIRE EST DIVISEE PAR x ET
C			   RENDUE IMPAIRE,
C		    IND=1  LA FONCTION IMPAIRE EST DIVISEE PAR x
C			   ET RENDUE PAIRE,
C		    IND=2  LA FONCTION PAIRE EST MULTIPLIEE PAR x
C			   ET RENDUE IMPAIRE,
C		    IND=3  LA FONCTION IMPAIRE EST MULTIPLIEE PAR x
C			   ET RENDUE PAIRE.
C
C	    CS	=TABLEAU INPUT
C	    CC	=TABLEAU OUTPUT
C
C##	Routine testee le 13/04/85
C
C......................................................................
C
	character*100 header
	data header/'$Header$'/

	DIMENSION CS(*),CC(*)
C
C		CAS IND=0: LA FONCTION PAIRE EST DIVISEE PAR x
C
	N1=N+1
C
	IF(IND.EQ.0) THEN
C
	CC(1)=-CS(1)
	DO L=2,N
	CC(L)=-CC(L-1)-CS(L)-CS(L)
	ENDDO
	CC(N1)=0
	RETURN
	ENDIF
C
C		CAS IND=3: LA FONCTION IMPAIRE EST MULTIPLIEE PAR x
C
	IF(IND.EQ.3) THEN
C
	CC(1)=-CS(1)
	DO L=2,N
	CC(L)=-(CS(L)+CS(L-1))*.5D+00
	ENDDO
	CC(N1)=-CS(N)
	RETURN
	ENDIF
C
C		CAS IND=2: LAFONCTION PAIRE EST MULTIPLIEE PAR x
C
	IF(IND.EQ.2) THEN
	DO L=1,N
	CC(L)=-(CS(L)+CS(L+1))*.5D+00
	ENDDO
	CC(N1)=0
	RETURN
	ENDIF
C
C		CAS IND=1: LA FONCTION IMPAIRE EST DIVISEE PAR x
C
	IF(IND.EQ.1) THEN
C
	CC(N)=-CS(N)-CS(N)
	DO L=1,N-1
	LL=N-L
	CC(LL)=-(CS(LL)+CS(LL)+CC(LL+1))
	ENDDO
	CC(N1)=0
	ENDIF
C
C  100	FORMAT(1X,'PIMP',1X,5D24.16)
C  101	FORMAT(1X,' ')
	RETURN
	END
