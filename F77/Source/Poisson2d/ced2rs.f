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
	SUBROUTINE CED2RS(N,IND,CC,CS,INDEX,F)

		implicit double precision (a-h,o-z)

C
C		SUBROUTINE POUR LE CALCUL DE LA DERIVEE 2ME D' UNE
C		FONCTION AVEC ECHANTILLONAGE RAREFIE A L'ORIGINE.
C
C	PARAMETRE DE LA SUBROUTINE:
C
C		N	= NOMBRE DES DEGRES DE LIBERTE.
C		IND	= PARAMETRE: IND=0 SI LA FONCTION EST PAIRE,
C			  IND=1 SI LA FONCTION EST IMPAIRE.
C		CC	= TABLEAU D'ENTRE' CONTENANT LES COEFFICIENTS
C			  DE TCHEBYTCHEV DE LA FONCTION DONT ON VEUT
C			  CALCULER LA DERIVEE 2ME.
C		CS	= TABLEAU DE SERVICE ET TABLEAU DE SORIE SI 
C			  INDEX=1. DANS CE CAS DANS CS IL Y A LES 
C			  COEFFICIENTS DE TCHEBYTCHEV DE LA DERIVEE 2ME
C			  DE LA FONCTION.
C		INDEX	= PARAMETRE: SI INDEX=1 ON A EN SORTIE DANS CS LES
C			  COEFFICIENTS DE LA DERIVEE, SI INDEX=2 ON A
C			  EN SORTIE DANS F LA DERIVEE 2ME DE LA FON-
C		          CTION. SUBROUTINE DISTRUPTIVE SI INDEX=2.
C
C		Subroutine testee le 24/4/85.
C
C....................................................................
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.1  1997/10/21  14:02:45  eric
c Initial revision
c
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	DIMENSION CC(*),CS(*),F(*)
	N1=N+1
	N11=N+2
	NN11=N1+N1
	N10=N-1
	SOM=0
	SOM1=0
C
C		CALCUL DE LA DERIVEE 2ME D'UNE FONCTION PAIRE.
C
	IF(IND.EQ.0) THEN
	C1=CC(N1)
	CC(N1)=C1*.5
C
	DO 1 LP=1,N
	LP2=LP+LP-1
	L=NN11-LP2
	L1=L-1
	CP=CC(N11-LP)*L1
	SOM=SOM+CP
	SOM1=SOM1+CP*L1**2
	CS(N1-LP)=SOM1-SOM*(L-3)**2
1	CONTINUE
	CS(N1)=0
	CC(N1)=C1
	IF(INDEX.EQ.1) RETURN
	CALL CIRARS(N,0,CS,CC,F)
	RETURN
	ENDIF
C
C		CALCUL DE LA DERIVEE 2ME D'UNE FONCTION IMPAIRE.
C
	IF(IND.EQ.1) THEN
C
	DO 2 LP=1,N
	LP2=LP+LP-2
	L=NN11-LP2
	L1=L-1
	CP=CC(N11-LP)*L1
	SOM=SOM+CP
	SOM1=SOM1+CP*L1**2
	CS(N1-LP)=SOM1-SOM*(L-3)**2
2	CONTINUE
	CS(N1)=0
	IF(INDEX.EQ.1) RETURN
	CALL CIRARS(N,1,CS,CC,F)
	ENDIF
	RETURN
	END
