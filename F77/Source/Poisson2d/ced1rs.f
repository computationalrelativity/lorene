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
       SUBROUTINE CED1RS(N,IND,CC,CS,INDEX,F)

		implicit double precision (a-h,o-z)

C
C***** CALCUL DE LA DERIVEE PREMIERE D'UNE FONCTION A ECHANTILLONAGE 
C      RAREFIE A L'ORIGINE.
C           
C           ARGUMENTS DE LA SSUBROUTINE:
C           
C           N   = NOMBRE DES DEGRES DE LIBERTEE.
C           IND = PARMETRE: SI IND=0 ON TRAITE LES FONCTIONS
C                 PAIRES, SI IND=1 LES FONCTION IMPAIRES, SI
C                 IND=2 ON CALCULE LES COEFFICIENTS DE TCHE-
C                  BITCHEF DE LA FONCTION 1/r*D/dr DES FONCTIONS 
C                  PAIRES
C           CC  = TABLEAU D'ENTREE: IL CONTIENT LES COEFF. DE
C                 TCHEBYTCHEV DE LA FONCTION D'ONT ON VEUT CAL-
C                 CULER LA DERIVEE.
C           CS  = TABLEAU DE TRAVAIL ET TABLEAU OUTPUT SI
C                 INDEX=1.
C           INDEX = PARMETRE: SI INDEX=1 ON A EN SOTIE DANS CS
C                 LES COEFFICIENTS DE LA DERIVEE DE LA FONCTION
C               SI INDEX=2 ON EN SORTIE DANS Y LA FONCTION DE-
C                RIVEE.
C           Y   = TABLEAU OUTPUT
C
C****    SUBROUTINE DISTRUPTIVE.
C
C           SUBROUTINE TESTEE LE 22/4/85.
C
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.1  1997/10/21  14:02:07  eric
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
       SOM=0
C
C           CALCUL DE 1/R*D/DR D'UNE FONCTION PAIRE.
C
       IF(IND.EQ.2) THEN
       C1=CC(N1)
       DO 4 L=1,N1
       CS(L)=0
4      CONTINUE
       CC(N1)=C1*.5
       SOM=0
       DO 1 LP=1,N,2
       L=N11-LP
       L1=L-1
       SOM=SOM  -  (L1+L1)*CC(L)
       CS(L1)=-SOM*4
1      CONTINUE
       SOM=0
       DO 2 LP=2,N,2
       L=N11-LP
       L1=L-1
       SOM=SOM  -  (L1+L1)*CC(L)
       CS(L1)=-4*SOM
2      CONTINUE
       CC(N1)=C1
       CS(N1)=0
       IF(INDEX.EQ.1) RETURN
       CALL CIRARS(N,0,CS,CC,F)
       RETURN
       ENDIF
C
C           CALCUL DE LA DERIVEE D'UNE FONCTION PAIRE.
C      
C
       IF(IND.EQ.0) THEN
C
       C1=CC(N1)
       CC(N1)=CC(N1)*.5
       DO 100 LP=1,N
       L=N1-LP
       SOM=SOM-CC(L+1)*(4*L)
       CS(L)=SOM
100	CONTINUE
c
       CC(N1)=C1
       CS(N1)=0
       IF(INDEX.EQ.1) RETURN
       CALL CIRARS(N,1,CS,CC,F)
       RETURN
       ENDIF
C
C           CALCUL DE LA DERIVEE D'UNE FONCTION IMPAIRE.
C
C
       IF(IND.EQ.1) THEN
C
       SOM=0
       DO 123 LP=1,N
       L=N1-LP
       SOM=SOM-CC(L)*(4*L-2)
       CS(L)=SOM
       CS(N1)=0
123	CONTINUE
       IF(INDEX.EQ.1) RETURN
       CALL CIRARS(N,0,CS,CC,F)
       ENDIF
       RETURN
       END
