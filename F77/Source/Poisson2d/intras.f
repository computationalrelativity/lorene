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

       SUBROUTINE INTRAS(N,CC,OUT)

		implicit double precision (a-h,o-z)

C
C ***** CALCUL DE L'INTEGRALE ENTRE 0 ET 1 D'UNE FONCTION ECHENTILLONEE
C ***** RAREMENT. (FONCTION PAIRE PAR RAPPORT TETA=PI/2)
C
C           ARGUMENTS DE LA ROUTINE:
C
C      N=   NOMBRE DES DEGRES DE LIBERTEE-1
C      CC   = COEFFICIENTS DE TCHEBYTCHEV DE LA FONCTION A INTEGRER.
C      OUT  =VALEURE DE L'INTEGRALE
C
C           Routine testee avec tous les tests de Routine
C           (test de repetition compris le 4/4/1986)
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
c Revision 1.2  1997/05/23  11:40:14  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:40:48  hyc
C Initial revision
C
C
C $Header$
C
C
	character*120 header
	data header/'$Header$'/

       DIMENSION CC(*)
       N1=N+1
       C1=CC(N1)
        CC(N1)=C1*.5
       SOM=0
       COS2=-1
       X1=1
       DO 5 L=2,N1
       COSL=-X1/(L+L-1)
       SOM=SOM+CC(L)*(COS2-COSL)
       COS2=COSL
  5    CONTINUE
       OUT=.5*(SOM+CC(1))
       CC(N1)=C1
  100  FORMAT(1X,10D24.16)
       RETURN
       END
