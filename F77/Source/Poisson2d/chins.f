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

	SUBROUTINE CHINS(N,F,CC,CS)

	IMPLICIT double PRECISION (A-H,O-Z)

C
C##	version du 19/11/93: derniere dimension des tableaux avec *
C
C
C **** INVERSION DE LA TF DE TCHEBYCHEV, F EST LE TABLEAU CONTENANT LE
C **** DE TCHEBYCHEV, CS EST L'OUTPUT. N DOIT ETRE UNE PUISSANCE DE 2.
C **** DIMENSION MINIMALE DES TABLEAUX=N+1
C **** ROUTINE DISTRUPTIVE
C
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:32:54  hyc
c *** empty log message ***
c
C Revision 1.1  1997/05/08 07:37:44  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	DIMENSION CC(*),CS(*),F(*)
	N1=N+1
	N2=N/2
	N21=N2+1
	N3=N2-1
	X0=0
	PI=2.*ACOS(X0)/N
C
C
C **** CALCUL DE LA FONCTION EN THETA=0 ET THETA=PI
C
	SOMM=0
	DO 4 L=3,N,2
	SOMM=SOMM+F(L)
4	CONTINUE
	SOM1=0
	DO 5 L=2,N ,2
	SOM1=SOM1+F(L)
5	CONTINUE
	FC=   (F(1)+F(N1))*.5
	F11=SOMM+SOM1    +FC
	FN21=SOMM-SOM1    +FC
C
C **** CALCUL DES COEFFICIENTS DE FOURIER DE LA FONCTION PONDEREE A PARTIR
C **** DES COEFFICIENTS DE TCHEBYCHEV.
C
	DO 1 L=1,N21
	CS(L)=0
1	CONTINUE
	DO 2 L=1,N2
	FL=F(L+L)
	CS(L)=CS(L)+FL
	CS(L+1)=CS(L+1)-FL
	SOMM=SOMM+FL
	CC(L)=F(L+L-1)
2	CONTINUE
	CC(N21)=F(N1)
	CS(1)=0
	CS(N21)=0
C
C *** INVERSION POUR OBTENIR LA FONCTION PONDEREE
C
	CALL TFINS(N,CC,CS,F)
C
C **** RESTITUTION DE LA FONCTION
C
	DO 3 L=1,N3
	N21L=N21+L
	N20L=N21-L
	F12=(F(N20L)-F(N21L))*(.5+.25/SIN((N2-L)*PI))
	CS(N21L)=F(N21L)+F12
	CS(N20L)=F(N20L)-F12
3	CONTINUE
	CS(1)=F11
	CS(N21)=F(N21)
	CS(N1)= FN21
100	FORMAT(10X,10D12.4)
101	FORMAT(1X,'IN')
	RETURN
	END
