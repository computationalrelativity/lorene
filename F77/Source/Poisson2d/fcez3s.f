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

	subroutine fcez3s(ndl,ndr,ndt,ndf,n64,in1,imp,c64,cc,cs,
     1	dent,den,denn)
c

	implicit double precision(a-h,o-z)

c		subroutine pour le calcul des transformees de fourier
c		e de tchebytchev pour l'echantillonage rarefie' a 3 dim.
c		dans le cas multizones
c		exploitant le simmetries des fonctions a transformer (sym-
c		metries existantes par ex. en coordonnes spheriques)
c	
c		le 3me indice est suppose' etre periodiques (coordonne azi-
c		muthale fi), le 2me varie entre 0 et pi. dans un developpement
c		en coordonnes spheriques, les coefficients de fourier d'une
c		fonction scalaire seront de fonctions symetriques en teta si
c		les coefficients  de fourier sont paires et antisymetriques
c		sont impairs. le developpement en teta est donc effectue' en  
c		polynomes de tchebytchev du 1er genre pour les fonctions
c		paires et du 2me genre (en serie de sinus) pour les fonctions
c		impaires. analoguement les coefficients cml(r) sont des
c		fonctions symetriques en r si m+l est paire et antisymetriques
c		dans le cas oppose'. parconsequant le developpement en r
c		des coefficients  cml(r) est effectue sur l'intervalle
c		0<r<1 en tenant compte de la parite'. la transformation 
c	n.b.	doit etre parconsequent ordonnee, c'est a dire il faut
c	---	d'abord proceder a la transformation de fourier sur
c		la variable fi (3me indice), puis a la transformation en teta 
c		(2me indice) et enfin la transformation en r (pemier indice).
c		cela peut etre genant. la subroutine fger3s evite cet inco-
c		venient.
c
c		le stockage des coefficients est le suivant (cfr.
c		la subroutine fuce3s). dans den(lr,lt,1) il y a le coefficient
c		correspondants a la frequence zero du developpement en cosinus
c		dans den(lr,lt,2),den(lr,lt,3) les cofficients cosinus et
c		sinus de la frequence 1, dans den(lr,lt,4), den(lr,lt,5)
c		les coefficients de la frequence 2 et ainsi de suite.
c			 den(lr,1,lm),den(lr,3,lm).... den(lr,2*n+1,lm)  
c		sont les coefficients du developpement sur les polynome
c		de tchebytchev du 1er ordre. den(lr,2,lm),den(lr,4,lm)...
c		den(lr,2*n,lm) les coefficients de tchebytchev du 2m ordre.
c		den(1,lr,lm),den(3,lr,lm).... sont les coefficients
c		du developpement en polynomes de chebytchev des fon-
c		ctions symetriques en r, et den(2,lr,lm),den(4,lr,lm)....
c		des fonctions antisymetriques.	
c		
c		subroutine completement craytinizee.
c
c		subroutine ayant teste'e avec le protocol usuel le jour
c		du segnuer 4/2/1987.
c
c		arguments de la subroutine:
c
c		ndl	= tableau, ndeg(3) contenant les de-
c			  grees de liberte des differents zones ou co-
c			  quilles.
c			  ndl(1) contient le nombre nzon des coquilles
c			  ndl(2),ndl(3),...ndl(nzon+1) le nombre des de-
c			  gres de liberte en r de la 1ere,2me,...nzon-eme
c			  coquille, ndl(nzon+2),ndl(nzon+3), les degres
c			  de lberte en thete et en phi.
c	ndr,ndt,ndfd	= dimensions des differents ta-
c			bleaux comme declare dans le programme appellant.
c		pour des raisons de craytinisation nndr,ndt,ndf ne 
c		doivent pas etre un multiple de 8.
c
c		n64	= parametre de la vectorization, par exemple
c			 n64=64 signifie que 64 fonctions a transformer
c			 sont vectorizee.
c
c		in1	= parametre, si in1=1 la transformee est
c			 effectuee sur le premier indice, si in1=2 sur le
c			 deuxieme, si in1=3 sur le 3me,.
c		imp	= parametre indicant la tenserioalite' de l'objet a
c			 transformer,(imp=0 tenseurs d'ordre 0,2,4,...,imp=1
c			 tenseurs d'ordre 1,3,5...), dans le cas ou il y
c			 aurait une symetrie par rapport le plan equato-riale
c			 ou une super symetrie invariance de la fonction
c			 parrapport la transformation x,y -> -x,-y voir
c			 les valeurs a donner a imp dans la routine fcer3s.
c
c		c64,cc,cs= tableaux de travail: dimension minime=
c			   (n64+1)*((max(ndeg(1),ndeg(2))+3)
c
c		dent	=tableau de travil a 3 dimensions. dimensions
c			 minimales de dent max(ndl(2),ndl(3),...ndl(nzon+1))
c			 ,nl(nzon+2),ndl(nzon+3)/2+1
c		den	=tableau de travil a 3 dimensions. dimensions
c			 minimales  max(ndl(2),ndl(3),...ndl(nzon+1))
c			 ,nl(nzon+2),ndl(nzon+3)
c		denn	=tableau a 4 dimensions contenant la fonction
c			 a transformer en imput, et la transformee en
c			 output. 
c
c	routine modifiee le 28/octobre 1994 - den(ndr,ndt,*),denn(ndr,ndt,ndf,*),
c
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.1  1997/10/23  08:07:19  eric
c Initial revision
c
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	dimension ndl(*)
	dimension ndeg(3),den(ndr,ndt,*),c64(*),cc(*),cs(*)
	dimension dent(ndr,ndt,*),denn(ndr,ndt,ndf,*)
c
	nzon=ndl(1)
	ny1=ndl(nzon+2)
	nf= ndl(nzon+3)
c
	ndeg(2)=ny1
	ndeg(3)=nf
c
	if(in1.eq.1) then
c
	do 10 lzon=1,nzon
	nr1=ndl(lzon+1)
	ndeg(1)=nr1
c
	do 1 lf=1,nf
	do 2 ly=1,ny1
	do 3 lr=1,nr1
	den(lr,ly,lf)=denn(lr,ly,lf,lzon)
   3	continue
   2	continue
   1	continue
c
	if(lzon.eq.1) then
	call fcer3s(ndeg,ndr,ndt,n64,in1,imp,c64,cc,cs,dent,den)
	else
	call fuce3s(ndeg,ndr,ndt,n64,2,1,c64,cc,cs,den)
	endif
c
	do 4 lf=1,nf
	do 5 ly=1,ny1
	do 6 lr=1,nr1
	denn(lr,ly,lf,lzon)=den(lr,ly,lf)
   6	continue
   5	continue
   4	continue
  10	continue
	return
	endif
c
	if(in1.gt.1) then
c
	do 20 lzon=1,nzon
	nr1=ndl(lzon+1)
	ndeg(1)=nr1
c
	do 11 lf=1,nf
	do 12 ly=1,ny1
	do 13 lr=1,nr1
	den(lr,ly,lf)=denn(lr,ly,lf,lzon)
   13	continue
   12	continue
   11	continue
c
	call fcer3s(ndeg,ndr,ndt,n64,in1,imp,c64,cc,cs,dent,den)
c
	do 14 lf=1,nf
	do 15 ly=1,ny1
	do 16 lr=1,nr1
	denn(lr,ly,lf,lzon)=den(lr,ly,lf)
   16	continue
   15	continue
   14	continue
   20	continue
	endif
	return
	end

