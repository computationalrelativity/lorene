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

	subroutine fuce3s(ndeg,ndimr,ndimy,nnn64,itch,in1,c64,cc,cs,den)

	implicit none
c
c## routine modifiee le 31.03.1994: ajout des save
c
c	implicit double precision(a-h,o-z)
c
c		routine pour le calcul des transformees de fourier
c		e de tchebytchev d'un tableau a 3 dimensions.
c
c		routine completement craytinizee.
c
c		arguments de la routine:
c
c			ndeg	= tableau, ndeg(3) contenant les de-
c				  grees de liberte des transformees
c				  a effectuer, ndeg(1) concerne le pre-
c				  mier indice de la matrice, ndeg(2)
c				  le 2me indice, ndeg(3) le 3me indice de la
c				  matrice.					
c				  ndeg doit imperativement etre de la
c				  forme 2**m*3**p*5**q pour les trans-
c				  formeees de fourier (m,p,q nombres 
c				  intiers)
c				  et 2**p*3**p*5*q+1 pour les transfor-
c				  mees de tchebytchev.
c
c		ndimy,ndimz  =dimesion du tableau yy(lr,ly,lz).
c		pour des raisons de craytinisation ndimy et ndimz ne doit pas
c		etre un multiple de 8.
c
c		nnn64	= parametre de la vectorization, par exemple
c			 nnn64=64 signifie que 64 fonctions a transformer
c			 sont vectorizee.
c
c		itch	=parametre, sii itch.eq.1 la transformee
c			 de fourier est effectuee, si itch=2 la rou-
c			 tine effectue la transformee de tchebytchev.
c
c		in1	= parametre, si in1=1 la transformee est
c			 effectuee sur le premier indice, si in1=2 sur le
c			 deuxieme, si in1=3 sur le 3me,.
c			 par exemple itch=2, in1=1 signifie
c			 qu'on effectue la transformation de tcheby-
c			 tchev sr le premier indice du tableau yy.
c			 si ndeg(in1) < 4 la trasformation n'est pas
c			 effectuee et en output on a le tableau imput
c
c		c64,cc,cs= tableaux de travail: dimension minime=
c			   (nnn64+1)*((max(ndeg(1),ndeg(2))+3)
c
c		den	=tableau a 3 dimensions contenant la fonction
c			 a transformer en imput, et la transformee en
c			 output. les coefficients de fourier (dans le cas
c			 d'une transformation de fourier) sont stockes
c			 dans la facon suivante: (exemple dans le cas
c			 de la transformation du 1er indice) dans den(1,l,m)
c			 il-y-a le coefficient de fourier cosinus de la 
c			 frequence zero, dans les coefficients paires 
c			 (den(2,l,m),den(4,l,m)...den(2*n,l,m) les coefi-
c			 cients en cosinus, dans les coefficients impaires
c			 den(3,l,m),den(5,l,m),...den(2*n-1) les termes en 
c			 sinus. parconsequent le termes den(2,l,m) et den(3,l,m)
c			 sont les termes du developpement de fourier correspon-
c			 dants a la meme frequence. en totale il-y-a 2*n
c			 degres de liberte.
c
c		routine ayant testee avec le protocol ordinaire le 10/12/1986
c		routine modifiee le 28/octobre 1994 - den(ndimr,ndimy,*)
c
c
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.1  1997/10/23  08:12:51  eric
c Initial revision
c
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	REAL*8 den,c64,cc,cs
c
	integer ndeg,ndy,ndimr,ndimy,ndr,ndz,neq,nr1,ny1,nz1,nr,ny,nr1z1,
     1	nz,nnn64,nnn65,nn64r,nn65r,ny1z1,multr,irestr,n6r,mult1,ireee2,
     1	ii,i2,lrr1,lrr2,ldr,n64r,n65r,n6565r,n365r,irr2r,nn64y,nn65y
c     
	integer multy,iresty,n6y,lzz1,lzz2,ldz,n64y,n65y,n6565y,n365y,irr2y,
     1	nn64z,nn65z,nr1y1,multz,irestz,n6z,n64z,n65z,lyy1,lyy2,ldy,irr2z,
     1	nr64r,nr65r,nr65z,ly3
c
     	integer n6565z,n365z,lr1,in1,lr2,lmu,jr,lr,jz,lz,ly,itch,n265,
     1	mu2,lz1,ly1,lrr,lm,lz2,ly2,jjz,lsy,lsr,lz3,jy,nr64y,n63r,
     1	nr65y,nr64z,jjy,n63y
c
	dimension den(ndimr,ndimy,*),c64(*),cc(*),cs(*),ndeg(3)
	data ndy,ndr,ndz/0,0,0/
	data neq/0/

	save ndr, ndy, ndz, neq, nr, ny, nz, nnn65, nn64r, nn65r, ny1z1
	save multr, irestr, n6r, mult1, ireee2, ii, i2, lrr1, lrr2, ldr
	save n64r, n65r, n6565r, n365r, irr2r, nn64y, nn65y, nr1z1, multy 
	save iresty, n6y, lzz1, lzz2, ldz, n64y, n65y, n6565y, n365y, irr2y
	save nn64z, nn65z, nr1y1, multz, irestz, n6z, n64z, n65z, lyy1
	save lyy2, ldy, irr2z, n6565z, n365z

	nr1=ndeg(1)
	ny1=ndeg(2)
	nz1=ndeg(3)
c
c		initialisation.
c
	if(ndy.eq.ny1.and.ndz.eq.nz1.and.ndr.eq.nr1.and.neq.eq.nnn64) go to 800
	ndr=nr1
	ndy=ny1
	ndz=nz1
	neq=nnn64
	nr=nr1-1
	nz=nz1-1
	ny=ny1-1
	nnn65=nnn64
	if((nnn64/8)*8.eq.nnn64)nnn65=nnn64+1
c
c			preparation des quantites ncessaires pour la transfor-
c		mation du 1er index du tableau.
c
	nn64r=nnn64
	nn65r=nnn65
	ny1z1=ny1*nz1
	multr=(ny1z1)/nn64r
	irestr=ny1z1-multr*nn64r
c
c		optimisation de nn64r: on cherche une valeur de nn64r qui
c		soit< nnn64 mais un multiple de ny1. cela en vue de reduire
c		le nombre d'operation dans le transfert des valeurs de den
c		dans c64.
c
	if(multr.gt.0.and.nnn64.gt.ny1) then
	n6r=(nnn64/ny1)*ny1
	mult1=ny1z1/n6r
	ireee2=ny1z1-mult1*n6r
	ii=0
	i2=0
	if(irestr.gt.0)ii=1
	if(ireee2.gt.0)i2=1
	if(mult1+i2.le.multr+ii) then
	multr=mult1
	nn64r=n6r
	nn65r=nn64r
	if((nn64r/8)*8.eq.nn64r)nn65r=nn64r+1
	irestr=ireee2
	endif
	endif

	if(multr.eq.0) then
	nn64r=ny1z1
	nn65r=nn64r
	if((nn64r/8)*8.eq.nn64r)nn65r=nn64r+1
	multr=1
	irestr=0
	endif
	lrr1=1
	lrr2=nn64r/ny1
	ldr=lrr2-lrr1
	if(irestr.gt.0) then
	n64r=irestr
	n65r=n64r
	if((n64r/8)*8.eq.n64r)n65r=n64r+1	
	endif
c
	n6565r=nn65r+nn65r
	n365r=n6565r+nn65r
	irr2r=nn64r-(nn64r/ny1)*ny1
c
c
c			preparation des quantites ncessaires pour la transfor-
c		mation du 2me index du tableau.
c
	nn64y=nnn64
	nn65y=nnn65
	nr1z1=nr1*nz1
	multy=(nr1z1)/nn64y
	iresty=nr1z1-multy*nn64y
c
c		optimisation de nn64y. (voir plus haut)
c
	if(multy.gt.0.and.nnn64.gt.nr1) then
	n6y=(nnn64/nr1)*nr1
	mult1=nr1z1/n6y
	ireee2=nr1z1-mult1*n6y
	ii=0
	i2=0
	if(iresty.gt.0)ii=1
	if(ireee2.gt.0)i2=1
	if(mult1+i2.le.multy+ii) then
	multy=mult1
	nn64y=n6y
	nn65y=nn64y
	if((nn64y/8)*8.eq.nn64y)nn65y=nn64y+1
	iresty=ireee2
	endif
	endif
c
	if(multy.eq.0) then
	nn64y=nr1z1
	nn65y=nn64y
	if((nn64y/8)*8.eq.nn64y)nn65y=nn64y+1
	multy=1
	iresty=0
	endif
	lzz1=1
	lzz2=nn64y/nr1
	ldz=lzz2-lzz1
	if(iresty.gt.0) then
	n64y=iresty
	n65y=n64y
	if((n64y/8)*8.eq.n64y)n65y=n64y+1	
	endif
c
	n6565y=nn65y+nn65y
	n365y=n6565y+nn65y
	irr2y=nn64y-(nn64y/nr1)*nr1
c
c		preparation des elements necessares pour la transformation
c		du 3me indice du tableau
c
	nn64z=nnn64
	nn65z=nnn65
	nr1y1=nr1*ny1
	multz=(nr1y1)/nn64z
	irestz=nr1y1-multz*nn64z
c
c		optimisation de nn64z (voir plus haut)
c
	if(multz.gt.0.and.nnn64.gt.nr1) then
	n6z=(nnn64/nr1)*nr1
	mult1=nr1y1/n6z
	ireee2=nr1y1-mult1*n6z
	ii=0
	i2=0
	if(irestz.gt.0)ii=1
	if(ireee2.gt.0)i2=1
	if(mult1+i2.le.multz+ii) then
	multz=mult1
	nn64z=n6z
	nn65z=nn64z
	if((nn64z/8)*8.eq.nn64z)nn65z=nn64z+1
	irestz=ireee2
	endif
	endif
c
	if(multz.eq.0) then
	nn64z=nr1y1
	nn65z=nn64z
	if((nn64z/8)*8.eq.nn64z)nn65z=nn64z+1
	multz=1
	irestz=0
	endif
c
	if(irestz.gt.0) then
	n64z=irestz
	n65z=n64z
	if((n64z/8)*8.eq.n64z)n65z=n64z+1	
	endif
c
	lyy1=1
	lyy2=nn64z/nr1
	ldy=lyy2-lyy1
	irr2z=nn64z-(nn64z/nr1)*nr1
	n6565z=nn65z+nn65z
	n365z=n6565z+nn65z
c
 800	continue
c
c		transformation du premier indice
	lr1=0
	ly1=0
	lz1=0
c
c 1111111111111111111111111111111111111111111111111111111111111111111111111
c 1111111111111111111111111111111111111111111111111111111111111111111111111
c
		if(in1.eq.1) then
c
	if(ndeg(1).lt.4) return
c
	if(irr2r.eq.0)then
c
c		on effectue la transformation dans le cas nn64r multiple de ny1.
c	
c

	lr1=lrr1
	lr2=lrr2
c
	do 40 lmu=1,multr
	jr=0
	do 1 lr=1,nr1
	jz=jr	
	do 2 lz=lr1,lr2
	do 3 ly=1,ny1
	c64(ly+jz)=den(lr,ly,lz)
   3	continue
	jz=jz+ny1
   2	continue
	jr=jr+nn65r
   1	continue
c
	if(itch.eq.1) call tfmxs(nr1,nn64r,c64,cc)
	if(itch.eq.2) call chexms(nr,nn64r,c64,cc,cs)
c
c		on stocke les coefficients de la transformation dans 
c		den(lr,ly,lz).
c
	if(itch.eq.1) then
	jr=0
	do 4 lr=1,2
	jz=jr	
	do 5 lz=lr1,lr2
	do 6 ly=1,ny1
  	den(lr,ly,lz)=cc(ly+jz)
   6	continue
	jz=jz+ny1
   5	continue
	jr=jr+n6565r
   4	continue
c
	do 7 lr=4,nr1,2
	jz=jr	
	do 8 lz=lr1,lr2
	do 9 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
   9	continue
	jz=jz+ny1
   8	continue
	jr=jr+n6565r
   7	continue
c
	jr=n365r
	do 10 lr=3,nr1,2
	jz=jr	
	do 11 lz=lr1,lr2
	do 12 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  12	continue
	jz=jz+ny1
  11	continue
	jr=jr+n6565r
  10	continue
		endif
c
		if(itch.eq.2) then
c
	jr=0
	do 13 lr=1,nr1
	jz=jr	
	do 14 lz=lr1,lr2
	do 15 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  15	continue
	jz=jz+ny1
  14	continue
	jr=jr+nn65r
  13	continue
		endif
c
	lr1=lr2+1
	lr2=lr1+ldr
  40	continue
c
c		le calcul est continue' si ny1*nr1 n'est pas un multiple
c		de n64y.
c
	if(irestr.gt.0) then
c
	jr=0
	do 16 lr=1,nr1
	jz=jr
	do 17 lz=lr1,nz1
	do 18 ly=1,ny1
	c64(ly+jz)=den(lr,ly,lz)
  18	continue
	jz=jz+ny1
  17	continue
	jr=jr+n65r
  16	continue
c
	if(itch.eq.1) call tfmxs(nr1,n64r,c64,cc)
	if(itch.eq.2) call chexms(nr,n64r,c64,cc,cs)
c
c		on reintroduit les coeff. de fourier dans den(lr,ly,lz)
c
	if(itch.eq.1)then
c
	n265=n65r+n65r
	jr=0
	do 19 lr=1,2
	jz=jr
	do 20 lz=lr1,nz1
	do 21 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  21	continue
  	jz=jz+ny1
  20	continue
	jr=jr+n265
  19	continue
c
	do 22 lr=4,nr1,2
	jz=jr
	do 23 lz=lr1,nz1
	do 24 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  24	continue
  	jz=jz+ny1
  23	continue
	jr=jr+n265
  22	continue
c
	jr=n265+n65r
	do 25 lr=3,nr1,2
	jz=jr
       	do 26 lz=lr1,nz1
	do 27 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  27	continue
	jz=jz+ny1
  26	continue
	jr=jr+n265
  25	continue
		endif
c
	if(itch.eq.2) then
	jr=0
	do 28 lr=1,nr1
	jz=jr
	do 29 lz=lr1,nz1
	do 30 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  30	continue
	jz=jz+ny1
  29	continue
	jr=jr+n65r
  28	continue
	endif
c
	endif
	return
	endif
c
c		calcul de la tf dans le cas ou n64y n'est pas un multiple de 
c		ny1.
c
	if(irr2r.gt.0) then
c		
	nr64r=nn64r
	n63y=ny1-nr64r+1
	nr65r=nn65r
	mu2=nr64r
	lz1=1
	ly1=1
	lrr=1
	do 99 lm=1,multr
	if(lrr.gt.ny1) lrr=1
	lrr1=0
	lz1=(mu2-nr64r)/ny1+1
	lz2=.99999+float(mu2)/ny1
	ly1=+mu2+n63y-lz1*ny1
	ly2=mu2-(lz2-1)*ny1
 	jr=0
c
	do 46 lr=1,nr1
	jjz=jr
	if(lz2.gt.lz1) then
	jz=jr-lrr+1
c
	do 41 ly=lrr,ny1
	c64(ly+jz)=den(lr,ly,lz1)
  41	continue
	jjz=ny1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 42 lz=lz1+1,lz2-1
	do 43 ly=1,ny1
	c64(ly+jz)=den(lr,ly,lz)
  43	continue
	jz=jz+ny1
  42	continue
	jjz=jz
	endif
c
	if(ly2.ge.ly1.or.lrr1.eq.1) then
c
	if(lrr1.eq.0) then
	jz=jjz+1-ly1
	do 44 ly=ly1,ly2
	c64(ly+jz)=den(lr,ly,lz2)
  44	continue
	endif
	if(lrr1.eq.1) then
	jz=jjz
  	do 45 ly=1,ly2
	c64(ly+jz)=den(lr,ly,lz2)
  45	continue
	endif
c
	lsy=ly2+1
	endif
	jr=jr+nr65r
  46	continue			
c
	if(itch.eq.1) call tfmxs(nr1,nr64r,c64,cc)
	if(itch.eq.2) call chexms(nr,nr64r,c64,cc,cs)
c
c		les valeures calclculees par tfmys sont memorisees dans 
c		den(lr,ly,lz).
c
	if(itch.eq.1) then
c
	n265=nr65r+nr65r
 	jr=0
	do 52 lr=1,2
	jjz=jr
	if(lz2.gt.lz1) then
c
	jz=jr-lrr+1
	do 47 ly=lrr,ny1
c
	den(lr,ly,lz1)=cc(ly+jz)
  47	continue
	jjz=ny1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 48 lz=lz1+1,lz2-1
	do 49 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  49	continue
	jz=jz+ny1
  48	continue
	jjz=jz
	endif
c
	if(ly2.ge.ly1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jz=jjz+1-ly1
	do 50 ly=ly1,ly2
	den(lr,ly,lz2)=cc(ly+jz)
  50	continue
	endif
c
	if(lrr1.eq.1) then
	jz=jjz
	do 51 ly=1,ly2
	den(lr,ly,lz2)=cc(ly+jz)
  51	continue
	endif
c
	lsy=ly2+1
	endif
	jr=jr+n265
  52	continue			
c
	do 58 lr=4,nr1,2
	jjz=jr
	if(lz2.gt.lz1) then
c
	jz=jr-lrr+1
	do 53 ly=lrr,ny1
	den(lr,ly,lz1)=cc(ly+jz)
  53	continue
	jjz=ny1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 54 lz=lz1+1,lz2-1
	do 55 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  55	continue
	jz=jz+ny1
  54	continue
	jjz=jz
	endif
c
	if(ly2.ge.ly1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jz=jjz+1-ly1
	do 56 ly=ly1,ly2
	den(lr,ly,lz2)=cc(ly+jz)
  56	continue
	endif
c
	if(lrr1.eq.1) then
	jz=jjz
	do 57 ly=1,ly2
	den(lr,ly,lz2)=cc(ly+jz)
  57	continue
	endif
c
	lsr=ly2+1
	endif
	jr=jr+n265
  58	continue			
c
	jr=n265+nr65r
	do 65 lr=3,nr1,2
	jjz=jr
	if(lz2.gt.lz1) then
c
	jz=jr-lrr+1
	do 60 ly=lrr,ny1
	den(lr,ly,lz1)=cc(ly+jz)
  60	continue
	jjz=ny1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 61 lz=lz1+1,lz2-1
	do 62 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  62	continue
	jz=jz+ny1
  61	continue
	jjz=jz
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jz=jjz+1-ly1
	do 63 ly=ly1,ly2
	den(lr,ly,lz2)=cc(ly+jz)
  63	continue
	endif
c
	if(lrr1.eq.1) then
	jz=jjz
	do 64 ly=1,ly2
	den(lr,ly,lz2)=cc(ly+jz)
  64	continue
	endif
c
	lsr=ly2+1
	endif
	jr=jr+n265
  65	continue			
c
		endif
c
		if(itch.eq.2) then
 	jr=0
	do 71 lr=1,nr1
	jjz=jr
	if(lz2.gt.lz1) then
c
	jz=jr-lrr+1
	do 66 ly=lrr,ny1
	den(lr,ly,lz1)=cc(ly+jz)
  66	continue
	jjz=ny1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 67 lz=lz1+1,lz2-1
	do 68 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  68	continue
	jz=jz+ny1
  67	continue
	jjz=jz
	endif
c
	if(ly2.ge.ly1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jz=jjz+1-ly1
	do 69 ly=ly1,ly2
	den(lr,ly,lz2)=cc(ly+jz)
  69	continue
	endif
c
	if(lrr1.eq.1) then
	jz=jjz
	do 70 ly=1,ly2
	den(lr,ly,lz2)=cc(ly+jz)
  70	continue
	endif
c
	lsy=ly2+1
	endif
	jr=jr+nr65r
  71	continue			
c
		endif
c
	lrr=lsy	
	mu2=mu2+nr64r
  99	continue
	endif
c
	if(irestr.eq.0) return
c
	jr=0
	ly1=lrr
	lz1=lz2
	lz3=lz2
	if(lrr.gt.ny1) then
	ly1=1
	lz3=lz2+1
	endif
	nr64r=irestr
	nr65r=nr64r
	if((nr64r/8)*8.eq.nr64r)nr65r=nr64r+1
	n265=nr65r+nr65r
	jr=-ly1+1
c
	do 76 lr=1,nr1
c
	jz=jr
	if(lrr.le.ny1) then
	do 73 ly=ly1,ny1
	c64(ly+jz)=den(lr,ly,lz2)
  73	continue
	lz3=1+lz2
	jz=jz+ny1
	endif
c
	if(lz3.le.nz1)then
	do 74 lz=lz3,nz1
	do 75 ly=1,ny1
	c64(ly+jz)=den(lr,ly,lz)
  75	continue
	jz=jz+ny1
  74	continue
	endif  
	jr=jr+nr65r
  76	continue
c	
	if(itch.eq.1) call tfmxs(nr1,nr64r,c64,cc)
	if(itch.eq.2) call chexms(nr,nr64r,c64,cc,cs)
c		
c		on stockes les coefficients de fourier dans den(lr,ly,lz).
c
		if(itch.eq.1) then
c
	jr=-ly1+1
c
	do 80 lr=1,2
	jz=jr
	if(lrr.le.ny1) then
	do 77 ly=ly1,ny1
	den(lr,ly,lz2)=cc(ly+jz)
  77	continue
	lz3=1+lz2
	jz=jz+ny1
	endif
c
	if(lz3.le.nz1)then
	do 78 lz=lz3,nz1
	do 79 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  79	continue
	jz=jz+ny1
  78	continue
	endif  
	jr=jr+n265
  80	continue
c
	do 84 lr=4,nr1,2
	jz=jr
	if(lrr.le.ny1) then
	do 81 ly=ly1,ny1
	den(lr,ly,lz2)=cc(ly+jz)
  81	continue
	lz3=1+lz2
	jz=jz+ny1
	endif
c
	if(lz3.le.nz1)then
	do 82 lz=lz3,nz1
	do 83 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  83	continue
	jz=jz+ny1
  82	continue
	endif  
	jr=jr+n265
  84	continue
c
	jr=-ly1+1+nr65r+n265
c
	do 88 lr=3,nr1,2
	jz=jr
	if(lrr.le.ny1) then
	do 85 ly=ly1,ny1
	den(lr,ly,lz2)=cc(ly+jz)
  85	continue
	lz3=1+lz2
	jz=jz+ny1
	endif
c
	if(lz3.le.nz1)then
	do 86 lz=lz3,nz1
	do 87 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  87	continue
	jz=jz+ny1
  86	continue
	endif  
	jr=jr+n265
  88	continue
		endif
c
		if(itch.eq.2) then
c
	jr=-ly1+1
c
	do 93 lr=1,nr1
c
	jz=jr
	if(lrr.le.ny1) then
	do 90 ly=ly1,ny1
	den(lr,ly,lz2)=cc(ly+jz)
  90	continue
	lz3=1+lz2
	jz=jz+ny1
	endif
c
	if(lz3.le.nz1)then
	do 91 lz=lz3,nz1
	do 92 ly=1,ny1
	den(lr,ly,lz)=cc(ly+jz)
  92	continue
	jz=jz+ny1
  91	continue
	endif  
	jr=jr+nr65r
  93	continue
c
		endif
c

	return
	endif
c

c		tranformee du 2me indice.
c
c222222222222222222222222222222222222222222222222222222222222222222222222.
c222222222222222222222222222222222222222222222222222222222222222222222222.
c
		if(in1.eq.2) then
c
	if(ndeg(2).lt.4) return
	if(ny.lt.2) return
c
	if(irr2y.eq.0)then
c
c		on effectue la transformation dans le cas nn64y multiple de nr1.
c	
	lz1=lzz1
	lz2=lzz2
c
	do 115 lmu=1,multy
	jy=0
	do  102 ly=1,ny1
	jz=jy	
	do 100 lz=lz1,lz2
	do 101 lr=1,nr1
	c64(lr+jz)=den(lr,ly,lz)
  101	continue
	jz=jz+nr1
  100	continue
	jy=jy+nn65y
  102	continue
c
	if(itch.eq.1) call tfmys(ny1,nn64y,c64,cc)
	if(itch.eq.2) call cheyms(ny,nn64y,c64,cc,cs)
c
c		on stocke les coefficients de la transformation dans 
c		den(lr,ly,lz).
c
c
	if(itch.eq.1) then
	jy=0
	do 105 ly=1,2
	jz=jy	
	do 103 lz=lz1,lz2
	do 104 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 104	continue
	jz=jz+nr1
 103	continue
	jy=jy+n6565y
 105	continue
c
	do 108 ly=4,ny1,2
	jz=jy	
	do 106 lz=lz1,lz2
	do 107 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 107	continue
	jz=jz+nr1
 106	continue
	jy=jy+n6565y
 108	continue
c
	jy=n365y
	do 111 ly=3,ny1,2
	jz=jy	
	do 109 lz=lz1,lz2
	do 110 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 110	continue
	jz=jz+nr1
 109	continue
	jy=jy+n6565y
 111	continue
		endif
c
		if(itch.eq.2) then
c
	jy=0
	do 114 ly=1,ny1
	jz=jy	
	do 112 lz=lz1,lz2
	do 113 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 113	continue
	jz=jz+nr1
 112	continue
	jy=jy+nn65y
 114	continue
		endif
c
	lz1=lz2+1
	lz2=lz1+ldz
  115	continue
c
c		le calcul est continue' si nr1*ny1 n'est pas un multiple
c		de n64y.
c
	if(iresty.gt.0) then
c
	jy=0
	do 118 ly=1,ny1
	jz=jy
	do 116 lz=lz1,nz1
	do 117 lr=1,nr1
	c64(lr+jz)=den(lr,ly,lz)
 117	continue
	jz=jz+nr1
 116	continue
	jy=jy+n65y
 118	continue
c
	if(itch.eq.1) call tfmys(ny1,n64y,c64,cc)
	if(itch.eq.2) call cheyms(ny,n64y,c64,cc,cs)
c
c		on reintroduit les coeff. de fourier dans den(lr,ly,lz)
c
	if(itch.eq.1)then
c
	n265=n65y+n65y
	jy=0
	do 121 ly=1,2
	jz=jy
	do 119 lz=lz1,nz1
	do 120 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 120	continue
	jz=jz+nr1
 119	continue
	jy=jy+n265
 121	continue
c
	do 124 ly=4,ny1,2
	jz=jy
	do 122 lz=lz1,nz1
	do 123 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 123	continue
	jz=jz+nr1
 122	continue
	jy=jy+n265
 124	continue
c
	jy=n265+n65y
	do 127 ly=3,ny1,2
	jz=jy
	do 125 lz=lz1,nz1
	do 126 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 126	continue
	jz=jz+nr1
 125	continue
	jy=jy+n265
 127	continue
		endif
c
	if(itch.eq.2) then
	jy=0
	do 130  ly=1,ny1
	jz=jy
	do 128 lz=lz1,nz1
	do 129 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 129	continue
	jz=jz+nr1
 128	continue
	jy=jy+n65y
 130	continue
	endif
c
	endif
	return
	endif
c
c		calcul de la tf dans le cas ou n64y n'est pas un multiple de 
c		nr1.
c
	if(irr2y.gt.0) then
c		
	nr64y=nn64y
	n63r=nr1-nr64y+1
	nr65y=nn65y
	mu2=nr64y
	lz1=1
	lr1=1
	lrr=1
	do 189 lm=1,multy
	if(lrr.gt.nr1) lrr=1
	lrr1=0
	lz1=(mu2-nr64y)/nr1+1
	lz2=.99999+float(mu2)/nr1
	lr1=+mu2+n63r-lz1*nr1
	lr2=mu2-(lz2-1)*nr1
 	jy=0
c
	do 136 ly=1,ny1
	jjz=jy
	if(lz2.gt.lz1) then
	jz=jy-lrr+1
c
	do 131 lr=lrr,nr1
	c64(lr+jz)=den(lr,ly,lz1)
 131	continue                               
	jjz=nr1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 132 lz=lz1+1,lz2-1
	do 133 lr=1,nr1
	c64(lr+jz)=den(lr,ly,lz)
 133	continue                               
	jz=jz+nr1
 132	continue                               
	jjz=jz
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
c
	if(lrr1.eq.0) then
	jz=jjz+1-lr1
	do 134 lr=lr1,lr2
	c64(lr+jz)=den(lr,ly,lz2)
 134	continue                               
	endif
	if(lrr1.eq.1) then
	jz=jjz
	do 135 lr=1,lr2
	c64(lr+jz)=den(lr,ly,lz2)
 135	continue                               
	endif
c
	lsr=lr2+1
	endif
	jy=jy+nr65y
 136	continue                               
c
	if(itch.eq.1) call tfmys(ny1,nr64y,c64,cc)
 	if(itch.eq.2) call cheyms(ny,nr64y,c64,cc,cs)
c
c		les valeures calclculees par tfmys sont memorisees dans 
c		den(lr,ly,lz).
c
	if(itch.eq.1) then
c
	n265=nr65y+nr65y
 	jy=0
	do 142 ly=1,2
	jjz=jy
	if(lz2.gt.lz1) then
c
	jz=jy-lrr+1
	do 137 lr=lrr,nr1
	den(lr,ly,lz1)=cc(lr+jz)
 137	continue                               
	jjz=nr1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 138 lz=lz1+1,lz2-1
	do 139 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 139	continue                               
	jz=jz+nr1
 138	continue                               
	jjz=jz
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jz=jjz+1-lr1
	do 140 lr=lr1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
 140	continue                               
	endif
c
	if(lrr1.eq.1) then
	jz=jjz
	do 141 lr=1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
 141	continue                               
	endif
c
	lsr=lr2+1
	endif
	jy=jy+n265
 142	continue                               
c
	do 149 ly=4,ny1,2 
	jjz=jy
	if(lz2.gt.lz1) then
c
	jz=jy-lrr+1
	do 143 lr=lrr,nr1
	den(lr,ly,lz1)=cc(lr+jz)
 143	continue                               
	jjz=nr1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 144 lz=lz1+1,lz2-1
	do 145 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 145	continue                               
	jz=jz+nr1
 144	continue                               
	jjz=jz
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jz=jjz+1-lr1
	do 147 lr=lr1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
 147	continue                               
	endif
c
	if(lrr1.eq.1) then
	jz=jjz
	do 148 lr=1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
 148	continue                               
	endif
c
	lsr=lr2+1
	endif
	jy=jy+n265
 149	continue                               
c
	jy=n265+nr65y
	do 155 ly=3,ny1,2
	jjz=jy
	if(lz2.gt.lz1) then
c
	jz=jy-lrr+1
	do 150 lr=lrr,nr1
	den(lr,ly,lz1)=cc(lr+jz)
 150	continue            
	jjz=nr1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 151 lz=lz1+1,lz2-1
	do 152 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 152	continue            
	jz=jz+nr1
 151	continue            
	jjz=jz
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jz=jjz+1-lr1
	do 153 lr=lr1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
 153	continue            
	endif
c
	if(lrr1.eq.1) then
	jz=jjz
	do 154 lr=1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
 154	continue            
	endif
c
	lsr=lr2+1
	endif
	jy=jy+n265
 155	continue            
c
		endif
c
		if(itch.eq.2) then
 	jy=0
	do 162 ly=1,ny1
	jjz=jy
	if(lz2.gt.lz1) then
c
	jz=jy-lrr+1
	do 157 lr=lrr,nr1
	den(lr,ly,lz1)=cc(lr+jz)
 157	continue          
	jjz=nr1+jz
	lrr1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 158 lz=lz1+1,lz2-1
	do 159 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 159	continue          
	jz=jz+nr1
 158	continue          
	jjz=jz
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jz=jjz+1-lr1
	do 160 lr=lr1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
 160	continue          
	endif
c
	if(lrr1.eq.1) then
	jz=jjz
	do 161 lr=1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
 161	continue          
	endif
c
	lsr=lr2+1
	endif
	jy=jy+nr65y
 162	continue          
c
		endif
c
	lrr=lsr	
	mu2=mu2+nr64y
 189	continue                               
	endif
c
	if(iresty.eq.0) return
c
	jy=0
	lr1=lrr
	lz1=lz2
	lz3=lz2
	if(lrr.gt.nr1) then
	lr1=1
	lz3=lz2+1
	endif
	nr64y=iresty
	nr65y=nr64y
	if((nr64y/8)*8.eq.nr64y)nr65y=nr64y+1
	n265=nr65y+nr65y
	jy=-lr1+1
c
	do 166 ly=1,ny1
c
	jz=jy
	if(lrr.le.nr1) then
	do 163 lr=lr1,nr1
	c64(lr+jz)=den(lr,ly,lz2)
 163	continue          
	lz3=1+lz2
	jz=jz+nr1
	endif
c
	if(lz3.le.nz1)then
	do 164 lz=lz3,nz1
	do 165 lr=1,nr1
	c64(lr+jz)=den(lr,ly,lz)
 165	continue          
	jz=jz+nr1
 164	continue          
	endif  
	jy=jy+nr65y
 166	continue          
c	
	if(itch.eq.1) call tfmys(ny1,nr64y,c64,cc)
	if(itch.eq.2) call cheyms(ny,nr64y,c64,cc,cs)
c		
c		on stockes les coefficients de fourier dans den(lr,ly,lz).
c
		if(itch.eq.1) then
c
	jy=-lr1+1
c
	do 170 ly=1,2
	jz=jy
	if(lrr.le.nr1) then
	do 167 lr=lr1,nr1
	den(lr,ly,lz2)=cc(lr+jz)
 167	continue          
	lz3=1+lz2
	jz=jz+nr1
	endif
c
	if(lz3.le.nz1)then
	do 168 lz=lz3,nz1
	do 169 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 169	continue          
	jz=jz+nr1
 168	continue          
	endif  
	jy=jy+n265
 170	continue          
c
	do 174 ly=4,ny1,2
	jz=jy
	if(lrr.le.nr1) then
	do 171 lr=lr1,nr1
	den(lr,ly,lz2)=cc(lr+jz)
 171	continue          
	lz3=1+lz2
	jz=jz+nr1
	endif
c
	if(lz3.le.nz1)then
	do 172 lz=lz3,nz1
	do 173 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 173	continue          
	jz=jz+nr1
 172	continue          
	endif  
	jy=jy+n265
 174	continue          
c
	jy=-lr1+1+nr65y+n265
c
	do 179 ly=3,ny1,2
	jz=jy
	if(lrr.le.nr1) then
	do 175 lr=lr1,nr1
	den(lr,ly,lz2)=cc(lr+jz)
 175	continue          
	lz3=1+lz2
	jz=jz+nr1
	endif
c
	if(lz3.le.nz1)then
	do 176 lz=lz3,nz1
	do 177 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 177	continue          
	jz=jz+nr1
 176	continue          
 	endif  
	jy=jy+n265
 179	continue          
		endif
c
		if(itch.eq.2) then
c
	jy=-lr1+1
	do 183 ly=1,ny1
c
	jz=jy
	if(lrr.le.nr1) then
	do 180 lr=lr1,nr1
	den(lr,ly,lz2)=cc(lr+jz)
 180	continue          
	lz3=1+lz2
	jz=jz+nr1
	endif
c
	if(lz3.le.nz1)then
	do 181 lz=lz3,nz1
	do 182 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
 182	continue          
	jz=jz+nr1
 181	continue          
	endif  
	jy=jy+nr65y
 183	continue          
c
		endif
c

	return
	endif
c
c
c		transformation du 3eme index.
c	
	if(in1.eq.3) then
c
c
c
c		tranformee du 3me indice.
c
c 333333333333333333333333333333333333333333333333333333333333333333333333***
c 333333333333333333333333333333333333333333333333333333333333333333333333***
c 
		if(in1.eq.3) then
	if(ndeg(3).lt.4) return

c
	if(irr2z.eq.0)then
c
c		on effectue la transformation dans le cas n64z multiple de nr1.
c	
	ly1=lyy1
	ly2=lyy2
c
	do 217 lmu=1,multz
	jz=0
	do  202 lz=1,nz1
	jy=jz	
	do 200 ly=ly1,ly2
	do 201 lr=1,nr1
	c64(lr+jy)=den(lr,ly,lz)
 201	continue                        
	jy=jy+nr1
 200	continue                        
	jz=jz+nn65z
 202	continue                        
c
	if(itch.eq.1) call tfmzs(nz1,nn64z,c64,cc)
	if(itch.eq.2) call chezms(nz,nn64z,c64,cc,cs)
c
c		on stocke les coefficients de la transformation dans 
c		den(lr,ly,lz).
c
		if(itch.eq.1) then
c
	jz=0
	do 205 lz=1,2
	jy=jz	
	do 203 ly=ly1,ly2
	do 204 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 204	continue                        
	jy=jy+nr1
 203	continue                        
	jz=jz+n6565z
 205	continue                        
c
	do 209 lz=4,nz1,2
	jy=jz	
	do 206 ly=ly1,ly2
	do 207 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 207	continue                        
	jy=jy+nr1
 206	continue                        
	jz=jz+n6565z
 209	continue                        
c	
	jz=n365z
	do 212 lz=3,nz1,2
	jy=jz
	do 210 ly=ly1,ly2
	do 211 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 211	continue                        
	jy=jy+nr1
 210	continue                        
	jz=jz+n6565z
 212	continue                        
c
		endif
c
	if(itch.eq.2) then
	jz=0
	do 216 lz=1,nz1
	jy=jz	
	do 213 ly=ly1,ly2
	do 214 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 214	continue                        
	jy=jy+nr1
 213	continue                        
 	jz=jz+nn65z
 216	continue                        
		endif
c
	ly1=ly2+1
	ly2=ly1+ldy
 217	continue                        
c
c		le calcul est continue' si nr1*ny1 n'est pas un multiple
c		de n64z.
c
	if(irestz.gt.0) then
	n265=n65z+n65z
	jz=0
	do 238 lz=1,nz1
	jy=jz
	do 218 ly=ly1,ny1
	do 219 lr=1,nr1
	c64(lr+jy)=den(lr,ly,lz)
 219	continue                        
	jy=jy+nr1
 218	continue                        
	jz=jz+n65z
 238	continue                        
c
	if(itch.eq.1) call tfmzs(nz1,n64z,c64,cc)
	if(itch.eq.2) call chezms(nz,n64z,c64,cc,cs)
c
c		on reintroduit les coeff. de fourier dans den(lr,ly,lz)
c
		if(itch.eq.1) then
	jz=0
	do 222 lz=1,2
	jy=jz	
	do 220 ly=ly1,ny1
	do 221 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 221	continue                        
	jy=jy+nr1
 220	continue                        
	jz=jz+n265
 222	continue                        
c
	do 225 lz=4,nz1,2
	jy=jz	
	do 223 ly=ly1,ny1
	do 224 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 224	continue                        
	jy=jy+nr1
 223	continue                        
	jz=jz+n265
 225	continue                        
c	
	jz=n265+n65z
	do 228 lz=3,nz1,2
	jy=jz
	do 226 ly=ly1,ny1
	do 227 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 227	continue                        
	jy=jy+nr1
 226	continue                        
	jz=jz+n265
 228	continue                        
c
		endif
c
		if(itch.eq.2) then
c
	jz=0
	do 231 lz=1,nz1
	jy=jz
	do 229 ly=ly1,ny1
	do 230 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 230	continue                        
	jy=jy+nr1
 229	continue                        
	jz=jz+n65z
 231	continue                        
	endif
c
	endif
	return
	endif
c
c		calcul de la tf dans le cas ou n64z n'est pas un multiple de 
c		nr1.
c
	if(irr2z.gt.0) then
c		
	nr64z=nn64z
	n63r=nr1-nr64z+1
	nr65z=nn65z
	mu2=nr64z
	ly1=1
	lr1=1
	lrr=1
c
	do 289 lm=1,multz
	if(lrr.gt.nr1) lrr=1
	lrr1=0
	ly1=(mu2-nr64z)/nr1+1
	ly2=.99999+float(mu2)/nr1
	lr1=+mu2+n63r-ly1*nr1
	lr2=mu2-(ly2-1)*nr1
 	jz=0
c
	do 237 lz=1,nz1
	jjy=jz
	if(ly2.gt.ly1) then
	jy=jz-lrr+1
c
	do 232 lr=lrr,nr1
	c64(lr+jy)=den(lr,ly1,lz)
 232	continue         
	jjy=nr1+jy
	lrr1=1
	endif
c
	if(ly2.gt.ly1+1) then
	jy=jjy
	do 233 ly=ly1+1,ly2-1
	do 234 lr=1,nr1
	c64(lr+jy)=den(lr,ly,lz)
 234	continue         
	jy=jy+nr1
 233	continue         
	jjy=jy
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
c
	if(lrr1.eq.0) then
	jy=jjy+1-lr1
	do 235 lr=lr1,lr2
	c64(lr+jy)=den(lr,ly2,lz)
 235	continue         
	endif
	if(lrr1.eq.1) then
	jy=jjy
	do 236 lr=1,lr2
	c64(lr+jy)=den(lr,ly2,lz)
 236	continue         
	endif
c
	lsr=lr2+1
	endif
	jz=jz+nr65z
 237	continue         
c
	if(itch.eq.1) call tfmzs(nz1,nr64z,c64,cc)
 	if(itch.eq.2) call chezms(nz,nr64z,c64,cc,cs)

c		les valeures calclculees par tfmzs sont memorisees dans 
c		den(lr,ly,lz).
c
		if(itch.eq.1) then
 	jz=0
	n265=nr65z+nr65z
c
	do 244 lz=1,2
	jjy=jz
	if(ly2.gt.ly1) then
c
	jy=jz-lrr+1
	do 239 lr=lrr,nr1
	den(lr,ly1,lz)=cc(lr+jy)
 239	continue
	jjy=nr1+jy
	lrr1=1
	endif
c
	if(ly2.gt.ly1+1) then
	jy=jjy
	do 240 ly=ly1+1,ly2-1
	do 241 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 241	continue                  
	jy=jy+nr1
 240	continue                  
	jjy=jy
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jy=jjy+1-lr1
	do 242 lr=lr1,lr2
	den(lr,ly2,lz)=cc(lr+jy)
 242	continue                  
	endif
c
	if(lrr1.eq.1) then
	jy=jjy
	do 243 lr=1,lr2
	den(lr,ly2,lz)=cc(lr+jy)
 243	continue                  
	endif
c
	lsr=lr2+1
	endif
 	jz=jz+n265
 244	continue                  
c
	do 250 lz=4,nz1,2
	jjy=jz
	if(ly2.gt.ly1) then
c
	jy=jz-lrr+1
	do 245 lr=lrr,nr1
	den(lr,ly1,lz)=cc(lr+jy)
 245	continue                  
	jjy=nr1+jy
	lrr1=1
	endif
c
	if(ly2.gt.ly1+1) then
	jy=jjy
	do 246 ly=ly1+1,ly2-1
	do 247 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 247	continue                  
	jy=jy+nr1
 246	continue
	jjy=jy
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jy=jjy+1-lr1
	do 248 lr=lr1,lr2
	den(lr,ly2,lz)=cc(lr+jy)
 248	continue                  
	endif
c
	if(lrr1.eq.1) then
	jy=jjy
	do 249 lr=1,lr2
	den(lr,ly2,lz)=cc(lr+jy)
 249	continue                  
	endif
c
	lsr=lr2+1
	endif
	jz=jz+n265
 250	continue                  
c
	jz=n265+nr65z
	do 256 lz=3,nz1,2
	jjy=jz
	if(ly2.gt.ly1) then
c
	jy=jz-lrr+1
	do 251 lr=lrr,nr1
	den(lr,ly1,lz)=cc(lr+jy)
 251	continue                     
	jjy=nr1+jy
	lrr1=1
	endif
c
	if(ly2.gt.ly1+1) then
	jy=jjy
	do 252 ly=ly1+1,ly2-1
	do 253 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 253	continue                     
	jy=jy+nr1
 252	continue                     
	jjy=jy
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jy=jjy+1-lr1
	do 254 lr=lr1,lr2
 	den(lr,ly2,lz)=cc(lr+jy)
 254	continue                     
	endif
c
	if(lrr1.eq.1) then
	jy=jjy
	do 255 lr=1,lr2
	den(lr,ly2,lz)=cc(lr+jy)
 255	continue                     
	endif
c
	lsr=lr2+1
	endif
	jz=jz+n265
  256	continue                     
c
		endif
c
		if(itch.eq.2) then
 	jz=0
	do 262 lz=1,nz1
	jjy=jz
	if(ly2.gt.ly1) then
c
	jy=jz-lrr+1
	do 257 lr=lrr,nr1
	den(lr,ly1,lz)=cc(lr+jy)
 257	continue              
	jjy=nr1+jy
	lrr1=1
	endif
c
	if(ly2.gt.ly1+1) then
	jy=jjy
	do 258 ly=ly1+1,ly2-1
	do 259 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 259	continue              
	jy=jy+nr1
 258	continue              
	jjy=jy
	endif
c
	if(lr2.ge.lr1.or.lrr1.eq.1) then
	if(lrr1.eq.0) then
	jy=jjy+1-lr1
	do 260 lr=lr1,lr2
	den(lr,ly2,lz)=cc(lr+jy)
 260	continue              
	endif
c
	if(lrr1.eq.1) then
	jy=jjy
	do 261 lr=1,lr2
	den(lr,ly2,lz)=cc(lr+jy)
 261	continue              
	endif
c
	lsr=lr2+1
	endif
	jz=jz+nr65z
 262	continue              
c
		endif
c
	lrr=lsr	
c
	mu2=mu2+nr64z
 289	continue
	endif
c
	if(irestz.eq.0) return
	jz=0
	lr1=lrr
	ly1=ly2
	ly3=ly2
	if(lrr.gt.nr1) then
	lr1=1
	ly3=ly2+1
	endif
	nr64z=irestz
	nr65z=nr64z
	if((nr64z/8)*8.eq.nr64z)nr65z=nr64z+1
	n265=nr65z+nr65z
	jz=-lr1+1
c
	do 293 lz=1,nz1
c
	jy=jz
	if(lrr.le.nr1) then
	do 290 lr=lr1,nr1
	c64(lr+jy)=den(lr,ly2,lz)
 290	continue          
	ly3=1+ly2
	jy=jy+nr1
	endif
c
	if(ly3.le.ny1)then
	do 291 ly=ly3,ny1
	do 292 lr=1,nr1
	c64(lr+jy)=den(lr,ly,lz)
 292	continue          
	jy=jy+nr1
 291	continue          
	endif  
	jz=jz+nr65z
 293	continue          
c	
	if(itch.eq.1) call tfmzs(nz1,nr64z,c64,cc)
	if(itch.eq.2) call chezms(nz,nr64z,c64,cc,cs)
c		
c		on stockes les coefficients de fourier dans den(lr,ly,lz).
c
		if(itch.eq.1) then
c
	jz=-lr1+1
	do 297 lz=1,2
	jy=jz
	if(lrr.le.nr1) then
	do 294 lr=lr1,nr1
	den(lr,ly2,lz)=cc(lr+jy)
 294	continue          
	ly3=1+ly2
	jy=jy+nr1
	endif
c
	if(ly3.le.ny1)then
	do 295 ly=ly3,ny1
	do 296 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 296	continue          
	jy=jy+nr1
 295	continue          
	endif  
 	jz=jz+n265
 297	continue          
c
	do 301 lz=4,nz1,2
	jy=jz
	if(lrr.le.nr1) then
	do 298 lr=lr1,nr1
	den(lr,ly2,lz)=cc(lr+jy)
 298	continue                       
	ly3=1+ly2
	jy=jy+nr1
	endif
c
	if(ly3.le.ny1)then
	do 299 ly=ly3,ny1
	do 300 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 300	continue                       
	jy=jy+nr1
 299	continue                       
	endif  
	jz=jz+n265
 301	continue                       
c
	jz=1-lr1+nr65z+n265
c
	do 305 lz=3,nz1,2
	jy=jz
	if(lrr.le.nr1) then
	do 302 lr=lr1,nr1
	den(lr,ly2,lz)=cc(lr+jy)
  302	continue                       
	ly3=1+ly2
	jy=jy+nr1
	endif
c
	if(ly3.le.ny1)then
 	do 303 ly=ly3,ny1
	do 304 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 304	continue                       
	jy=jy+nr1
 303	continue                       
	endif  
	jz=jz+n265
 305	continue                       
		endif
c
		if(itch.eq.2) then
	jz=-lr1+1
c
	do 309 lz=1,nz1
c
	jy=jz
	if(lrr.le.nr1) then
	do 306 lr=lr1,nr1
	den(lr,ly2,lz)=cc(lr+jy)
 306	continue                       
	ly3=1+ly2
	jy=jy+nr1
	endif
c
	if(ly3.le.ny1)then
	do 307 ly=ly3,ny1
	do 308 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jy)
 308	continue                       
	jy=jy+nr1
 307	continue                       
	endif  
	jz=jz+nr65z
 309	continue                       
c
	endif
c
		endif
	endif
1000	format(1x,'f',10d12.4)
1010	format(1x,' ')
2000	format(1x,'fou',20i5)
	return
	end

