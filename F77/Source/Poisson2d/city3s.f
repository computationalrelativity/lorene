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

	subroutine city3s(ndeg,ndimr,ndimy,nnn64,itch,idr,ind,c64,cc,cs,den)

c
c## version du 08.10.1993 : suppression des variables declarees mais non 
c			    utilisees
c
		implicit none
c
c		routine pour le calcul des transformees de fourier
c		et de tchebytchev d'un tableau a 3 dimensions.
c
c		routine completement craytinizee.
c
c		arguments de la routine:
c
c			ndeg	= tableau, ndeg(3) contenant les de-
c				  gres de liberte des transformees
c				  a effectuer, ndeg(1) concerne le pre-
c				  mier indice de la matrice, ndeg(2)
c				  le 2me indice, ndeg(3) le 3me indice de la
c				  matrice.					
c				  ndeg doit imperativement etre de la
c				  forme 2**m*3**p*5**q pour les trans-
c				  formees de fourier (m,p,q nombres 
c				  entiers)
c				  et 2**p*3**p*5*q+1 pour les transfor-
c				  mees de tchebytchev.
c
c		ndimy,ndimz  =dimesion du tableau yy(lr,ly,lz).
c		pour des raisons de craytinisation ndimy et ndimz ne doivent pas
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
c
c		idr	=prametre: si idr=1 la derivee premiere est effectuee
c			 si idr=2 la derivee 2me est executee.
c
c		in1	= parametre, si in1=1 la transformee est
c			 effectuee sur le premier indice, si in1=2 sur le
c			 deuxieme, si in1=3 sur le 3me,.
c			 par exemple itch=2, in1=1 signifie
c			 qu'on effectue la transformation de tcheby-
c			 tchev sur le premier indice du tableau yy.
c
c		ind	= parametre: en output on a les coefficients
c			  de fourier (ou de tchebytchev selon itch) des deri-
c			  vees (premieres ou deuxiemes selon idr) si ind=1.
c			  si ind=2 on a les derivees des fonctions.
c
c		c64,cc,cs= tableaux de travail: dimension minime=
c			   (nnn64+1)*((max(ndeg(1),ndeg(2))+3)
c
c		den	=tableau a 3 dimensions contenant la fonction
c			 a transformer en input, et la transformee en
c			 output. les coefficients de fourier (dans le cas
c			 d'une transformation de fourier) sont stockes
c			 dans la facon suivante: (exemple dans le cas
c			 de la transformation du 1er indice) dans den(1,l,m)
c			 il-y-a le coefficient de fourier cosinus de la 
c			 frequence zero, dans les coefficients pairs 
c			 (den(2,l,m),den(4,l,m)...den(2*n,l,m) les coefi-
c			 cients en cosinus, dans les coefficients impairs
c			 den(3,l,m),den(5,l,m),...den(2*n-1) les termes en 
c			 sinus. par consequent les termes den(2,l,m) et den(3,l,m)
c			 sont les termes du developpement de fourier correspon-
c			 dants a la meme frequence. en totale il-y-a 2*n
c			 degres de liberte.
c
c		routine ayant testee avec le protocol ordinaire le 10/1/1987
c		routine modifiee le 28/octobre 1994 - den(ndimr,ndimy,*)
c
C
C $Id$
C $Log$
C Revision 1.2  2005/11/02 08:59:33  m_bejger
C Commented out multiple initializations of data variables, line 129,130
C (issue known to produce internal compiler error with
C gcc-gfortran 4.0.0, 4.0.1)
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.2  1997/05/23  11:31:15  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:21:45  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	integer
     1	ndimr,ndimy,nr1,ny1,nz1,ndy,ndz,ndr,neq,nnn64,nr,ny,nz,nnn65
     1	,nn64y,nn65y,ireee2,multy,ndeg,nr1z1,iresty,n6y,mult1,ii
     1	,i2,ldy,lyy2,lyy1,n64y,n65y,n6565y,iyy2y,ly1,ly2,l
     1	,lmu,jy,ly,jz,lz,lr,n365y,itch,lyy,ny64y,lz3,initi
     1  ,ny65y,mu2,lz1,lr1,lm,lz2,lr2,n63r,jjz,lsr,n265,idr1
     1  ,n256,ide1,ide2,n25,idr2,ind,idr
c
	real*8 cc,c64,den,cs
c
c	data ndy,ndr,ndz/0,0,0/
c	data neq/0/
c
	dimension ide1(257),den(ndimr,ndimy,*),c64(*),cc(*),cs(*),ndeg(3)
	dimension ide2(257)
	data initi,neq,ndy,ndr,ndz/0,0,0,0,0/

	save	n256,n25,ide1,ide2,initi,nr1,ny1,nz1,ndr,ndy,ndz,neq
	save	nr,nz,ny,nnn65,nn64y,nn65y,nr1z1,multy,iresty,n6y
	save	mult1,ireee2,ii,i2,lyy1,lyy2,ldy,n6565y,n365y,iyy2y
	save	n64y,n65y

c
c		initialisation.
c
	if(initi.ne.314) then
c
	n256=256
	n25=n256/2
	do 700 l=1,n25
	ide1(l+l)=-l
	ide1(l+l+1)=l
	ide2(l+l)=-l*l
	ide2(l+l+1)=-l*l
  700	continue
	initi=314
	endif
c
	nr1=ndeg(1)
	ny1=ndeg(2)
	nz1=ndeg(3)
c
	if(ndy.eq.ny1.and.ndz.eq.nz1.and.ndr.eq.nr1.and.neq.eq.nnn64) go to 800
c
	if(nr1.gt.n256.or.ny1.gt.n256.or.nz1.gt.n256) then
	print 2000,nr1,ny1,nz1
2000	format(1x,'DIMENSIONS INSUFF. DANS LA SUB. CITY3S: NR1,NY1,NZ1='
     1	,3i3)
	call exit
	endif
c
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
c			preparation des quantites necessaires pour la transfor-
c		mation du 1er index du tableau.
c
	nn64y=nnn64
	nn65y=nnn65
	nr1z1=nr1*nz1
	multy=(nr1z1)/nn64y
	iresty=nr1z1-multy*nn64y
c
c		optimisation de nn64y: on cherche une valeur de nn64y qui
c		soit < nnn64 mais un multiple de ny1. cela en vue de reduire
c		le nombre d'operations dans le transfert des valeurs de den
c		dans c64.
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

	if(multy.eq.0) then
	nn64y=nr1z1
	nn65y=nn64y
	if((nn64y/8)*8.eq.nn64y) nn65y=nn64y+1
	multy=1
	iresty=0
	endif
	lyy1=1
	lyy2=nn64y/nr1
	ldy=lyy2-lyy1
c
	if(iresty.gt.0) then
	n64y=iresty
	n65y=n64y
	if((n64y/8)*8.eq.n64y)n65y=n64y+1	
	endif
c
	n6565y=nn65y+nn65y
	n365y=n6565y+nn65y
	iyy2y=nn64y-(nn64y/nr1)*nr1
c
 800	continue
c
c22222222222222222222222222222222222222222222222222222222222222222222222222222222
c222222222222222222222222222222222222222222222222222222222222222222222222222
c
c		transformation du 2eme indice
c
c	calcul des coefficients de la derivee premiere
c
	idr1=1
	idr2=1
	if(idr.eq.6) idr1=-1
	if(idr.eq.4) idr2=-1
c
	if(iyy2y.eq.0)then
c
c		on effectue la transformation dans le cas nn64y multiple de nr1.
c
	ly1=lyy1
	ly2=lyy2
c
	do 40 lmu=1,multy
c
	jy=0
	do 13 ly=1,ny1
	jz=jy	
	do 14 lz=ly1,ly2
	do 15 lr=1,nr1
  	cc(lr+jz)=den(lr,ly,lz)
  15	continue
	jz=jz+nr1
  14	continue
	jy=jy+nn65y
  13	continue
c
	if(idr.eq.0) call citxms(ny,nn64y,itch,cc,cs,c64)
c
	if(idr.eq.1) then
	if(itch.eq.0) call cd1mrs(ny,nn64y,2,cc,c64,1,cs)
	if(itch.eq.1) call derfms(ny,nn64y,cc,cs,c64)
	endif
c
	if(idr.eq.2) call cd2mrs(ny,nn64y,itch,cc,cs,c64)
c
	if(idr.eq.3.or.idr.eq.4) call dir2ms(ny,itch,idr2,nn64y,cc,cs,c64)
	if(idr.eq.5.or.idr.eq.6) call dimras(ny,itch,idr1,nn64y,cc,c64)
	if(idr.eq.7) call cd1mrs(nr,nn64y,itch,cc,c64,1,cs)
c
c		les valeurs de la fonction sont stockees dans den
c
	if(idr.ne.0.and.ind.eq.2) then
c
	call citxms(ny,nn64y,itch,c64,cs,cc)
	jy=0
	do 1 ly=1,ny1
	jz=jy	
	do 2 lz=ly1,ly2
	do 3 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
   3	continue
	jz=jz+nr1
   2	continue
	jy=jy+nn65y
   1	continue
		go to 801
c
	endif
c
	jy=0
	do 413 ly=1,ny1
	jz=jy	
	do 415 lz=ly1,ly2
	do 414 lr=1,nr1
	den(lr,ly,lz)=c64(lr+jz)
  414	continue       
  	jz=jz+nr1
  415	continue       
  	jy=jy+nn65y
  413	continue       
  801	continue
c
	ly1=ly2+1
	ly2=ly1+ldy
  40	continue
c
	if(iresty.eq.0) return
c
c		le calcul est continue' si ny1*nr1 n'est pas un multiple
c		de n64y(cas ind=1).
c
	if(iresty.gt.0) then
c
	jy=0
	do 28 ly=1,ny1
	jz=jy
	do 29 lz=ly1,nz1
	do 30 lr=1,nr1
	cc(lr+jz)=den(lr,ly,lz)
  30	continue
	jz=jz+nr1
  29	continue
	jy=jy+n65y
  28	continue
c
	if(idr.eq.0) call citxms(ny,n64y,itch,cc,cs,c64)
c
	if(idr.eq.1) then
	if(itch.eq.0) call cd1mrs(ny,n64y,2,cc,c64,1,cs)
	if(itch.eq.1) call derfms(ny,n64y,cc,cs,c64)
	endif
c
	if(idr.eq.2) call cd2mrs(ny,n64y,itch,cc,cs,c64)
c
	if(idr.eq.3.or.idr.eq.4) call dir2ms(ny,itch,idr2,n64y,cc,cs,c64)
	if(idr.eq.5.or.idr.eq.6) call dimras(ny,itch,idr1,n64y,cc,c64)
	if(idr.eq.7) call cd1mrs(ny,n64y,itch,cc,c64,1,cs)
c
	if(ind.eq.2.and.idr.ne.0) then
c
	call citxms(ny,n64y,itch,c64,cs,cc)
	jy=0
	do 416 ly=1,ny1
	jz=jy	
	do 417 lz=ly1,nz1
	do 418 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
  418	continue       
	jz=jz+nr1
  417	continue       
	jy=jy+n65y
  416	continue       
	return
	endif
c
	jy=0
	do 16 ly=1,ny1
	jz=jy
	do 17 lz=ly1,nz1
	do 18 lr=1,nr1
	den(lr,ly,lz)=c64(lr+jz)
  18	continue
	jz=jz+nr1
  17	continue
	jy=jy+n65y
  16	continue
c
	endif
	return
	endif
c
c		fin du calcul pour le cas nn64z  multiple de ny1 (cas ind=1)
c
c...........................................................................
c
c		calcul de la tf dans le cas ou n64y n'est pas un multiple de 
c		ny1.cas(ind=1).
c
	if(iyy2y.gt.0) then
c		
c		les coefficients de fourier sont stockes dans cc
c
	ny64y=nn64y
	n63r=nr1-ny64y+1
	ny65y=nn65y
	mu2=ny64y
	lz1=1
	lr1=1
	lyy=1
c
	do 99 lm=1,multy
	if(lyy.gt.nr1) lyy=1
	lyy1=0
	lz1=(mu2-ny64y)/nr1+1
	lz2=.99999+float(mu2)/nr1
	lr1=+mu2+n63r-lz1*nr1
	lr2=mu2-(lz2-1)*nr1
c
 	jy=0
	do 71 ly=1,ny1
	jjz=jy
	if(lz2.gt.lz1) then
c
	jz=jy-lyy+1
	do 66 lr=lyy,nr1
	cc(lr+jz)=den(lr,ly,lz1)
  66	continue
	jjz=nr1+jz
	lyy1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 67 lz=lz1+1,lz2-1
	do 68 lr=1,nr1
	cc(lr+jz)=den(lr,ly,lz)
  68	continue
	jz=jz+nr1
  67	continue
	jjz=jz
	endif
c
	if(lr2.ge.lr1.or.lyy1.eq.1) then
	if(lyy1.eq.0) then
	jz=jjz+1-lr1
	do 69 lr=lr1,lr2
	cc(lr+jz)=den(lr,ly,lz2)
  69	continue
	endif
c
	if(lyy1.eq.1) then
	jz=jjz
	do 70 lr=1,lr2
	cc(lr+jz)=den(lr,ly,lz2)
  70	continue
	endif
c
	lsr=lr2+1
	endif
	jy=jy+ny65y
  71	continue			
c
	if(idr.eq.0) call citxms(ny,ny64y,itch,cc,cs,c64)
c
	if(idr.eq.1) then
	if(itch.eq.0) call cd1mrs(ny,ny64y,2,cc,c64,1,cs)
	if(itch.eq.1) call derfms(ny,ny64y,cc,cs,c64)
	endif
c
	if(idr.eq.2) call cd2mrs(ny,ny64y,itch,cc,cs,c64)
c
	if(idr.eq.3.or.idr.eq.4) call dir2ms(ny,itch,idr2,ny64y,cc,cs,c64)
	if(idr.eq.5.or.idr.eq.6) call dimras(ny,itch,idr1,ny64y,cc,c64)
	if(idr.eq.7) call cd1mrs(ny,ny64y,itch,cc,c64,1,cs)
c
	if(ind.eq.2.and.idr.ne.0) then
c
	call citxms(ny,ny64y,itch,c64,cs,cc)
 	jy=0
	do 424 ly=1,ny1
	jjz=jy
	if(lz2.gt.lz1) then
	jz=jy-lyy+1
c
	do 419 lr=lyy,nr1
	den(lr,ly,lz1)=cc(lr+jz)
  419	continue    
	jjz=nr1+jz
	lyy1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 420 lz=lz1+1,lz2-1
	do 421 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
  421	continue    
	jz=jz+nr1
  420	continue    
	jjz=jz
	endif
c
	if(lr2.ge.lr1.or.lyy1.eq.1) then
c
	if(lyy1.eq.0) then
	jz=jjz+1-lr1
	do 422 lr=lr1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
  422	continue    
	endif
	if(lyy1.eq.1) then
	jz=jjz
  	do  423 lr=1,lr2
	den(lr,ly,lz2)=cc(lr+jz)
  423	continue    
	endif
c
	lsr=lr2+1
	endif
	jy=jy+ny65y
  424	continue    
		go to 802
	endif
 	jy=0
c
	do 46 ly=1,ny1
	jjz=jy
	if(lz2.gt.lz1) then
	jz=jy-lyy+1
c
	do 41 lr=lyy,nr1
	den(lr,ly,lz1)=c64(lr+jz)
  41	continue
	jjz=nr1+jz
	lyy1=1
	endif
c
	if(lz2.gt.lz1+1) then
	jz=jjz
	do 42 lz=lz1+1,lz2-1
	do 43 lr=1,nr1
	den(lr,ly,lz)=c64(lr+jz)
  43	continue
	jz=jz+nr1
  42	continue
	jjz=jz
	endif
c
	if(lr2.ge.lr1.or.lyy1.eq.1) then
c
	if(lyy1.eq.0) then
	jz=jjz+1-lr1
	do 44 lr=lr1,lr2
	den(lr,ly,lz2)=c64(lr+jz)
  44	continue
	endif
	if(lyy1.eq.1) then
	jz=jjz
  	do 45 lr=1,lr2
	den(lr,ly,lz2)=c64(lr+jz)
  45	continue
	endif
c
	lsr=lr2+1
	endif
	jy=jy+ny65y
  46	continue			
 802	continue
c
	lyy=lsr	
	mu2=mu2+ny64y
  99	continue
c
	if(iresty.eq.0) return
c
c		si n64r n'est pas un multiple exacte de nr1*nz1 le calcul
c		continue
c
	jy=0
	lr1=lyy
	lz1=lz2
	lz3=lz2
	if(lyy.gt.nr1) then
	lr1=1
	lz3=lz2+1
	endif
	ny64y=iresty
	ny65y=ny64y
	if((ny64y/8)*8.eq.ny64y)ny65y=ny64y+1
	n265=ny65y+ny65y
c
c		
c		les coefficients de fourier sont stockes dans cc
c
	jy=-lr1+1
c
	do 93 ly=1,ny1
c
	jz=jy
	if(lyy.le.nr1) then
	do 90 lr=lr1,nr1
	cc(lr+jz)=den(lr,ly,lz2)
  90	continue
	lz3=1+lz2
	jz=jz+nr1
	endif
c
	if(lz3.le.nz1)then
	do 91 lz=lz3,nz1
	do 92 lr=1,nr1
	cc(lr+jz)=den(lr,ly,lz)
  92	continue
	jz=jz+nr1
  91	continue
	endif  
	jy=jy+ny65y
  93	continue
c	
	if(idr.eq.0) call citxms(ny,ny64y,itch,cc,cs,c64)
c
	if(idr.eq.1) then
	if(itch.eq.0) call cd1mrs(ny,ny64y,2,cc,c64,1,cs)
	if(itch.eq.1) call derfms(ny,ny64y,cc,cs,c64)
	endif
c
	if(idr.eq.2) call cd2mrs(ny,ny64y,itch,cc,cs,c64)
c
	if(idr.eq.3.or.idr.eq.4) call dir2ms(ny,itch,idr2,ny64y,cc,cs,c64)
	if(idr.eq.5.or.idr.eq.6) call dimras(ny,itch,idr1,ny64y,cc,c64)
	if(idr.eq.7) call cd1mrs(ny,ny64y,itch,cc,c64,1,cs)
	if(ind.eq.2.and.idr.ne.0) then
c
	call citxms(ny,ny64y,itch,c64,cs,cc)
	jy=-lr1+1
	do 76 ly=1,ny1
	jz=jy
	if(lyy.le.nr1) then
	do 73 lr=lr1,nr1
	den(lr,ly,lz2)=cc(lr+jz)
  73	continue
	lz3=1+lz2
	jz=jz+nr1
	endif
c
	if(lz3.le.nz1)then
	do 74 lz=lz3,nz1
	do 75 lr=1,nr1
	den(lr,ly,lz)=cc(lr+jz)
  75	continue
	jz=jz+nr1
  74	continue
	endif  
	jy=jy+ny65y
  76	continue
	return
	endif
c		
	jy=-lr1+1
	do 428 ly=1,ny1
	jz=jy
	if(lyy.le.nr1) then
	do 425 lr=lr1,nr1
	den(lr,ly,lz2)=c64(lr+jz)
  425	continue
	lz3=1+lz2
	jz=jz+nr1
	endif
c
	if(lz3.le.nz1)then
	do 426 lz=lz3,nz1
	do 427 lr=1,nr1
	den(lr,ly,lz)=c64(lr+jz)
  427	continue
	jz=jz+nr1
  426	continue
	endif  
	jy=jy+ny65y
  428	continue
	endif
	return
 1000	format(1x,10e12.4)
 1010	format(1x,' ')
	end

