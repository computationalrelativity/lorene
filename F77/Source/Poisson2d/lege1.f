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



       subroutine lege1(n,ndim,m1,wp)
c
	implicit none
c
c           routine pour le calcul par rcourance ascendente
c           des fonctions associees de llegendre d'ordre m=m1-1
c           de dgree compris entre m et n dans l'interval
c           0.le.teta.le.pi.
c
C
C $Id$
C $Log$
C Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
C Initial revision
C
c Revision 1.2  1997/05/23  11:39:50  hyc
c *** empty log message ***
c
C Revision 1.1  1997/03/17 20:40:51  hyc
C Initial revision
C
C
C $Header$
C
C
	character*100 header
	data header/'$Header$'/

	REAL*8
     1	x1,x0,pig,pi,teta,cosen,prod,facto1,facto2,fac1,fac2
     1	,wp,x,se,sen

	integer ndim,n1,n2,n,n21,m,m1,m2,m22,jmax,l,isin,j,nj1
     1	,nm21,nmn,ncon

       dimension wp(ndim,*),sen(259),cosen(259)
c
	save ncon,n2,n21,n1,jmax,x0,pig,pi,cosen
c
	data ncon /0/
c     
	if(n.gt.258) then
	write(*,*)
     1	'dimensions insuffisantes dasns la routine lege1,n=',n
	call exit
	endif
c
	if(n.eq.ncon) go to 999
	ncon=n
	 n1=n+1
       n2=n/2
       n21=n2+1
       jmax=n1
       x0=0
       pig=acos(x0)
       pi=pig/n
       do 1 l=1,n1
       teta=(l-1)*pi+pig
       cosen(l)=cos(teta)
  1    continue
  999	continue
c
       m=m1-1
       m2=m+2
       m22=m+m
c
       prod=1
       facto1=prod
       facto2=prod*(m22+1)
       isin=(-1)**m
       fac1=facto1*isin
       fac2=facto2*isin
       do 3 l=1,n1
       wp(l,1)=fac1
       wp(l,2)=cosen(l)*fac2
  3    continue
c
       if(m.gt.0) then
       do 4 l=1,n1
       x=cosen(l)
       se=sqrt(1-x**2)**m
       wp(l,1)=wp(l,1)*se
       wp(l,2)=wp(l,2)*se
       sen(l)=se
  4    continue
             endif
c
c           calcul des polyn. par recourance ascendente.
c
       do 5 j=3,jmax
       nj1=m+j-2
       nm21=2*nj1+1
       nmn=nj1+m
       x1=1
       fac1=x1/(nj1+1-m)
  300  format(1x,3i5,10e11.3)
       do 6 l=1,n1
       wp(l,j)=(nm21*cosen(l)*wp(l,j-1)-nmn*wp(l,j-2))*fac1
  6    continue
  5    continue
  100  format(1x,10d12.4)
  101  format(1x,' ')
       return
       end
