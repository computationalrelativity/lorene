/*
 *  Method of class Et_bin_bhns_extr to construct spherical harmonics
 *
 *    (see file et_bin_bhns_extr.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Joshua A. Faber
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

char et_bin_bhns_extr_ylm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/12/29 16:30:46  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>
#include <math.h>

// Lorene headers
#include "map.h"
#include "tenseur.h"
#include "et_bin_bhns_extr.h"

void Et_bin_bhns_extr::get_ylm(int nylm, Cmp** ylmvec) {
  
  //  IMPORTANT NOTE:
  // For Y_lm with m>=1, we have the real and imaginary parts, 
  // not Y_{l,m} and Y_{l,-m}.  This changes the normalization
  // properties.  We multiply all harmonics with m>=1 by sqrt(2.0),
  // in this routine, so that the volume integral of all terms squared
  // is the same.  

  // This is not sufficient to produce the proper Boundary condition, since
  // the spherical harmonic decomposition of Cmp arrays is not given in the
  // "balanced" representation.  In order to normalize properly, we multiply
  // all fields in get_integrals below by anoth factor of sqrt(2.) when
  // m>=1.

  cout << "Constructing ylm" << endl;

  int nz = mp.get_mg()->get_nzone() ;
  int nr = mp.get_mg()->get_nr(0) ;  
  int np = mp.get_mg()->get_np(0) ;
  int nt = mp.get_mg()->get_nt(0) ;

  for (int l=0 ; l<nz ; l++) {    

    Mtbl Xabs (mp.x) ;
    Mtbl Yabs (mp.y) ;
    Mtbl Zabs (mp.z) ;

    for (int k=0 ; k<np ; k++) {
      for (int j=0 ; j<nt ; j++) {
        for (int i=0 ; i<nr ; i++) {

	  double xval=Xabs(l,k,j,i);
	  double yval=Yabs(l,k,j,i);
	  double zval=Zabs(l,k,j,i);
	  double rval=sqrt(xval*xval+yval*yval+zval*zval);

	  //	  cout <<l<<" "<<k<<" "<<j<<" "<<i<<endl;

	  //l=0,m=0
	  ylmvec[0]->set(l,k,j,i)=1.0*sqrt(1.0/4.0/M_PI);
	  //	  cout << " 0 " << endl;
	  // l=1 included?
	  if (nylm>1 ) {
	    if (nylm <4) {abort();} else {
	      //l=1,m=0
	  ylmvec[1]->set(l,k,j,i)=zval*sqrt(3.0/4.0/M_PI);
	  //l=1,m=1
	  ylmvec[2]->set(l,k,j,i)=-1.0*xval*sqrt(3.0/4.0/M_PI);
	  ylmvec[3]->set(l,k,j,i)=-1.0*yval*sqrt(3.0/4.0/M_PI);
	    }
	  }
	  // l=2 included?
	  if (nylm>4 ) {
	    if (nylm <9) {abort();} else {
	  //l=2,m=0
	  ylmvec[4]->set(l,k,j,i)=(3.0*zval*zval-rval*rval)*sqrt(5.0/16.0/M_PI);
	  //l=2,m=1
	  ylmvec[5]->set(l,k,j,i)=-1.0*zval*xval*sqrt(15.0/4.0/M_PI);
	  ylmvec[6]->set(l,k,j,i)=-1.0*zval*yval*sqrt(15.0/4.0/M_PI);
	  //l=2,m=2
	  ylmvec[7]->set(l,k,j,i)=(xval*xval-yval*yval)*sqrt(15.0/16.0/M_PI);
	  ylmvec[8]->set(l,k,j,i)=2.0*xval*yval*sqrt(15.0/16.0/M_PI);
	    }
	  }
	  // l=3 included?
	  if (nylm>9 ) {
	    if (nylm <16) {abort();} else {
	  //l=3,m=0
	  ylmvec[9]->set(l,k,j,i)=(5.0*pow(zval,3)-3.0*zval*rval*rval)*
	    sqrt(7.0/16.0/M_PI);
	  //l=3,m=1
	  ylmvec[10]->set(l,k,j,i)=-1.0*(5.0*zval*zval-rval*rval)*xval*
	    sqrt(21.0/32.0/M_PI);
	  ylmvec[11]->set(l,k,j,i)=-1.0*(5.0*zval*zval-rval*rval)*yval*
	    sqrt(21.0/32.0/M_PI);
	  //l=3,m=2
	  ylmvec[12]->set(l,k,j,i)=zval*(xval*xval-yval*yval)*
	    sqrt(105./16.0/M_PI);
	  ylmvec[13]->set(l,k,j,i)=zval*(2.0*xval*yval)*
	    sqrt(105./16.0/M_PI);
	  //l=3,m=3
	  ylmvec[14]->set(l,k,j,i)=-1.0*(pow(xval,3)-3.0*xval*yval*yval)*
	    sqrt(35.0/32.0/M_PI);
	  ylmvec[15]->set(l,k,j,i)=-1.0*(3.0*xval*xval*yval-pow(yval,3))*
	    sqrt(35.0/32.0/M_PI);
	    }
	  }
	  // l=4 included?
	  if (nylm>16 ) {
	    if (nylm <25) {abort();} else {
	  //l=4,m=0
	  ylmvec[16]->set(l,k,j,i)=(35.0*pow(zval,4)-30.0*zval*zval*rval*rval+3*pow(rval,4))*
	    sqrt(9.0/256.0/M_PI);
	  //l=4,m=1
	  ylmvec[17]->set(l,k,j,i)=-1.0*(7.0*pow(zval,3)-3*zval*rval*rval)*xval*
	    sqrt(45.0/32.0/M_PI);
	  ylmvec[18]->set(l,k,j,i)=-1.0*(7.0*pow(zval,3)-3*zval*rval*rval)*yval*
	    sqrt(45.0/32.0/M_PI);
	  //l=4,m=2
	  ylmvec[19]->set(l,k,j,i)=(7.0*zval*zval-rval*rval)*(xval*xval-yval*yval)*
	    sqrt(45./64.0/M_PI);
	  ylmvec[20]->set(l,k,j,i)=(7.0*zval*zval-rval*rval)*(2.0*xval*yval)*
	    sqrt(45./64.0/M_PI);
	  //l=4,m=3
	  ylmvec[21]->set(l,k,j,i)=-1.0*zval*(pow(xval,3)-3.0*xval*yval*yval)*
	    sqrt(315.0/32.0/M_PI);
	  ylmvec[22]->set(l,k,j,i)=-1.0*zval*(3.0*xval*xval*yval-pow(yval,3))*
	    sqrt(315.0/32.0/M_PI);
	  //l=4,m=4
	  ylmvec[23]->set(l,k,j,i)=(pow(xval,4)-6*xval*xval*yval*yval+pow(yval,4))*
	    sqrt(315.0/256.0/M_PI);
	  ylmvec[24]->set(l,k,j,i)=4.0*xval*yval*(xval*xval-yval*yval)*
	    sqrt(315.0/256.0/M_PI);
	    }
	  }
	  // l=5 included?
	  if (nylm>25 ) {
	    if (nylm <36) {abort();} else {
	  //l=5,m=0
	  ylmvec[25]->set(l,k,j,i)=(63.0*pow(zval,5)-70.0*pow(zval,3)*rval*rval+15*zval*pow(rval,4))*
	    sqrt(11.0/256.0/M_PI);
	  //l=5,m=1
	  ylmvec[26]->set(l,k,j,i)=-1.0*(21.0*pow(zval,4)-14*zval*zval*rval*rval+pow(rval,4))*xval*
	    sqrt(165.0/256.0/M_PI);
	  ylmvec[27]->set(l,k,j,i)=-1.0*(21.0*pow(zval,4)-14*zval*zval*rval*rval+pow(rval,4))*yval*
	    sqrt(165.0/256.0/M_PI);
	  //l=5,m=2
	  ylmvec[28]->set(l,k,j,i)=(3.0*pow(zval,3)-zval*rval*rval)*(xval*xval-yval*yval)*
	    sqrt(1155./64.0/M_PI);
	  ylmvec[29]->set(l,k,j,i)=(3.0*pow(zval,3)-zval*rval*rval)*(2.0*xval*yval)*
	    sqrt(1155./64.0/M_PI);
	  //l=5,m=3
	  ylmvec[30]->set(l,k,j,i)=-1.0*(9.0*zval*zval-rval*rval)*(pow(xval,3)-3.0*xval*yval*yval)*
	    sqrt(385.0/512.0/M_PI);
	  ylmvec[31]->set(l,k,j,i)=-1.0*(9.0*zval*zval-rval*rval)*(3.0*xval*xval*yval-pow(yval,3))*
	    sqrt(385.0/512.0/M_PI);
	  //l=5,m=4
	  ylmvec[32]->set(l,k,j,i)=zval*(pow(xval,4)-6*xval*xval*yval*yval+pow(yval,4))*
	    sqrt(3465.0/256.0/M_PI);
	  ylmvec[33]->set(l,k,j,i)=zval*4.0*xval*yval*(xval*xval-yval*yval)*
	    sqrt(3465.0/256.0/M_PI);
	  //l=5,m=5
	  ylmvec[34]->set(l,k,j,i)=-1.0*(pow(xval,5)-10.0*pow(xval,3)*yval*yval+5.0*xval*pow(yval,4))*
	    sqrt(693.0/512.0/M_PI);
	  ylmvec[35]->set(l,k,j,i)=-1.0*(5.0*pow(xval,4)*yval-10.0*xval*xval*pow(yval,3)+pow(yval,5))*
	    sqrt(693.0/512.0/M_PI);
	    }
	  }
	  // l=6 included?
	  if (nylm>36 ) {
	    if (nylm <49) {abort();} else {
	  //l=6,m=0
	  ylmvec[36]->set(l,k,j,i)=(231.0*pow(zval,6)-315.0*pow(zval,4)*rval*rval+105.0*zval*zval*pow(rval,4)-5.0*pow(rval,6))*
	    sqrt(13.0/1024.0/M_PI);
	  //l=6,m=1
	  ylmvec[37]->set(l,k,j,i)=-1.0*(33.0*pow(zval,5)-30.0*pow(zval,3)*rval*rval+5.0*zval*pow(rval,4))*xval*
	    sqrt(273.0/256.0/M_PI);
	  ylmvec[38]->set(l,k,j,i)=-1.0*(33.0*pow(zval,5)-30.0*pow(zval,3)*rval*rval+5.0*zval*pow(rval,4))*yval*
	    sqrt(273.0/256.0/M_PI);
	  //l=6,m=2
	  ylmvec[39]->set(l,k,j,i)=(33.0*pow(zval,4)-18.0*zval*zval*rval*rval+pow(rval,4))*(xval*xval-yval*yval)*
	    sqrt(1365./2048.0/M_PI);
	  ylmvec[40]->set(l,k,j,i)=(33.0*pow(zval,4)-18.0*zval*zval*rval*rval+pow(rval,4))*(2.0*xval*yval)*
	    sqrt(1365./2048.0/M_PI);
	  //l=6,m=3
	  ylmvec[41]->set(l,k,j,i)=-1.0*(11.0*pow(zval,3)-3.0*zval*rval*rval)*(pow(xval,3)-3.0*xval*yval*yval)*
	    sqrt(1365.0/512.0/M_PI);
	  ylmvec[42]->set(l,k,j,i)=-1.0*(11.0*pow(zval,3)-3.0*zval*rval*rval)*(3.0*xval*xval*yval-pow(yval,3))*
	    sqrt(1365.0/512.0/M_PI);
	  //l=6,m=4
	  ylmvec[43]->set(l,k,j,i)=(11.0*zval*zval-rval*rval)*(pow(xval,4)-6*xval*xval*yval*yval+pow(yval,4))*
	    sqrt(819.0/1024.0/M_PI);
	  ylmvec[44]->set(l,k,j,i)=(11.0*zval*zval-rval*rval)*4.0*xval*yval*(xval*xval-yval*yval)*
	    sqrt(819.0/1024.0/M_PI);
	  //l=6,m=5
	  ylmvec[45]->set(l,k,j,i)=-1.0*zval*(pow(xval,5)-10.0*pow(xval,3)*yval*yval+5.0*xval*pow(yval,4))*
	    sqrt(9009.0/512.0/M_PI);
	  ylmvec[46]->set(l,k,j,i)=-1.0*zval*(5.0*pow(xval,4)*yval-10.0*xval*xval*pow(yval,3)+pow(yval,5))*
	    sqrt(9009.0/512.0/M_PI);
	  //l=6,m=6
	  ylmvec[47]->set(l,k,j,i)=(pow(xval,6)-15.0*pow(xval,4)*yval*yval+15.0*xval*xval*pow(yval,4)-pow(yval,6))*
	    sqrt(3003.0/2048.0/M_PI);
	  ylmvec[48]->set(l,k,j,i)=(6.0*pow(xval,5)*yval-20.0*pow(xval,3)*pow(yval,3)+6.0*xval*pow(yval,5))*
	    sqrt(3003.0/2048.0/M_PI);
	    }
	  }
	  if(nylm >49) {
	    cout << "l>6 not implemented!!!!!!!"<< endl;
	    abort();
	  }
	}
      }
    }
  }

}

void Et_bin_bhns_extr::get_integrals(int nylm, double* intvec, Cmp** ylmvec,
				     Cmp source) {

  // As mentioned in the comment before get_ylm, our "balanced" representation
  // of the Y_lm Cmp's does not agree with the normalization used to define
  // the spherical harmonic decomposition of Cmp arrays.  Thus, we multiply
  // all terms with m>=1 by another factor of sqrt(2.0) in order to
  // produce the correct decomposition.

  int nz=mp.get_mg()->get_nzone() ;

  Map_af mapping (mp);
  
  const double* a1 = mapping.get_alpha() ;
  const double* b1 = mapping.get_beta() ;
  
  double rlim=a1[nz-1]+b1[nz-1];
  
  int ll=0;
  int mm=0;
  int ncount=0;
  for (int n=0; n<nylm; n++) {

    Cmp ylmsource=(*ylmvec[n]*source);
    int symcheck=1;
    for (int l=0; l<nz; l++) {
      int symc=ylmsource.va.base.get_base_t(l);
      if(symc!=2304 && symc!=1280)symcheck=0;
    }
    if(symcheck==1) {
      intvec[n]=ylmsource.integrale()/(2.0*ll+1.0)/sqrt(2.0*M_PI)/pow(rlim,ll+1);
    } else {
      intvec[n]=0;
    }
    if(mm>=1)intvec[n]*=sqrt(2.0);

    int lnew=0;
    int mnew=0;
    int nnew=0;
    if(mm<ll) {
      if(mm==0) {
	lnew=ll;
	mnew=mm+1;
	nnew=0;
      }
      if(mm>0&&ncount==0) {
	lnew=ll;
	mnew=mm;
	nnew=1;
      }
      if(mm>0&&ncount==1) {
	lnew=ll;
	mnew=mm+1;
	nnew=0;
      }
    }
    if(mm==ll) {
      if(mm==0) {
	lnew=ll+1;
	mnew=0;
	nnew=0;
      }
      if(mm>0&&ncount==0) {
	lnew=ll;
	mnew=mm;
	nnew=1;
      }
      if(mm>0&&ncount==1) {
	lnew=ll+1;
	mnew=0;
	nnew=0;
      }
    }
    ll=lnew;
    mm=mnew;
    ncount=nnew;
  }
}
