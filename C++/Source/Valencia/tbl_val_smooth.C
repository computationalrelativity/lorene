/*
 *  Smoothes the junction with an eventual atmosphere.
 *
 */

/*
 *   Copyright (c) 2004 Jerome Novak
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

char tbl_val_smooth_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/12/30 16:14:01  j_novak
 * Changed the name of a shadowed variable.
 *
 * Revision 1.2  2004/12/03 13:24:01  j_novak
 * Minor modif.
 *
 * Revision 1.1  2004/11/26 17:02:19  j_novak
 * Added a function giving a smooth transition to the atmosphere.
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "tbl_val.h"

//Local prototypes
void radial_smoothing(double* , const double* , int , double) ;
//****************************************************************************

void Tbl_val::smooth_atmosphere(double atmosphere_thr) {

  const Gval_spher* gspher = dynamic_cast<const Gval_spher*>(gval) ;
  assert(gspher != 0x0) ;
  int ndim = gspher->get_ndim() ;
  int fant = gspher->get_fantome() ;
  int nr = get_dim(0) + 2*fant;

  switch (ndim) {
  case 1: {
    radial_smoothing(t, gspher->zr->t, nr, atmosphere_thr) ;
    break ;
  } 
  case 2: {
    int nt = get_dim(1)  + 2*fant ;
    for (int j=0; j<nt; j++) 
      radial_smoothing(t+j*nr, gspher->zr->t, nr, atmosphere_thr) ;
    break ;
  }
  case 3: {
    int nt = get_dim(1)  + 2*fant ;
    int np = get_dim(2)  + 2*fant ;
    for (int j=0; j<nt; j++) 
      for (int k=0; k<np; k++) 
	radial_smoothing(t+k*nt*nr+j*nr, gspher->zr->t, nr, atmosphere_thr) ;
    break ;
  }
    
  default: {
    cerr << "Tbl_val::smooth_atmosphere : strange number of dimensions!" 
	 << endl ;
    abort() ;
    break ;
  }
  }
  return ;
}

void radial_smoothing(double* tab, const double* rr, int n, double rho) {

  assert((tab!= 0x0)&&(rr!=0x0)) ;
  assert (rho >= 0.) ;
  
  if (fabs(tab[n-1]) > rho) // no atmosphere here
    return ;

  double* t = tab + (n-1) ;
  int indice = -1 ;
  bool atmos = true ;
  bool jump = false ;
  for (int i=0; ((i<n)&&(atmos)); i++) {
    if (atmos) atmos = ( fabs(*t) < rho) ;
    t-- ;
    if (atmos) {
      jump = ( fabs(*t) > rho ) ;
      if (jump) // discontinuity found
	indice = n - i - 2 ;
    }
  }
  if (indice == -1) return ;
  int np = 2*(n-indice-2)/3 ;
  int nm = indice / 100 + 3 ;
  assert(n > nm+np) ;
  if (indice < n - np+1) { // enough points to interpolate

    // The inteprolation is done using a cubic polynomial
    //---------------------------------------------------

    int ileft = indice - nm + 2 ;
    int iright = indice + np - 1 ;
    double alpha = ( rr[ileft - 2] - rr[ileft - 1]) /
      ( rr[ileft -1] - rr[ileft]) ;
    double der_l = ( alpha*(alpha+2.)*tab[ileft] 
		     - (1.+alpha)*(1.+alpha)*tab[ileft-1] 
		     + tab[ileft-2] ) /
      ( (1.+alpha)*(rr[ileft - 1] - rr[ileft - 2]) ) ;
    double f_l = tab[ileft] ;
    double f_r = tab[iright] ;
    double tau = rr[ileft] - rr[iright] ;
    double alp = der_l / (tau*tau) + 2.*(f_r - f_l)/(tau*tau*tau) ;
    double bet = 0.5*(der_l + alp*tau*(rr[iright] - 3*rr[ileft])) / tau ;
    for (int i=ileft; i<iright; i++) {
      tab[i] = f_r + (alp*rr[i]+bet)*(rr[i] - rr[iright])*(rr[i] - rr[iright]);
    }
  }
  else { // too few points to interpolate -> linear extrapolation ...
    int ileft = indice ;
    double alpha = ( rr[ileft - 2] - rr[ileft - 1]) /
      ( rr[ileft -1] - rr[ileft]) ;
    double der_l = ( alpha*(alpha+2.)*tab[ileft] 
		     - (1.+alpha)*(1.+alpha)*tab[ileft-1] 
		     + tab[ileft-2] ) /
      ( (1.+alpha)*(rr[ileft - 1] - rr[ileft - 2]) ) ;
    for (int i=ileft; i<n; i++) {
      tab[i] = tab[ileft] + (rr[i] - rr[ileft])*der_l ;
    }    
  }
  return ;
}
