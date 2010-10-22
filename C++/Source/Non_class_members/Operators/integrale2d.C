/*
 *  Performs a 2D integration from 0 to infty of the l=0 part of a Scalar
 *
 */

/*
 *   Copyright (c) 2010 Jerome Novak
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

char integrale2d_C[] = "$Header$" ;

/*
 * $Id $
 * $Log $
 *
 * $Header$
 *
 */

// Lorene headers
#include "tensor.h"

double integrale2d(const Scalar& fi) {

  const Map* mp = &fi.get_mp() ;
  const Map_af* mp_aff =  dynamic_cast<const Map_af*>(mp) ;
  assert( mp_aff != 0x0 ) ;
  assert( fi.check_dzpuis(4) || fi.check_dzpuis(3) )  ;

  double lambda ;
  if (fi.get_etat() == ETATZERO) lambda = 0. ;
  else {
    const Base_val& base = fi.get_spectral_base() ;
    Scalar tmp(*mp_aff) ;
    tmp.annule_hard() ;
    tmp.set_dzpuis(fi.get_dzpuis()) ;
    tmp += fi ;
    tmp.set_spectral_base(base) ;
    const Mg3d& mg = *mp_aff->get_mg() ;
    int nz0 = mg.get_nzone() ;
    tmp.mult_r_dzpuis(2) ;
    tmp.set_spectral_va().coef() ;
    tmp.set_spectral_va().ylm_i() ;
    Mtbl_cf& mtc = *tmp.set_spectral_va().c_cf ;
    for (int lz=0; lz<nz0; lz++) {
      int np0 = mg.get_np(lz)+2 ;
      int nt0 = mg.get_nt(lz) ;
      int nr0 = mg.get_nr(lz) ;
      for (int k=0; k<np0; k++)
	for (int j=0; j<nt0; j++)
	  for (int i=0; i<nr0; i++) {
	    double resu = 0. ; 
	    if ((j==0)&&(k==0))
	      resu = mtc(lz, k, j, i) ;
	    mtc.set(lz, k, j, i) = resu ;
	  }
    }
    if (tmp.get_spectral_va().c != 0x0) {
      delete tmp.set_spectral_va().c ;
      tmp.set_spectral_va().c = 0x0 ;
	}

    Scalar integ(*mp_aff) ;
    mp_aff->primr(tmp, integ, true) ;
    lambda = -integ.val_grid_point(0,0,0,0) ;
  }
    
  return lambda ;

}
