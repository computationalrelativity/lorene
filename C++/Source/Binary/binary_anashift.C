/*
 * Method of class Binary to set some analytical form to the shift vector.
 *
 * (see file binary.h for documentation).
 */

/*
 *   Copyright (c) 2004 Francois Limousin
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char binary_anashift_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2005/02/17 17:34:50  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.6  2004/03/25 10:29:01  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.5  2004/02/27 10:01:32  f_limousin
 * Correct sign of shift_auto to agree with the new convention
 * for shift.
 *
 * Revision 1.4  2004/01/22 10:09:41  f_limousin
 * First executable version
 *
 * Revision 1.3  2004/01/20 15:21:23  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "binary.h"
#include "unites.h"

void Binary::analytical_shift(){
    
  using namespace Unites ;
        
    for (int i=0; i<2; i++) {

	// Radius of the star:
	double a0 = et[i]->ray_eq() ; 
    
	// Mass ratio
	double p_mass = et[i]->mass_g() / et[1-i]->mass_g() ; 
    
	// G M Omega R / (1+p) 
	double www = ggrav * et[i]->mass_g() * omega 
		    * separation() / (1. + p_mass) ;  
    
	const Map& mp = et[i]->get_mp() ; 
	Scalar tmp(mp) ;  
	Scalar tmp_ext(mp) ;  
	int nzet = et[i]->get_nzet() ; 
	int nzm1 = mp.get_mg()->get_nzone() - 1 ; 
    
	Vector w_shift (mp, CON, mp.get_bvect_cart()) ;
	Scalar khi_shift (mp) ;

	// Computation of w_shift 
	// ----------------------
	// X component
	// -----------

	w_shift.set(1) = 0 ; 

	// Y component
	// -----------

        // For the incompressible case :
	tmp = - 6  * www / a0 * ( 1 - (mp.r)*(mp.r) / (3*a0*a0) ) ; 

	tmp.annule(nzet, nzm1) ; 
	tmp_ext = - 4 * www / mp.r ;
	tmp_ext.annule(0, nzet-1) ; 
    
	w_shift.set(2) = tmp + tmp_ext ; 

	// Z component
	// -----------
	w_shift.set(3) = 0 ; 

	w_shift.std_spectral_base() ; 
	    
	Scalar mp_x (mp) ;
	Scalar mp_y (mp) ;
	Scalar mp_z (mp) ;
	mp_x = mp.x ;
	mp_y = mp.y ;
	mp_z = mp.z ;

	Scalar Wjxj = w_shift(1) * mp_x + w_shift(2) * mp_y + 
	             w_shift(3) * mp_z ;

	int nzone = mp.get_mg()->get_nzone() ;
	int nr = mp.get_mg()->get_nr(0);
	int nt = mp.get_mg()->get_nt(0);
	int np = mp.get_mg()->get_np(0);
	
	// In the last domain, Wjwj does not depend on r but it 
	// has a value nan at infinity : we put at infinity his value for 
	// another r.

	for (int k=0; k<=np-1; k++)
	    for (int j=0; j<=nt-1; j++){
		Wjxj.set_grid_point(nzone-1, k, j, nr-1) = 
		                  Wjxj.val_grid_point(nzone-1, k, j, nr-2) ; 
	    }
	

	Wjxj.std_spectral_base() ;

	// Computation of khi_shift
	// ------------------------

	tmp = 2 * www / a0 * (mp.y) * ( 1 - 3 * (mp.r)*(mp.r) / (5*a0*a0) ) ;
	tmp.annule(nzet, nzm1) ; 
	tmp_ext = 0.8 * www * a0*a0 * (mp.sint) * (mp.sinp) 
					    / (mp.r * mp.r) ;   
	tmp_ext.annule(0, nzet-1) ; 

	khi_shift = tmp + tmp_ext ; 

	// Sets the standard spectral bases for a scalar field
	khi_shift.std_spectral_base() ; 	    
    

	// Computation of shift auto.
	// --------------------------
	
	const Metric_flat& flat (mp.flat_met_cart()) ;
	Vector temp(mp, CON, mp.get_bvect_cart()) ;
 
	temp = khi_shift.derive_con(flat) + Wjxj.derive_con(flat) ;
	temp.dec_dzpuis(2) ;

	// See Eq (92) from Gourgoulhon et al.(2001) and with the new 
	// convention for shift = - N^i
	
	et[i]->set_beta_auto() = - 7./8. * w_shift + 1./8. * temp ;
	 
	et[i]->set_beta_auto().std_spectral_base() ;
	et[i]->set_beta_auto().change_triad(mp.get_bvect_spher()) ;

    }

}
