/*
 *  Resolution of the divergence-free vector Poisson equation
 *
 *    (see file vector.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

char vector_df_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2004/12/23 10:23:06  j_novak
 * Improved method 5 in the case when some components are zero.
 * Changed Vector_divfree::poisson to deduce eta from the equation. 
 *
 * Revision 1.7  2004/03/31 12:33:21  f_limousin
 * Minor modifs.
 *
 * Revision 1.6  2004/03/29 16:24:26  j_novak
 * *** empty log message ***
 *
 * Revision 1.5  2004/03/29 16:11:24  j_novak
 * Minor modifs. not to have compilation warnings.
 *
 * Revision 1.4  2004/03/11 13:14:09  f_limousin
 * Add method Vector_divfree::poisson with parameters.
 *
 * Revision 1.3  2003/10/27 09:03:01  j_novak
 * Added an assert on the triad
 *
 * Revision 1.2  2003/10/20 19:44:43  e_gourgoulhon
 * First full version
 *
 * Revision 1.1  2003/10/20 14:45:27  e_gourgoulhon
 * Not ready yet...
 *
 *
 * $Header$
 *
 */

// C headers
//#include <>

// Lorene headers
#include "param_elliptic.h"
#include "param.h"
#include "cmp.h"


Vector_divfree Vector_divfree::poisson() const {

  // All this has a meaning only for spherical components:
#ifndef NDEBUG 
  const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
  assert(bvs != 0x0) ; 
#endif
  // Solution for the r-component
  // ----------------------------
	
  Scalar f_r(*mp) ;
  if (cmp[0]->get_etat() != ETATZERO) {
      int nz = mp->get_mg()->get_nzone() ;

      //------------------------------------------
      // The elliptic operator for the r-component
      //------------------------------------------
      
      Param_elliptic param_fr(*cmp[0]) ;
      for (int lz=0; lz<nz; lz++) 
	  param_fr.set_poisson_vect_r(lz) ;
      
      f_r = cmp[0]->sol_elliptic(param_fr) ;
  }
  else
      f_r.set_etat_zero() ;
  
  // Deduction of eta from the div=0 cond. + the Poisson eq.
  //--------------------------------------------------------
  Scalar source_eta = - *(cmp[0]) ;
  source_eta.mult_r_dzpuis(2) ;
  f_r.set_spectral_va().ylm() ;
  Scalar tmp = 2*f_r + f_r.lapang() ;
  tmp.div_r_dzpuis(2) ;
  source_eta += tmp ;
  source_eta = source_eta.primr() ;
  
  Scalar eta_resu = (f_r + source_eta).poisson_angu() ;
  
  // Solution for mu
  // ---------------
	
  Scalar mu_resu = mu().poisson() ;
	
  // Final result
  // ------------
	
  Vector_divfree resu(*mp, *triad, *met_div) ; 

  resu.set_vr_eta_mu(f_r, eta_resu, mu_resu) ; 
	
	return resu ;

}


Vector_divfree Vector_divfree::poisson(Param& par ) const {

   // All this has a meaning only for spherical components:
#ifndef NDEBUG 
    const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
    assert(bvs != 0x0) ; 
#endif
    
    int nitermax = par.get_int() ; 
    int& niter = par.get_int_mod() ; 
    double relax = par.get_double() ; 
    double precis = par.get_double(1) ;     
    Cmp& ss_khi = par.get_cmp_mod(0) ;
    Cmp& ss_mu = par.get_cmp_mod(1) ;
      
    // Solution for the r-component
    // ----------------------------
    
    Scalar source_r = *(cmp[0]) ; 
    source_r.mult_r() ; 
    
    Param par_khi ;
    par_khi.add_int(nitermax, 0) ;
    par_khi.add_int_mod(niter, 0) ;
    par_khi.add_double(relax, 0) ;
    par_khi.add_double(precis, 1) ;
    par_khi.add_cmp_mod(ss_khi, 0) ;

    Scalar khi (*mp) ;
    khi.set_etat_zero() ;

    source_r.poisson(par_khi, khi) ; 
    khi.div_r() ;   // khi now contains V^r
    
    // Solution for mu
    // ---------------
    
    Param par_mu ;
    par_mu.add_int(nitermax, 0) ;
    par_mu.add_int_mod(niter, 0) ;
    par_mu.add_double(relax, 0) ;
    par_mu.add_double(precis, 1) ;
    par_mu.add_cmp_mod(ss_mu, 0) ;

    Scalar mu_resu (*mp) ;
    mu_resu.set_etat_zero() ;
 
    mu().poisson(par_mu, mu_resu) ;

    // Final result
    // ------------
    
    Vector_divfree resu(*mp, *triad, *met_div) ; 
    
    resu.set_vr_mu(khi, mu_resu) ; 
    
    return resu ;

}

