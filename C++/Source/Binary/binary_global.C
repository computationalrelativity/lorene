/*
 * Methods of class Binary to compute global quantities
 *
 * (see file binary.h for documentation)
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


char binary_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2004/07/21 11:46:24  f_limousin
 * Add function mass_adm_vol() to compute the ADM mass of the system
 * with a volume integral instead of a surface one.
 *
 * Revision 1.7  2004/05/25 14:25:53  f_limousin
 * Add the virial theorem for conformally flat configurations.
 *
 * Revision 1.6  2004/03/31 12:44:54  f_limousin
 * Minor modifs.
 *
 * Revision 1.5  2004/03/25 10:29:01  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2004/02/27 10:25:30  f_limousin
 * Modif. to avoid an error in compilation.
 *
 * Revision 1.3  2004/02/27 10:03:04  f_limousin
 * The computation of mass_adm() and mass_komar() is now OK !
 *
 * Revision 1.2  2004/01/20 15:21:36  f_limousin
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
#include "map.h"
#include "unites.h"
#include "coord.h"

		    //---------------------------------//
		    //		ADM mass	       //
		    //---------------------------------//

double Binary::mass_adm() const {
    
  using namespace Unites ;
  if (p_mass_adm == 0x0) {	    // a new computation is requireed
    
    p_mass_adm = new double ; 
	    
    *p_mass_adm = 0 ; 
    
    const Metric& gamij = (et[0]->get_gamma()) ;
    const Metric& flat = (et[0]->get_flat()) ;
    const Map_af map0 (et[0]->get_mp()) ;
    
    Sym_tensor gamij_cov = gamij.cov() ;
    const Tensor& dcov_gamij = gamij_cov.derive_cov(flat) ;
    
    Vector dgamma_1 =  contract( flat.con(), 1, contract(contract(flat.con()
       			   , 0, dcov_gamij, 2), 0, 2), 0) ;
    Scalar integrant_1 (dgamma_1(1)) ;
    
    Vector dgamma_2 =  contract( flat.con(), 1, contract(contract(flat.con()
			   , 0, dcov_gamij, 0), 0, 1), 0) ;
    Scalar integrant_2 (dgamma_2(1)) ;

    *p_mass_adm = map0.integrale_surface_infini (integrant_1 - integrant_2)/(4*qpig) ;

		
    }	// End of the case where a new computation was necessary
    
    return *p_mass_adm ; 
    
}

double Binary::mass_adm_vol() const {

  using namespace Unites ;

  double mass_adm ;
  mass_adm = 0. ;

  for (int i=0; i<=1; i++) {	    // loop on the stars

      const Scalar& psi4 = et[i]->get_psi4() ;
      const Scalar& ener_euler = et[i]->get_ener_euler() ;
      const Scalar& kcar_auto = et[i]->get_kcar_auto() ;
      const Scalar& kcar_comp = et[i]->get_kcar_comp() ;
      const Metric& gtilde = et[i]->get_gtilde() ;
      const Metric& flat = et[i]->get_flat() ;
      const Sym_tensor& hij_auto = et[i]->get_hij_auto() ;
      const Scalar& beta = et[i]->get_beta_auto() + et[i]->get_beta_comp() ;
      const Scalar& logn = et[i]->get_logn_auto() + et[i]->get_logn_comp() ;

      const Sym_tensor& gtilde_cov = gtilde.cov() ;
      const Sym_tensor& gtilde_con = gtilde.con() ;
      const Scalar& lnpsi = 0.5 * (beta - logn) ;
      Scalar nnn = exp(logn) ;
      nnn.std_spectral_base() ;
      Scalar qq = exp(beta) ;
      qq.std_spectral_base() ;

      const Tensor& dcov_hij_auto = hij_auto.derive_cov(flat) ;
      const Tensor& dcov_gtilde = gtilde_cov.derive_cov(flat) ;
      const Tensor& dcov_lnpsi = lnpsi.derive_cov(flat) ;
      const Tensor& dcov_logn = logn.derive_cov(flat) ;
      Tensor dcovdcov_qq = qq.derive_cov(flat).derive_cov(flat) ;
      dcovdcov_qq.inc_dzpuis() ;
      Tensor dcovdcov_nnn = nnn.derive_cov(flat).derive_cov(flat) ;
      dcovdcov_nnn.inc_dzpuis() ;
      
      
      Scalar source = pow(psi4, 1.25) * (ener_euler + (kcar_auto + kcar_comp)
					 / (4.*qpig)) ;

      source += - pow(psi4, 0.25) / (4.*qpig) * 
	  (0.0625 * contract(gtilde_con, 0, 1, contract(
			      dcov_hij_auto, 0, 1, dcov_gtilde, 0, 1), 0, 1) - 
           0.125 * contract(gtilde_con, 0, 1, contract(dcov_hij_auto, 
			      0, 1, dcov_gtilde, 0, 2), 0, 1) + 
	   contract(hij_auto, 0, 1, 2*dcov_lnpsi * dcov_lnpsi +
		    4*dcov_lnpsi * dcov_logn - dcovdcov_qq / qq +
		    dcovdcov_nnn / nnn, 0, 1) ) ;

      source.std_spectral_base() ;

      mass_adm += source.integrale() ;
	  
  }

  return mass_adm ;

}

		    //---------------------------------//
		    //		Komar mass	       //
		    //---------------------------------//

double Binary::mass_kom() const {
    
  using namespace Unites ;

  if (p_mass_kom == 0x0) {	    // a new computation is requireed
    
    p_mass_kom = new double ; 
      
    *p_mass_kom = 0 ; 
    
    const Tensor& logn = et[0]->get_logn() ;
    const Metric&  gamij = (et[0]->get_gamma()) ;
    Map_af map0 (et[0]->get_mp()) ; 
    
    Vector vect = logn.derive_con(gamij) ;
    Scalar integrant (vect(1)) ;
    
    *p_mass_kom = map0.integrale_surface_infini (integrant) / qpig ;
    
  }	// End of the case where a new computation was necessary
    
  return *p_mass_kom ; 
    
}


		    //---------------------------------//
		    //	 Total angular momentum        //
		    //---------------------------------//

const Tbl& Binary::angu_mom() const {
    
    if (p_angu_mom == 0x0) {	    // a new computation is requireed
	
      p_angu_mom = new Tbl(3) ; 
      
      p_angu_mom->annule_hard() ;	// fills the double array with zeros
      
      const Sym_tensor& kij_auto = et[0]->get_tkij_auto() ;
      const Sym_tensor& kij_comp = et[0]->get_tkij_comp() ;
      const Tensor& psi4 = et[0]->get_psi4() ;
      const Map_af map0 (kij_auto.get_mp()) ;

      Sym_tensor kij = (kij_auto + kij_comp) / psi4 ;
      kij.change_triad(map0.get_bvect_cart()) ;
 
      // X component
      // -----------

      Vector vect_x(et[0]->get_mp(), CON, map0.get_bvect_cart()) ;      
       
      for (int i=1; i<=3; i++) {

	  Scalar kij_1 = kij(3, i) ;
	  Scalar kij_2 = kij(2, i) ;
	  	  
	  kij_1.mult_rsint() ;
	  Valeur vtmp = kij_1.get_spectral_va().mult_sp() ;
	  kij_1.set_spectral_va() = vtmp ; 
	  
	  kij_2.mult_r() ;
	  vtmp = kij_2.get_spectral_va().mult_ct() ;
	  kij_2.set_spectral_va() = vtmp ; 

	  vect_x.set(i) = kij_1 - kij_2 ;
      }
 
      vect_x.change_triad(map0.get_bvect_spher()) ;
      Scalar integrant_x (vect_x(1)) ;
      
      p_angu_mom->set(0) = map0.integrale_surface_infini (integrant_x) 
	                  / (8*M_PI) ;
      
      // Y component
      // -----------
      
      Vector vect_y(et[0]->get_mp(), CON, map0.get_bvect_cart()) ;      
   
      for (int i=1; i<=3; i++) {

	  Scalar kij_1 = kij(1, i) ;
	  Scalar kij_2 = kij(3, i) ;	  
	  
	  kij_1.mult_r() ;
	  Valeur vtmp = kij_1.get_spectral_va().mult_ct() ;
	  kij_1.set_spectral_va() = vtmp ; 
	  
	  kij_2.mult_rsint() ;
	  vtmp = kij_2.get_spectral_va().mult_cp() ;
	  kij_2.set_spectral_va() = vtmp ; 
	 
	  vect_y.set(i) = kij_1 - kij_2 ;
      }

      vect_y.change_triad(map0.get_bvect_spher()) ;
      Scalar integrant_y (vect_y(1)) ;
      
      p_angu_mom->set(1) = map0.integrale_surface_infini (integrant_y) 
	                  / (8*M_PI) ;
      
      // Z component
      // -----------

      Vector vect_z(et[0]->get_mp(), CON, map0.get_bvect_cart()) ;      

      for (int i=1; i<=3; i++) {

	  Scalar kij_1 = kij(2, i) ;
	  Scalar kij_2 = kij(1, i) ;	  	  

	  kij_1.mult_rsint() ;
	  Valeur vtmp = kij_1.get_spectral_va().mult_cp() ;
	  kij_1.set_spectral_va() = vtmp ; 

	  kij_2.mult_rsint() ;
	  vtmp =  kij_2.get_spectral_va().mult_sp() ;
	  kij_2.set_spectral_va() = vtmp ;

	  vect_z.set(i) = kij_1 - kij_2 ;
      }
       
      vect_z.change_triad(map0.get_bvect_spher()) ;
      Scalar integrant_z (vect_z(1)) ;
      
      p_angu_mom->set(2) = map0.integrale_surface_infini (integrant_z) 
	                 ;// (8*M_PI) ;
      
      
    }	// End of the case where a new computation was necessary
    
    return *p_angu_mom ; 
    
}



		    //---------------------------------//
		    //		Total energy	       //
		    //---------------------------------//

double Binary::total_ener() const {
    /*
    if (p_total_ener == 0x0) {	    // a new computation is requireed
	
	p_total_ener = new double ; 
	    
	    *p_total_ener = mass_adm() - star1.mass_b() - star2.mass_b() ; 
	    
    }	// End of the case where a new computation was necessary
    
    */
    return *p_total_ener ; 
    
}


		    //---------------------------------//
		    //	 Error on the virial theorem   //
		    //---------------------------------//

double Binary::virial_cp() const {
    
    if (p_virial == 0x0) {	    // a new computation is requireed
	
	p_virial = new double ; 
	    
	    *p_virial = 1. - mass_kom() / mass_adm() ; 
	    
	}
    
    return *p_virial ; 
    
}
