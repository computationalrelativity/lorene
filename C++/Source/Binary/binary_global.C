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
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "binary.h"

		    //---------------------------------//
		    //		ADM mass	       //
		    //---------------------------------//

double Binary::mass_adm() const {
    /* 

    if (p_mass_adm == 0x0) {	    // a new computation is requireed
	
	p_mass_adm = new double ; 
	    
#include "unites.h"
    if (this == 0x0) {	// To avoid any compilation warning
      cout << f_unit << msol << km << mevpfm3 << qpig ;
    }
    
    *p_mass_adm = 0 ; 
    
    const Metric&  gamij = (et[0]->get_gamij()) ;
    const Metric& flat = (et[0]->get_flat()) ;
    Map_af map0 (et[0]->get_mp()) ; 
    
    Sym_tensor gamij_cov = gamij.cov() ;
    gamij_cov.std_spectral_base() ;
    Tensor dcov_gamij = gamij_cov.derive_cov(flat) ;
    
    Vector dgamma_1 (map0, CON, map0.get_bvect_spher()) ; 
    dgamma_1 =  contract( flat.con(), 1, contract(flat.con()
       			   , 0, dcov_gamij, 0).scontract(0, 2), 0) ;
    Scalar integrant_1 (dgamma_1(1)) ;
    
    Vector dgamma_2 (map0, CON, map0.get_bvect_spher()) ; 
    dgamma_2 =  contract( flat.con(), 1, contract(flat.con()
			   , 0, dcov_gamij, 1).scontract(0, 2), 0) ;
    Scalar integrant_2 (dgamma_2(1)) ;

    *p_mass_adm = map0.integrale_surface_infini (integrant_1 - integrant_2)/(4*qpig) ;

		
    }	// End of the case where a new computation was necessary
    */
    return *p_mass_adm ; 
    
}


		    //---------------------------------//
		    //		Komar mass	       //
		    //---------------------------------//

double Binary::mass_kom() const {
    /*
  if (p_mass_kom == 0x0) {	    // a new computation is requireed
    
    p_mass_kom = new double ; 
      
#include "unites.h"
    if (this == 0x0) {	// To avoid any compilation warning
      cout << f_unit << msol << km << mevpfm3 ;
    }
    
    *p_mass_kom = 0 ; 
    
    const Tensor& logn = et[0]->get_logn() ;
    const Metric&  gamij = (et[0]->get_gamij()) ;
    Map_af map0 (et[0]->get_mp()) ; 
    
    Vector vect (map0, CON, map0.get_bvect_spher()) ;  
    vect = logn.derive_con(gamij) ;
    Scalar integrant (vect(1)) ;
    
    *p_mass_kom = map0.integrale_surface_infini (integrant) / qpig ;
    
  }	// End of the case where a new computation was necessary
    */
  return *p_mass_kom ; 
    
}


		    //---------------------------------//
		    //	 Total angular momentum        //
		    //---------------------------------//

const Tbl& Binary::angu_mom() const {
    /*
    if (p_angu_mom == 0x0) {	    // a new computation is requireed
	
      p_angu_mom = new Tbl(3) ; 
      
      p_angu_mom->annule_hard() ;	// fills the double array with zeros
      
      
      Map_af map0 (et[0]->get_mp()) ; 
      const Sym_tensor& kij_auto = et[0]->get_tkij_auto() ;
      const Sym_tensor& kij_comp = et[0]->get_tkij_auto() ;
      
      Sym_tensor kij = kij_auto + kij_comp ;
      kij.change_triad(map0.get_bvect_cart()) ;
 
      // X component
      // -----------

      Vector vect_x(map0, CON, map0.get_bvect_cart()) ;
       
      Scalar kij_mod1(map0) ;
      Scalar kij_mod2(map0) ;
     
      for (int i=1; i<=3; i++) {

	  kij_mod1 = kij(3, i) ;
	  kij_mod2 = kij(2, i) ;
	  	  
	  kij_mod1.mult_rsint() ;
	  kij_mod1.get_spectral_va().mult_sp() ;
	  
	  kij_mod2.mult_r() ;
	  kij_mod2.get_spectral_va().mult_cp() ;

	  vect_x.set(i) = kij_mod1 - kij_mod2 ;
      
      }
      vect_x.std_spectral_base() ;
      vect_x.change_triad(map0.get_bvect_spher()) ;
      Scalar integrant_x (vect_x(0)) ;
      
      p_angu_mom->set(0) = map0.integrale_surface_infini (integrant_x) 
	                  / (8*M_PI) ;
      
      // Y component
      // -----------
      
      Vector vect_y(map0, CON, map0.get_bvect_cart()) ;
   
      for (int i=1; i<=3; i++) {

	  kij_mod1 = kij(1, i) ;
	  kij_mod2 = kij(3, i) ;	  
	  
	  kij_mod1.mult_r() ;
	  kij_mod1.get_spectral_va().mult_ct() ;
	  
	  kij_mod2.mult_rsint() ;
	  kij_mod2.get_spectral_va().mult_cp() ;
	 
	  vect_y.set(i) = kij_mod1 - kij_mod2 ;
 
      }
      vect_y.std_spectral_base() ;
      vect_y.change_triad(map0.get_bvect_spher()) ;
      Scalar integrant_y (vect_y(0)) ;
      
      p_angu_mom->set(1) = map0.integrale_surface_infini (integrant_y) 
	                  / (8*M_PI) ;
      
      // Z component
      // -----------
      
      Vector vect_z(map0, CON, map0.get_bvect_cart()) ;
  
      for (int i=0; i<=2; i++) {

	  kij_mod1 = kij(2, i) ;
	  kij_mod2 = kij(1, i) ;
	  	  
	  kij_mod1.mult_rsint() ;
	  kij_mod1.get_spectral_va().mult_cp() ;

	  kij_mod2.mult_rsint() ;
	  kij_mod2.get_spectral_va().mult_sp() ;
	  
	vect_z.set(i) = kij_mod1 - kij_mod2 ;
      }
       
      vect_z.std_spectral_base() ;
      vect_z.change_triad(map0.get_bvect_spher()) ;
      Scalar integrant_z (vect_z(0)) ;
      
      p_angu_mom->set(2) = map0.integrale_surface_infini (integrant_z) 
	                 / (8*M_PI) ;
      
      (*p_angu_mom).change_triad(map0.get_bvect_spher()) ;
      
    }	// End of the case where a new computation was necessary
    */
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

double Binary::virial() const {
    /*
    if (p_virial == 0x0) {	    // a new computation is requireed
	
	p_virial = new double ; 
	    
	    *p_virial = 1. - mass_kom() / mass_adm() ; 
	    
	}
    */
    return *p_virial ; 
    
}


/*	     //----------------------------------------------//
	     //	 Virial error by Gourgoulhon and Bonazzola   //
	     //----------------------------------------------//

double Bin_ns_ncp::virial_gb() const {
    
    if (p_virial_gb == 0x0) {	    // a new computation is requireed
	
	p_virial_gb = new double ; 
	    
	if (star1.is_relativistic()) {	// Relativistic case
					// -----------------
	
#include "unites.h"
    if (this == 0x0) {	// To avoid any compilation warning
	cout << f_unit << msol << km << mevpfm3 ;
    }

	    assert( star2.is_relativistic() ) ; 
	    
	    *p_virial_gb = 0 ;

	    double vir_pres = 0. ;
	    double vir_extr = 0. ;
	    double vir_grav = 0. ;

	    for (int i=0; i<=1; i++) {  // loop on the stars

	      const Cmp& a2 = (et[i]->get_a_car())() ;
	      const Cmp& se = (et[i]->get_s_euler())() ;
	      const Cmp& ak2_auto = (et[i]->get_akcar_auto())() ;
	      const Cmp& ak2_comp = (et[i]->get_akcar_comp())() ;

	      const Tenseur& dnu_auto = et[i]->get_d_logn_auto() ;
	      const Tenseur& dnu_comp = et[i]->get_d_logn_comp() ;
	      const Tenseur& dbe_auto = et[i]->get_d_beta_auto() ;
	      const Tenseur& dbe_comp = et[i]->get_d_beta_comp() ;

	      Cmp source = 2. * a2 * sqrt(a2) * se ;
	      vir_pres += source.integrale() ;

	      source = 1.5 * sqrt(a2) * (ak2_auto + ak2_comp) / qpig ;
	      vir_extr += source.integrale() ;

	      Tenseur sprod1 = flat_scalar_prod(dbe_auto, dbe_auto+dbe_comp) ;
	      Tenseur sprod2 = flat_scalar_prod(dnu_auto, dnu_auto+dnu_comp) ;
	      Tenseur sprod3 = flat_scalar_prod(dbe_auto, dnu_auto+dnu_comp) ;

	      source = sqrt(a2) * ( sprod1() - sprod2() - 2.*sprod3() )/qpig ;
	      vir_grav += source.integrale() ;

	    }  // End of the loop on the stars


	    *p_virial_gb = (vir_pres + vir_extr + vir_grav) / mass_adm() ;
	    
	}
	else {		// Newtonian case 
			// --------------
			
	    *p_virial_gb = virial() ; 
	    
		
	}   // End of the Newtonian case 

    }	// End of the case where a new computation was necessary
    
    return *p_virial_gb ; 
    
}

*/
