/*
 * Methods of class Bin_ns_ncp to compute global quantities
 *
 * (see file bin_ns_ncp.h for documentation)
 */

/*
 *   Copyright (c) 2003 Francois Limousin
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


char bin_ns_ncp_global_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "bin_ns_ncp.h"

		    //---------------------------------//
		    //		ADM mass	       //
		    //---------------------------------//

double Bin_ns_ncp::mass_adm() const {
    
    if (p_mass_adm == 0x0) {	    // a new computation is requireed
	
	p_mass_adm = new double ; 
	    
#include "unites.h"
    if (this == 0x0) {	// To avoid any compilation warning
      cout << f_unit << msol << km << mevpfm3 << qpig ;
    }
	    *p_mass_adm = 0 ; 
	    
	    for(int i=0; i<=1; i++) {

	      const Metrique& gamma_auto = (et[i]->get_metgamma_auto()) ;
	      const Metrique& flat = (et[i]->get_flat()) ;
	      Map_af mapping (et[i]->get_mp()) ; 
	    
	      Tenseur dgamma_1 (mapping, 1, COV, mapping.get_bvect_cart()) ; 
	      dgamma_1.set_etat_qcq() ;
	      dgamma_1 = contract(contract(contract(operator*(flat.con(), operator*(flat.con(), (gamma_auto.con()).derive_cov(flat))), 3, 6), 2, 3), 1, 2) ;
	      dgamma_1.change_triad(mapping.get_bvect_spher()) ;	    
	      Cmp integrant_1 (dgamma_1(0)) ;
	    

	      Tenseur dgamma_2 (mapping, 1, COV, mapping.get_bvect_cart()) ; 
	      dgamma_2.set_etat_qcq() ;
	      dgamma_2 = contract(contract(contract(operator*(flat.con(), operator*(flat.con(), (gamma_auto.con()).derive_cov(flat))), 3, 6), 2, 4), 1, 2) ;
	      dgamma_2.change_triad(mapping.get_bvect_spher()) ;
	      Cmp integrant_2 (dgamma_2(0)) ;

	      *p_mass_adm += mapping.integrale_surface_infini (integrant_1 - integrant_2)/(4*qpig) ;

	    }
		
    }	// End of the case where a new computation was necessary
    
    return *p_mass_adm ; 
    
}


		    //---------------------------------//
		    //		Komar mass	       //
		    //---------------------------------//

double Bin_ns_ncp::mass_kom() const {
    
    if (p_mass_kom == 0x0) {	    // a new computation is requireed
	
	p_mass_kom = new double ; 
	    
#include "unites.h"
    if (this == 0x0) {	// To avoid any compilation warning
	cout << f_unit << msol << km << mevpfm3 ;
    }
	
	    *p_mass_kom = 0 ; 

	    for(int i=0; i<=1; i++) {
	    
	      const Tenseur& logn_auto = et[i]->get_logn_auto() ;
	      const Metrique& gamma = et[i]->get_met_gamma() ;
	      Map_af mapping (et[i]->get_mp()) ; 

	      Tenseur vecteur (mapping, 1, COV, mapping.get_bvect_spher()) ;  
	      vecteur.set_etat_qcq() ;
	      vecteur = logn_auto.derive_con(gamma) ;
	      vecteur.change_triad(mapping.get_bvect_spher()) ;
	      Cmp integrant (vecteur(0)) ;

	      *p_mass_kom += mapping.integrale_surface_infini (integrant) / qpig ;
		
	    }	// End of the case where a new computation was necessary
    }
    return *p_mass_kom ; 
    
}


		    //---------------------------------//
		    //	 Total angular momentum        //
		    //---------------------------------//

const Tbl& Bin_ns_ncp::angu_mom() const {
    
    if (p_angu_mom == 0x0) {	    // a new computation is requireed
	
	p_angu_mom = new Tbl(3) ; 
	
	p_angu_mom->annule_hard() ;	// fills the double array with zeros
	    
	for(int i=0; i<=1; i++) {

	  Map_af mapping (et[i]->get_mp()) ; 
	  const Tenseur& kij_auto = et[i]->get_tkij_auto() ;
	  
	  Cmp xx(mapping) ;
	  Cmp yy(mapping) ;
	  Cmp zz(mapping) ;

	  xx = mapping.xa ;
	  yy = mapping.ya ;
	  zz = mapping.za ;
	  
	  // X component
	  // -----------

	  Tenseur vecteur_x(mapping, 1, CON, mapping.get_bvect_cart()) ;
	  vecteur_x.set_etat_qcq() ;
	  vecteur_x = yy*kij_auto(2) - zz*kij_auto(1) ;
	  vecteur_x.change_triad(mapping.get_bvect_spher()) ;
	  Cmp integrant_x (vecteur_x(0)) ;

	  p_angu_mom->set(0) += mapping.integrale_surface_infini (integrant_x) / (8*M_PI) ;

 	  // Y component
	  // -----------

	  Tenseur vecteur_y(mapping, 1, CON, mapping.get_bvect_cart()) ;
	  vecteur_y.set_etat_qcq() ;
	  vecteur_y = zz*kij_auto(0) - xx*kij_auto(2) ;
	  vecteur_y.change_triad(mapping.get_bvect_spher()) ;
	  Cmp integrant_y (vecteur_y(0)) ;

	  p_angu_mom->set(0) += mapping.integrale_surface_infini (integrant_y) / (8*M_PI) ;
 
	  // Z component
	  // -----------

	  Tenseur vecteur_z(mapping, 1, CON, mapping.get_bvect_cart()) ;
	  vecteur_z.set_etat_qcq() ;
	  vecteur_z = xx*kij_auto(1) - yy*kij_auto(0) ;
	  vecteur_z.change_triad(mapping.get_bvect_spher()) ;
	  Cmp integrant_z (vecteur_z(0)) ;

	  p_angu_mom->set(0) += mapping.integrale_surface_infini (integrant_z) / (8*M_PI) ;

	}

    }	// End of the case where a new computation was necessary
    
    return *p_angu_mom ; 
    
}



		    //---------------------------------//
		    //		Total energy	       //
		    //---------------------------------//

double Bin_ns_ncp::total_ener() const {
    
    if (p_total_ener == 0x0) {	    // a new computation is requireed
	
	p_total_ener = new double ; 
	    
	    *p_total_ener = mass_adm() - star1.mass_b() - star2.mass_b() ; 
	    
    }	// End of the case where a new computation was necessary
    
    
    return *p_total_ener ; 
    
}


		    //---------------------------------//
		    //	 Error on the virial theorem   //
		    //---------------------------------//

double Bin_ns_ncp::virial() const {
    
    if (p_virial == 0x0) {	    // a new computation is requireed
	
	p_virial = new double ; 
	    
	    *p_virial = 1. - mass_kom() / mass_adm() ; 
	    
	}

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


	     //------------------------------------------------//
	     //	 Virial error by Friedman, Uryu, and Shibata   //
	     //------------------------------------------------//

double Bin_ns_ncp::virial_fus() const {
    
    if (p_virial_fus == 0x0) {	    // a new computation is requireed
	
	p_virial_fus = new double ; 
	    
	if (star1.is_relativistic()) {	// Relativistic case
					// -----------------

#include "unites.h"
    if (this == 0x0) {	// To avoid any compilation warning
	cout << f_unit << msol << km << mevpfm3 ;
    }
	
	    assert( star2.is_relativistic() ) ; 
	    
	    *p_virial_fus = 0 ;

	    double vir_pres = 0. ;
	    double vir_extr = 0. ;
	    double vir_grav = 0. ;

	    for (int i=0; i<=1; i++) {  // loop on the stars

	      const Cmp& lapse = (et[i]->get_nnn())() ;
	      const Cmp& a2 = (et[i]->get_a_car())() ;
	      const Cmp& se = (et[i]->get_s_euler())() ;
	      const Cmp& ak2_auto = (et[i]->get_akcar_auto())() ;
	      const Cmp& ak2_comp = (et[i]->get_akcar_comp())() ;

	      const Tenseur& dnu_auto = et[i]->get_d_logn_auto() ;
	      const Tenseur& dnu_comp = et[i]->get_d_logn_comp() ;
	      const Tenseur& dbe_auto = et[i]->get_d_beta_auto() ;
	      const Tenseur& dbe_comp = et[i]->get_d_beta_comp() ;

	      Cmp source = 2. * lapse * a2 * sqrt(a2) * se ;
	      vir_pres += source.integrale() ;

	      source = 1.5 * lapse * sqrt(a2) * (ak2_auto + ak2_comp) / qpig ;
	      vir_extr += source.integrale() ;

	      Tenseur sprod = flat_scalar_prod( dbe_auto, dbe_auto+dbe_comp )
		- flat_scalar_prod( dnu_auto, dnu_auto+dnu_comp ) ;

	      source = lapse * sqrt(a2) * sprod() / qpig ;
	      vir_grav += source.integrale() ;

	    }  // End of the loop on the stars


	    *p_virial_fus = (vir_pres + vir_extr + vir_grav) / mass_adm() ;
	    
	}
	else {		// Newtonian case 
			// --------------
			
	    *p_virial_fus = virial() ; 
	    
		
	}   // End of the Newtonian case 

    }	// End of the case where a new computation was necessary
    
    return *p_virial_fus ; 
    
}

*/
