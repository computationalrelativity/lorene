/*
 * Methods of class Binaire to compute global quantities
 *
 * (see file binaire.h for documentation)
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
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


char binaire_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2001/12/20 14:18:40  k_taniguchi
 * Addition of the Komar mass, the virial error by Gourgoulhon and Bonazzola, and the virial error by Friedman, Uryu, and Shibata.
 *
 * Revision 1.3  2001/12/16 16:21:38  e_gourgoulhon
 * #include "unites.h" is now local to Binaire::mass_adm(), in order
 * not to make Lorene's units global variables.
 *
 * Revision 1.2  2001/12/14 09:45:14  k_taniguchi
 * Correction of missing 16 pi G factor in the ADM mass
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2000/03/08  12:26:33  eric
 * Ajout de l'appel a std_base_scal() sur le Cmp source dans le cas
 * relativiste (masse ADM).
 *
 * Revision 2.2  2000/02/23  11:26:00  keisuke
 * Changement de "virial relation".
 *
 * Revision 2.1  2000/02/18  15:48:55  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/02/18  14:53:09  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "binaire.h"

		    //---------------------------------//
		    //		ADM mass	       //
		    //---------------------------------//

double Binaire::mass_adm() const {
    
    if (p_mass_adm == 0x0) {	    // a new computation is requireed
	
	p_mass_adm = new double ; 
	    
	if (star1.is_relativistic()) {	// Relativistic case
					// -----------------
#include "unites.h"
    if (this == 0x0) {	// To avoid any compilation warning
	cout << f_unit << msol << km << mevpfm3 ;
    }
	
	    assert( star2.is_relativistic() ) ;
	    
	    *p_mass_adm = 0 ; 
	    
	    for (int i=0; i<=1; i++) {	    // loop on the stars
	    
		const Cmp& a2 = (et[i]->get_a_car())() ; 
		const Cmp& ee = (et[i]->get_ener_euler())() ; 
		const Cmp& ak2_auto = (et[i]->get_akcar_auto())() ;
		const Cmp& ak2_comp = (et[i]->get_akcar_comp())() ;
	    
		Cmp source = pow(a2, 1.25) * ee 
		  + pow(a2, 0.25) * (ak2_auto + ak2_comp) / (4.*qpig) ; 
			   
		source.std_base_scal() ; 
			   
		*p_mass_adm += source.integrale() ; 
	    
	    }    
	
	}
	else {		// Newtonian case 
			// --------------
			
	    *p_mass_adm = star1.mass_b() + star2.mass_b() ; 
	    
	}
		
    }	// End of the case where a new computation was necessary
    
    return *p_mass_adm ; 
    
}


		    //---------------------------------//
		    //		Komar mass	       //
		    //---------------------------------//

double Binaire::mass_kom() const {
    
    if (p_mass_kom == 0x0) {	    // a new computation is requireed
	
	p_mass_kom = new double ; 
	    
	if (star1.is_relativistic()) {	// Relativistic case
					// -----------------

#include "unites.h"
    if (this == 0x0) {	// To avoid any compilation warning
	cout << f_unit << msol << km << mevpfm3 ;
    }
	
	    assert( star2.is_relativistic() ) ; 
	    
	    *p_mass_kom = 0 ; 
	    
	    for (int i=0; i<=1; i++) {	    // loop on the stars
	    
	        const Cmp& lapse = (et[i]->get_nnn())() ;
	        const Cmp& a2 = (et[i]->get_a_car())() ; 
		const Cmp& ee = (et[i]->get_ener_euler())() ; 
		const Cmp& se = (et[i]->get_s_euler())() ; 
		const Cmp& ak2_auto = (et[i]->get_akcar_auto())() ;
		const Cmp& ak2_comp = (et[i]->get_akcar_comp())() ;

		const Tenseur& dnu_auto = et[i]->get_d_logn_auto() ;
		const Tenseur& dnu_comp = et[i]->get_d_logn_comp() ;
		const Tenseur& dbe_auto = et[i]->get_d_beta_auto() ;
		const Tenseur& dbe_comp = et[i]->get_d_beta_comp() ;

		Tenseur dndb = flat_scalar_prod(dnu_auto,
						dbe_auto + dbe_comp) ;
		Tenseur dndn = flat_scalar_prod(dnu_auto,
						dnu_auto + dnu_comp) ;

		Cmp source = lapse * ( a2 * (ee + se)
				       + (ak2_auto + ak2_comp)/qpig
				       - dndb()/qpig + dndn()/qpig ) ;
			   
		source.std_base_scal() ; 
			   
		*p_mass_kom += source.integrale() ; 
	    
	    }    
	
	}
	else {		// Newtonian case 
			// --------------
			
	    *p_mass_kom = star1.mass_b() + star2.mass_b() ; 
	    
	}
		
    }	// End of the case where a new computation was necessary
    
    return *p_mass_kom ; 
    
}


		    //---------------------------------//
		    //	 Total angular momentum        //
		    //---------------------------------//

const Tbl& Binaire::angu_mom() const {
    
    if (p_angu_mom == 0x0) {	    // a new computation is requireed
	
	p_angu_mom = new Tbl(3) ; 
	
	p_angu_mom->annule_hard() ;	// fills the double array with zeros
	    
	for (int i=0; i<=1; i++) {	    // loop on the stars
	    
	    const Map& mp = et[i]->get_mp() ; 
	    
	    Cmp xx(mp) ; 
	    Cmp yy(mp) ; 
	    Cmp zz(mp) ; 
	    
	    xx = mp.xa ;
	    yy = mp.ya ;
	    zz = mp.za ;
	    
	    const Cmp& vx = (et[i]->get_u_euler())(0) ; 
	    const Cmp& vy = (et[i]->get_u_euler())(1) ; 
	    const Cmp& vz = (et[i]->get_u_euler())(2) ; 

	    Cmp rho(mp) ; 
	    
	    if ( et[i]->is_relativistic() ) {
		const Cmp& a2 = (et[i]->get_a_car())() ; 
		const Cmp& ee = (et[i]->get_ener_euler())() ; 
		const Cmp& pp = (et[i]->get_press())() ; 
		rho = pow(a2, 2.5) * (ee + pp) ; 
	    }
	    else {
		rho = (et[i]->get_nbar())() ;
	    }
	    

	    Base_val** base = (et[i]->get_mp()).get_mg()->std_base_vect_cart() ;

	    // X component
	    // -----------
	    
	    Cmp source = rho * ( yy * vz  -  zz * vy ) ;
	     
	    (source.va).set_base( *(base[2]) ) ;    // same basis as V^z
	    
//##	    p_angu_mom->set(0) += source.integrale() ; 

	    p_angu_mom->set(0) += 0 ; 
	    
	    // y component
	    // -----------
	    
	    source = rho * ( zz * vx  -  xx * vz ) ;
	    
	    (source.va).set_base( *(base[2]) ) ;    // same basis as V^z
	    
//##	    p_angu_mom->set(1) += source.integrale() ; 
	    p_angu_mom->set(1) += 0 ; 
	
	    
	    // Z component
	    // -----------
	    
	    source = rho * ( xx * vy - yy * vx ) ;
	    
	    source.std_base_scal() ;	// same basis as V^x (standard scalar
					//    field)
	     
	    p_angu_mom->set(2) += source.integrale() ; 
	    
	    delete base[0] ; 
	    delete base[1] ; 
	    delete base[2] ; 
	    delete [] base ; 
	    
	}  // End of the loop on the stars

    }	// End of the case where a new computation was necessary
    
    return *p_angu_mom ; 
    
}




		    //---------------------------------//
		    //		Total energy	       //
		    //---------------------------------//

double Binaire::total_ener() const {
    
    if (p_total_ener == 0x0) {	    // a new computation is requireed
	
	p_total_ener = new double ; 
	    
	if (star1.is_relativistic()) {	// Relativistic case
					// -----------------
	
	    assert( star2.is_relativistic() ) ; 
	    
	    *p_total_ener = mass_adm() - star1.mass_b() - star2.mass_b() ; 
	    
	}
	else {		// Newtonian case 
			// --------------
			
	    *p_total_ener = 0 ; 
	    
	    for (int i=0; i<=1; i++) {	    // loop on the stars
	    
		const Cmp e_int = (et[i]->get_ener())()
				    - (et[i]->get_nbar())()  ; 

		const Cmp& rho = (et[i]->get_nbar())() ;

		// Fluid velocity with respect to the inertial frame
		const Tenseur& vit = et[i]->get_u_euler() ; 
		
		Cmp vit2 = flat_scalar_prod(vit, vit)() ; 
		
		// Gravitational potential 
		const Cmp nu = (et[i]->get_logn_auto())() 
			       + (et[i]->get_logn_comp())() ;
	    
		Cmp source = e_int + .5 * rho * vit2 + .5 * rho * nu ; 
			   
		*p_total_ener += source.integrale() ; 
	    
	    
	    }   // End of the loop on the stars
	
	}   // End of Newtonian case	
    
    }	// End of the case where a new computation was necessary
    
    
    return *p_total_ener ; 
    
}


		    //---------------------------------//
		    //	 Error on the virial theorem   //
		    //---------------------------------//

double Binaire::virial() const {
    
    if (p_virial == 0x0) {	    // a new computation is requireed
	
	p_virial = new double ; 
	    
	if (star1.is_relativistic()) {	// Relativistic case
					// -----------------
	
	    assert( star2.is_relativistic() ) ; 
	    
	    *p_virial = 1. - mass_kom() / mass_adm() ; 
	    
	}
	else {		// Newtonian case 
			// --------------
			
	    *p_virial = 0 ; 
	    
	    
	    double vir_mat = 0 ; 
	    double vir_grav = 0 ; 
	    
	    for (int i=0; i<=1; i++) {	    // loop on the stars
	    
		const Cmp& pp = (et[i]->get_press())()  ; 

		const Cmp& rho = (et[i]->get_nbar())() ;

		// Fluid velocity with respect to the inertial frame
		const Tenseur& vit = et[i]->get_u_euler() ; 
		
		Cmp vit2 = flat_scalar_prod(vit, vit)() ; 
		
		// Gravitational potential 
		const Cmp nu = (et[i]->get_logn_auto())() 
			       + (et[i]->get_logn_comp())() ;
	    
		Cmp source = 3*pp + rho * vit2 ;
		
		vir_mat +=  source.integrale() ;
		 
		source =  .5 * rho * nu ; 

		vir_grav +=  source.integrale() ;
	    
	    }  // End of the loop on the stars
	
	    *p_virial = ( vir_mat + vir_grav ) / fabs(vir_grav) ;
		
	}   // End of the Newtonian case 

    }	// End of the case where a new computation was necessary
    
    return *p_virial ; 
    
}


	     //----------------------------------------------//
	     //	 Virial error by Gourgoulhon and Bonazzola   //
	     //----------------------------------------------//

double Binaire::virial_gb() const {
    
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

double Binaire::virial_fus() const {
    
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

