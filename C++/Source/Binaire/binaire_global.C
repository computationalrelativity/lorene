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
#include "unites.h"

		    //---------------------------------//
		    //		ADM mass	       //
		    //---------------------------------//

double Binaire::mass_adm() const {
    
    if (p_mass_adm == 0x0) {	    // a new computation is requireed
	
	p_mass_adm = new double ; 
	    
	if (star1.is_relativistic()) {	// Relativistic case
					// -----------------
	
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
	    
	    *p_virial = - 1000 ; 
	    
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


