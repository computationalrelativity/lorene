/*
 * Methods for computing global quantities within the class Etoile_rot
 *
 * (see file etoile.h for documentation)
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char et_rot_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.5  2000/11/19  18:52:09  eric
 * grv2() operationnelle.
 *
 * Revision 1.4  2000/10/12  15:34:55  eric
 * Calcul de la masse grav, de GRV3 et du moment quadrupolaire.
 *
 * Revision 1.3  2000/08/31  11:25:58  eric
 * *** empty log message ***
 *
 * Revision 1.2  2000/08/25  12:28:16  eric
 * *** empty log message ***
 *
 * Revision 1.1  2000/07/20  15:32:56  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <math.h>

// Headers Lorene
#include "etoile.h"

			//--------------------------//
			//	Baryon mass	    //
			//--------------------------//

double Etoile_rot::mass_b() const {

    if (p_mass_b == 0x0) {    // a new computation is required
	
	if (relativistic) {

	    Cmp dens = a_car() * bbb() * gam_euler() * nbar() ;
	    
	    dens.std_base_scal() ; 

	    p_mass_b = new double( dens.integrale() ) ;


	}
	else{  // Newtonian case 
	    assert(nbar.get_etat() == ETATQCQ) ; 

	    p_mass_b = new double( nbar().integrale() ) ;

	}

    }
    
    return *p_mass_b ; 

} 


			//----------------------------//
			//	Gravitational mass    //
			//----------------------------//

double Etoile_rot::mass_g() const {

    if (p_mass_g == 0x0) {    // a new computation is required
	
	if (relativistic) {

	    Tenseur source = nnn * (ener_euler + s_euler) 
				+ 2 * bbb * (ener_euler + press)
				    * tnphi * uuu ; 
	    source = a_car * bbb * source ;

	    source.set_std_base() ;

	    p_mass_g = new double( source().integrale() ) ;


	}
	else{  // Newtonian case 
	    p_mass_g = new double( mass_b() ) ;   // in the Newtonian case
						    //  M_g = M_b
	}
    }
    
    return *p_mass_g ; 

} 
		
			//----------------------------//
			//	Angular momentum      //
			//----------------------------//

double Etoile_rot::angu_mom() const {

    if (p_angu_mom == 0x0) {    // a new computation is required
	
	Cmp dens = uuu() ; 

	dens.mult_r() ;			//  Multiplication by
	dens.va = (dens.va).mult_st() ;	//    r sin(theta)

	if (relativistic) {
	    dens = a_car() * b_car() * (ener_euler() + press()) 
			* dens ; 
	}
	else {    // Newtonian case 
	    dens = nbar() * dens ; 
	}

	dens.std_base_scal() ; 

	p_angu_mom = new double( dens.integrale() ) ;

    }
    
    return *p_angu_mom ; 

}


			//----------------------------//
			//	     T/W	      //
			//----------------------------//

double Etoile_rot::tsw() const {

    if (p_tsw == 0x0) {    // a new computation is required
	
	double tcin = 0.5 * omega * angu_mom() ;
	
	if (relativistic) {
	    
	    Cmp dens = a_car() * bbb() * gam_euler() * ener() ;
	    dens.std_base_scal() ; 
	    double mass_p = dens.integrale() ; 
	    
	    p_tsw = new double( tcin / ( mass_p + tcin - mass_g() ) ) ;  	
	   
	}
	else {	    // Newtonian case 
	    Cmp dens = 0.5 * nbar() * logn() ;
	    dens.std_base_scal() ; 
	    double wgrav = dens.integrale() ; 
	    p_tsw = new double( tcin / fabs(wgrav) ) ;  
	}


    }
    
    return *p_tsw ; 

}


			//----------------------------//
			//	     GRV2	      //
			//----------------------------//

double Etoile_rot::grv2() const {

    if (p_grv2 == 0x0) {    // a new computation is required
	
		// To get qpig:	
		#include "unites.h"	
		// To avoid some compilation warnings
		if (p_grv2 != 0x0) {
	    	cout << f_unit << msol << km << mevpfm3 << endl ;
		}

        Tenseur sou_m =  2 * qpig * a_car * (press + (ener_euler+press)
        						* uuu*uuu ) ;
        						
        Tenseur sou_q =  1.5 * ak_car
        				 - flat_scalar_prod(logn.gradient_spher(),
						     				logn.gradient_spher() ) ;	

		p_grv2 = new double( double(1) - lambda_grv2(sou_m(), sou_q()) ) ; 	

    }
    
    return *p_grv2 ; 

}


			//----------------------------//
			//	     GRV3	      //
			//----------------------------//

double Etoile_rot::grv3(ostream* ost) const {

    if (p_grv3 == 0x0) {    // a new computation is required

	// To get qpig:	
	#include "unites.h"	    
	// To avoid some compilation warnings
	if (p_grv3 != 0x0) {
	    cout << f_unit << msol << km << mevpfm3 << endl ; 
	}    


	Tenseur source(mp) ; 
	
	// Gravitational term [cf. Eq. (43) of Gourgoulhon & Bonazzola
	// ------------------	    Class. Quantum Grav. 11, 443 (1994)]

	if (relativistic) {
	    Tenseur alpha = dzeta - logn ; 
	    Tenseur beta = log( bbb ) ; 
	    beta.set_std_base() ; 
	    
	    source = 0.75 * ak_car 
		     - flat_scalar_prod(logn.gradient_spher(),
					logn.gradient_spher() )
		     + 0.5 * flat_scalar_prod(alpha.gradient_spher(),
					      beta.gradient_spher() ) ; 
	    
	    Cmp aa = alpha() - 0.5 * beta() ; 
	    Cmp daadt = aa.srdsdt() ;	// 1/r d/dth
	    
	    // What follows is valid only for a mapping of class Map_radial : 
	    const Map_radial* mpr = dynamic_cast<const Map_radial*>(&mp) ; 
	    if (mpr == 0x0) {
		cout << "Etoile_rot::grv3: the mapping does not belong"
		     << " to the class Map_radial !" << endl ; 
		abort() ; 
	    }
		
	    // Computation of 1/tan(theta) * 1/r daa/dtheta
	    if (daadt.get_etat() == ETATQCQ) {
		Valeur& vdaadt = daadt.va ; 
		vdaadt = vdaadt.ssint() ;	// division by sin(theta)
		vdaadt = vdaadt.mult_ct() ;	// multiplication by cos(theta)
	    }
	    
	    Cmp temp = aa.dsdr() + daadt ; 
	    temp = ( bbb() - a_car()/bbb() ) * temp ; 
	    temp.std_base_scal() ; 
	    
	    // Division by r 
	    Valeur& vtemp = temp.va ; 
	    vtemp = vtemp.sx() ;    // division by xi in the nucleus
				    // Id in the shells
				    // division by xi-1 in the ZEC
	    vtemp = (mpr->xsr) * vtemp ; // multiplication by xi/r in the nucleus
					 //		  by 1/r in the shells
					 //		  by r(xi-1) in the ZEC

	    // In the ZEC, a multiplication by r has been performed instead
	    //   of the division: 			
	    temp.set_dzpuis( temp.get_dzpuis() + 2 ) ;  
	    
	    source = bbb() * source() + 0.5 * temp ; 

	}
	else{
	    source = - 0.5 * flat_scalar_prod(logn.gradient_spher(),
					      logn.gradient_spher() ) ; 
	}
	
	source.set_std_base() ; 

	double int_grav = source().integrale() ; 

	// Matter term
	// -----------
	
	if (relativistic) {    
	    source  = qpig * a_car * bbb * s_euler ;
	}
	else{
	    source = qpig * ( 3 * press + nbar * uuu * uuu ) ; 
	}

	source.set_std_base() ; 

	double int_mat = source().integrale() ; 

	// Virial error
	// ------------
	if (ost != 0x0) {
	    *ost << "Etoile_rot::grv3 : gravitational term : " << int_grav 
		 << endl ;
	    *ost << "Etoile_rot::grv3 : matter term :        " << int_mat 
		 << endl ;
	}

	p_grv3 = new double( (int_grav + int_mat) / int_mat ) ; 	 

    }
    
    return *p_grv3 ; 

}


			//----------------------------//
			//	     R_circ	      //
			//----------------------------//

double Etoile_rot::r_circ() const {

    if (p_r_circ == 0x0) {    // a new computation is required
	
	// Index of the point at phi=0, theta=pi/2 at the surface of the star:
	const Mg3d* mg = mp.get_mg() ; 
	assert(mg->get_type_t() == SYM) ; 
	int l_b = nzet - 1 ; 
	int i_b = mg->get_nr(l_b) - 1 ; 
	int j_b = mg->get_nt(l_b) - 1 ; 
	int k_b = 0 ; 
    
	p_r_circ = new double( bbb()(l_b, k_b, j_b, i_b) * ray_eq() ) ; 

    }
    
    return *p_r_circ ; 

}


			//----------------------------//
			//	   Flattening	      //
			//----------------------------//

double Etoile_rot::aplat() const {

    if (p_aplat == 0x0) {    // a new computation is required
	
	p_aplat = new double( ray_pole() / ray_eq() ) ; 	 

    }
    
    return *p_aplat ; 

}


			//----------------------------//
			//	     Z_eq_f	      //
			//----------------------------//

double Etoile_rot::z_eqf() const {

    if (p_z_eqf == 0x0) {    // a new computation is required
	
	// Index of the point at phi=0, theta=pi/2 at the surface of the star:
	const Mg3d* mg = mp.get_mg() ; 
	assert(mg->get_type_t() == SYM) ; 
	int l_b = nzet - 1 ; 
	int i_b = mg->get_nr(l_b) - 1 ; 
	int j_b = mg->get_nt(l_b) - 1 ; 
	int k_b = 0 ; 
    
	double u_eq = uuu()(l_b, k_b, j_b, i_b) ; 
	double n_eq = nnn()(l_b, k_b, j_b, i_b) ; 
	double nphi_eq = nphi()(l_b, k_b, j_b, i_b) ; 
	
	p_z_eqf = new double( sqrt((1.-u_eq)/(1.+u_eq)) 
				/ (n_eq + nphi_eq * r_circ() )
			      - 1. ) ;
    }
    
    return *p_z_eqf ; 

}
			//----------------------------//
			//	     Z_eq_b	      //
			//----------------------------//

double Etoile_rot::z_eqb() const {

    if (p_z_eqb == 0x0) {    // a new computation is required
	
	// Index of the point at phi=0, theta=pi/2 at the surface of the star:
	const Mg3d* mg = mp.get_mg() ; 
	assert(mg->get_type_t() == SYM) ; 
	int l_b = nzet - 1 ; 
	int i_b = mg->get_nr(l_b) - 1 ; 
	int j_b = mg->get_nt(l_b) - 1 ; 
	int k_b = 0 ; 
    
	double u_eq = uuu()(l_b, k_b, j_b, i_b) ; 
	double n_eq = nnn()(l_b, k_b, j_b, i_b) ; 
	double nphi_eq = nphi()(l_b, k_b, j_b, i_b) ; 
	
	p_z_eqb = new double(  sqrt((1.+u_eq)/(1.-u_eq)) 
				/ (n_eq - nphi_eq * r_circ() )
			      - 1. )  ;

    }
    
    return *p_z_eqb ; 

}


			//----------------------------//
			//	     Z_pole	      //
			//----------------------------//

double Etoile_rot::z_pole() const {

    if (p_z_pole == 0x0) {    // a new computation is required
	
	double n_pole = nnn().val_point(ray_pole(), 0., 0.) ; 
	
	p_z_pole = new double(  1. / n_pole - 1. ) ; 

    }
    
    return *p_z_pole ; 

}


			//----------------------------//
			//     Quadrupole moment      //
			//----------------------------//

double Etoile_rot::mom_quad() const {

    if (p_mom_quad == 0x0) {    // a new computation is required
	
	// To get qpig:	
	#include "unites.h"	    
	// To avoid some compilation warnings
	if (p_mom_quad != 0x0) {
	    cout << f_unit << msol << km << mevpfm3 << endl ; 
	}    

	// Source for of the Poisson equation for nu
	// -----------------------------------------

	Tenseur source(mp) ; 
	
	if (relativistic) {
	    Tenseur beta = log(bbb) ; 
	    beta.set_std_base() ; 
	    source =  qpig * a_car *( ener_euler + s_euler ) 
			+ ak_car - flat_scalar_prod(logn.gradient_spher(), 
			    logn.gradient_spher() + beta.gradient_spher()) ; 
	}
	else {
	    source = qpig * nbar ; 
	}
	source.set_std_base() ; 	

	// Multiplication by -r^2 P_2(cos(theta))
	//  [cf Eq.(7) of Salgado et al. Astron. Astrophys. 291, 155 (1994) ]
	// ------------------------------------------------------------------
	
	// Multiplication by r^2 : 
	// ----------------------
	Cmp& csource = source.set() ; 
	csource.mult_r() ; 
	csource.mult_r() ; 
	if (csource.check_dzpuis(2)) {
	    csource.inc2_dzpuis() ; 
	}
		
	// Muliplication by cos^2(theta) :
	// -----------------------------
	Cmp temp = csource ; 
	
	// What follows is valid only for a mapping of class Map_radial : 
	assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ; 
		
	if (temp.get_etat() == ETATQCQ) {
	    Valeur& vtemp = temp.va ; 
	    vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
	    vtemp = vtemp.mult_ct() ;	// multiplication by cos(theta)
	}
	
	// Muliplication by -P_2(cos(theta)) :
	// ----------------------------------
	source = 0.5 * source() - 1.5 * temp ; 
	
	// Final result
	// ------------

	p_mom_quad = new double( source().integrale() / qpig ) ; 	 

    }
    
    return *p_mom_quad ; 

}




