/*
 * Methods for computing global quantities within the class Etoile_bin
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.13  2016/12/05 16:17:53  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.12  2014/10/13 08:52:55  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.11  2004/12/30 17:24:42  f_limousin
 * Remove #include "et_bin_ncp.h"
 *
 * Revision 1.10  2004/12/30 17:22:25  f_limousin
 * Remove all functions belonging to the class Et_bin_ncp
 *
 * Revision 1.9  2003/06/20 14:02:44  f_limousin
 * pow(gamma(), 1./2.) is changed to  pow(a_car(), 3./2.)
 *
 * Revision 1.8  2003/02/06 17:25:56  f_limousin
 * Add global quantities for the class et_bin_ncp
 *
 * Revision 1.7  2003/01/17 13:32:27  f_limousin
 * Add comments
 *
 * Revision 1.6  2002/12/17 09:48:21  k_taniguchi
 * Suppress every desaliasing multiplications using "%".
 *
 * Revision 1.5  2002/12/16 17:32:47  k_taniguchi
 * Suppress the things I did in the previous version.
 *
 * Revision 1.4  2002/12/16 16:59:39  k_taniguchi
 * Set some Cmp to the state of "std_base_scal()".
 *
 * Revision 1.3  2002/12/16 14:36:39  k_taniguchi
 * Introduce a new Cmp for the calculation of gravitational mass.
 *
 * Revision 1.2  2002/12/10 15:45:25  k_taniguchi
 * Change the multiplication "*" to "%".
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/07/06  10:02:22  eric
 * *** empty log message ***
 *
 * Revision 2.3  2000/07/06  09:40:37  eric
 * Ajout de la fonction xa_barycenter().
 *
 * Revision 2.2  2000/02/02  09:23:12  eric
 * 1ere version operationnelle dans le cas relativiste.
 *
 * Revision 2.1  2000/02/01  16:00:13  eric
 * Le calcul de mass_b est implemente dans le cas relativiste.
 *
 * Revision 2.0  2000/01/31  15:57:21  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C

// Headers Lorene
#include "etoile.h"

			//--------------------------//
			//	Baryon mass	    //
			//--------------------------//

namespace Lorene {
double Etoile_bin::mass_b() const {

    if (p_mass_b == 0x0) {    // a new computation is required
	
	if (relativistic) {

	    Cmp sqrt_acar = sqrt(a_car()) ;
	    sqrt_acar.std_base_scal() ;
	  
	    // See Eq (69) from Gourgoulhon et al. (2001)
	    // Cmp dens = a_car() % sqrt_acar % gam_euler() % nbar() ;

	    Cmp dens = a_car() * sqrt( a_car() ) * gam_euler() * nbar() ;
	    
	    dens.std_base_scal() ; 

	    p_mass_b = new double( dens.integrale() ) ;

	}
	else{
	    assert(nbar.get_etat() == ETATQCQ) ; 

	    p_mass_b = new double( nbar().integrale() ) ;

	}

    }
    
    return *p_mass_b ; 

} 


			//----------------------------//
			//	Gravitational mass    //
			//----------------------------//

double Etoile_bin::mass_g() const {

    if (p_mass_g == 0x0) {    // a new computation is required
	
	if (relativistic) {


	    Cmp sqrt_acar = sqrt(a_car()) ;
	    sqrt_acar.std_base_scal() ;
	    Cmp dens = a_car() * sqrt( a_car() ) * nnn()
		* ( ener_euler() + s_euler() ) ;

	    dens.std_base_scal() ; 

	    p_mass_g = new double( dens.integrale() ) ;

	}
	else{
	    p_mass_g = new double( mass_b() ) ;   // in the Newtonian case
						    //  M_g = M_b
	}
    }
    
    return *p_mass_g ; 

} 


		
			//----------------------------------//
			//  X coordinate of the barycenter  //
			//----------------------------------//

double Etoile_bin::xa_barycenter() const {

    if (p_xa_barycenter == 0x0) {    // a new computation is required
	
	Cmp xxa(mp) ; 
	xxa = mp.xa ;	// Absolute X coordinate
	xxa.std_base_scal() ;

	Cmp dens = a_car() * sqrt( a_car() ) * gam_euler() * nbar() * xxa ; 
	
	dens.std_base_scal() ; 

	p_xa_barycenter = new double( dens.integrale() / mass_b() ) ;
	
    }
    
    return *p_xa_barycenter ; 

}



}
