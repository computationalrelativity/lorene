/*
 *  Basic methods for class Bin_ns_bh
 *
 */

/*
 *   Copyright (c) 2002  Philippe Grandclement, Keisuke Taniguchi,
 *              Eric Gourgoulhon
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

char bin_ns_bh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2004/03/25 10:28:58  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2003/02/13 16:40:25  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 * Revision 1.3  2002/12/19 14:51:19  e_gourgoulhon
 * Added the new functions set_omega and set_x_axe
 *
 * Revision 1.2  2002/12/18 10:31:15  e_gourgoulhon
 * irrot : int -> bool
 *
 * Revision 1.1  2002/12/17 13:10:11  e_gourgoulhon
 * Methods for class Bin_ns_bh
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <math.h>

// Lorene headers
#include "bin_ns_bh.h"
#include "utilitaires.h"
#include "unites.h"

  			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
Bin_ns_bh::Bin_ns_bh(Map& mp_ns, int nzet, const Eos& eos, bool irrot_ns,
	        Map_af& mp_bh)
		: ref_triad(0., "Absolute frame Cartesian basis"),
		  star(mp_ns, nzet, true, eos, irrot_ns, ref_triad),
		  hole(mp_bh),
		  omega(0),
		  x_axe(0) {

     // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;
}

// Copy constructor
// ----------------
Bin_ns_bh::Bin_ns_bh(const Bin_ns_bh& bibi)
                : ref_triad(0., "Absolute frame Cartesian basis"),
		 star(bibi.star),
		 hole(bibi.hole),
		 omega(bibi.omega),
		 x_axe(bibi.x_axe) {

     // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;
}

// Constructor from a file
// -----------------------
Bin_ns_bh::Bin_ns_bh(Map& mp_ns, const Eos& eos, Map_af& mp_bh, FILE* fich)
		: ref_triad(0., "Absolute frame Cartesian basis"),
		  star(mp_ns, eos, ref_triad, fich),
		  hole(mp_bh, fich) {

    // omega and x_axe are read in the file:
    fread_be(&omega, sizeof(double), 1, fich) ;
    fread_be(&x_axe, sizeof(double), 1, fich) ;

    assert(hole.get_omega() == omega) ;

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;
}

			    //------------//
			    // Destructor //
			    //------------//

Bin_ns_bh::~Bin_ns_bh(){

    del_deriv() ;

}

			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Bin_ns_bh::del_deriv() const {

    if (p_mass_adm != 0x0) delete p_mass_adm ;
    if (p_mass_kom != 0x0) delete p_mass_kom ;
    if (p_angu_mom != 0x0) delete p_angu_mom ;
    if (p_total_ener != 0x0) delete p_total_ener ;
    if (p_virial != 0x0) delete p_virial ;
    if (p_virial_gb != 0x0) delete p_virial_gb ;
    if (p_virial_fus != 0x0) delete p_virial_fus ;
    if (p_ham_constr != 0x0) delete p_ham_constr ;
    if (p_mom_constr != 0x0) delete p_mom_constr ;

    set_der_0x0() ;
}




void Bin_ns_bh::set_der_0x0() const {

    p_mass_adm = 0x0 ;
    p_mass_kom = 0x0 ;
    p_angu_mom = 0x0 ;
    p_total_ener = 0x0 ;
    p_virial = 0x0 ;
    p_virial_gb = 0x0 ;
    p_virial_fus = 0x0 ;
    p_ham_constr = 0x0 ;
    p_mom_constr = 0x0 ;

}

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Binaire
// -----------------------------

void Bin_ns_bh::operator=(const Bin_ns_bh& bibi) {

    assert( bibi.ref_triad == ref_triad ) ;

    star = bibi.star ;
    hole = bibi.hole ;

    omega = bibi.omega ;
    x_axe = bibi.x_axe ;

    // ref_triad remains unchanged

    del_deriv() ;  // Deletes all derived quantities

}

void Bin_ns_bh::set_omega(double omega_i) {

        omega = omega_i ;
	hole.set_omega(omega) ;

	del_deriv() ;

}

void Bin_ns_bh::set_x_axe(double x_axe_i) {

        x_axe = x_axe_i ;

	del_deriv() ;

}



			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Bin_ns_bh::sauve(FILE* fich) const {

    star.sauve(fich) ;
    hole.sauve(fich) ;

    fwrite_be(&omega, sizeof(double), 1, fich) ;
    fwrite_be(&x_axe, sizeof(double), 1, fich) ;

}

// Printing
// --------
ostream& operator<<(ostream& ost, const Bin_ns_bh& bibi)  {
    bibi >> ost ;
    return ost ;
}


ostream& Bin_ns_bh::operator>>(ostream& ost) const {

  using namespace Unites ;

    ost << endl ;
    ost << "Neutron star - black hole binary system" << endl ;
    ost << "=======================================" << endl ;
    ost << endl <<
	"Orbital angular velocity : " << omega * f_unit << " rad/s" << endl ;
    ost <<
	"Absolute coordinate X of the rotation axis : " << x_axe / km
	    << " km" << endl ;
    ost << endl << "Neutron star : " << endl ;
    ost <<         "============   " << endl ;
    ost << star << endl ;


    ost << "Black hole : " << endl ;
    ost << "==========   " << endl ;
    ost << "Coordinate radius of the throat : " << hole.get_rayon() / km << " km" << endl ;
    ost << "Absolute abscidia of the throat center : " << (hole.get_mp()).get_ori_x() / km
        << " km" << endl ;
    return ost ;
}
