/*
 *  Methods of class Black_hole
 *
 *    (see file blackhole.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2007 Keisuke Taniguchi
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

char blackhole_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:52:45  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:02  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2008/07/03 14:55:26  k_taniguchi
 * Addition of the angular momentum.
 *
 * Revision 1.2  2008/05/15 19:25:07  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:18:47  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <cmath>

// Lorene headers
#include "blackhole.h"
#include "utilitaires.h"
#include "unites.h"

                    //----------------------//
                    //     Constructors     //
                    //----------------------//

// Standard constructor
// --------------------
namespace Lorene {
Black_hole::Black_hole(Map& mp_i, bool kerrschild_i, double massbh)
      : mp(mp_i),
	kerrschild(kerrschild_i),
	mass_bh(massbh),
	lapconf(mp_i),
	lapconf_rs(mp_i),
	lapconf_bh(mp_i),
	lapse(mp_i),
	shift(mp_i, CON, mp_i.get_bvect_cart()),
	shift_rs(mp_i, CON, mp_i.get_bvect_cart()),
	shift_bh(mp_i, CON, mp_i.get_bvect_cart()),
	confo(mp_i),
	taij(mp_i, CON, mp_i.get_bvect_cart()),
	taij_rs(mp_i, CON, mp_i.get_bvect_cart()),
	taij_quad(mp_i),
	taij_quad_rs(mp_i),
	flat(mp_i, mp_i.get_bvect_cart()) {

    // Pointers of derived quantities are initialized to zero
    set_der_0x0() ;

    // The metric quantities are initialized to the flat one
    lapconf = 1. ;
    lapconf.std_spectral_base() ;
    lapconf_rs = 0. ;
    lapconf_rs.std_spectral_base() ;
    lapconf_bh = 1. ;
    lapconf_bh.std_spectral_base() ;
    lapse = 1. ;
    lapse.std_spectral_base() ;
    shift.set_etat_zero() ;
    shift_rs.set_etat_zero() ;
    shift_bh.set_etat_zero() ;
    confo = 1. ;
    confo.std_spectral_base() ;

    taij.set_etat_zero() ;
    taij_rs.set_etat_zero() ;
    taij_quad = 0. ;
    taij_quad.std_spectral_base() ;
    taij_quad_rs = 0. ;
    taij_quad_rs.std_spectral_base() ;

}

// Copy constructor
// ----------------
Black_hole::Black_hole(const Black_hole& bh)
      : mp(bh.mp),
	kerrschild(bh.kerrschild),
	mass_bh(bh.mass_bh),
	lapconf(bh.lapconf),
	lapconf_rs(bh.lapconf_rs),
	lapconf_bh(bh.lapconf_bh),
	lapse(bh.lapse),
	shift(bh.shift),
	shift_rs(bh.shift_rs),
	shift_bh(bh.shift_bh),
	confo(bh.confo),
	taij(bh.taij),
	taij_rs(bh.taij_rs),
	taij_quad(bh.taij_quad),
	taij_quad_rs(bh.taij_quad_rs),
	flat(bh.flat) {

    set_der_0x0() ;

}

// Constructor from a file
// -----------------------
Black_hole::Black_hole(Map& mp_i, FILE* fich)
      : mp(mp_i),
	lapconf(mp_i),
	lapconf_rs(mp_i, *(mp_i.get_mg()), fich),
	lapconf_bh(mp_i),
	lapse(mp_i),
	shift(mp_i, CON, mp_i.get_bvect_cart()),
	shift_rs(mp_i, mp_i.get_bvect_cart(), fich),
	shift_bh(mp_i, CON, mp_i.get_bvect_cart()),
	confo(mp_i, *(mp_i.get_mg()), fich),
	taij(mp_i, CON, mp_i.get_bvect_cart()),
	taij_rs(mp_i, CON, mp_i.get_bvect_cart()),
	taij_quad(mp_i),
	taij_quad_rs(mp_i),
	flat(mp_i, mp_i.get_bvect_cart()) {

    // Background
    // ----------
    fread(&kerrschild, sizeof(bool), 1, fich) ;
    fread(&mass_bh, sizeof(double), 1, fich) ;

    // All other fields are initialized to zero or some values
    // -------------------------------------------------------
    lapconf = lapconf_rs ;
    lapconf.std_spectral_base() ;
    lapconf_bh = 0. ;
    lapconf_bh.std_spectral_base() ;

    lapse = lapconf / confo ;

    shift = shift_rs ;
    shift.std_spectral_base() ;
    shift_bh.set_etat_zero() ;

    taij.set_etat_zero() ;
    taij_rs.set_etat_zero() ;
    taij_quad = 0. ;
    taij_quad.std_spectral_base() ;
    taij_quad_rs = 0. ;
    taij_quad_rs.std_spectral_base() ;

    // Pointers of derived quantities are initialized to zero
    // ------------------------------------------------------
    set_der_0x0() ;

}


                    //--------------------//
                    //     Destructor     //
                    //--------------------//

Black_hole::~Black_hole() {

    del_deriv() ;

}


                    //------------------------------------------//
                    //     Management of derived quantities     //
                    //------------------------------------------//

void Black_hole::del_deriv() const {

    if (p_mass_irr != 0x0) delete p_mass_irr ;
    if (p_mass_adm != 0x0) delete p_mass_adm ;
    if (p_mass_kom != 0x0) delete p_mass_kom ;
    if (p_rad_ah != 0x0) delete p_rad_ah ;
    if (p_spin_am_bh != 0x0) delete p_spin_am_bh ;
    if (p_angu_mom_bh != 0x0) delete p_angu_mom_bh ;

    set_der_0x0() ;

}

void Black_hole::set_der_0x0() const {

    p_mass_irr = 0x0 ;
    p_mass_adm = 0x0 ;
    p_mass_kom = 0x0 ;
    p_rad_ah = 0x0 ;
    p_spin_am_bh = 0x0 ;
    p_angu_mom_bh = 0x0 ;

}


                    //--------------------//
                    //     Assignment     //
                    //--------------------//

// Assignment to another Black_hole
// --------------------------------
void Black_hole::operator=(const Black_hole& bh) {

    assert( &(bh.mp) == &mp ) ;     // Same mapping

    kerrschild = bh.kerrschild ;
    mass_bh = bh.mass_bh ;
    lapconf = bh.lapconf ;
    lapconf_rs = bh.lapconf_rs ;
    lapconf_bh = bh.lapconf_bh ;
    lapse = bh.lapse ;
    shift = bh.shift ;
    shift_rs = bh.shift_rs ;
    shift_bh = bh.shift_bh ;
    confo = bh.confo ;
    taij = bh.taij ;
    taij_rs = bh.taij_rs ;
    taij_quad = bh.taij_quad ;
    taij_quad_rs = bh.taij_quad_rs ;
    flat = bh.flat ;

    del_deriv() ;     // Deletes all derived quantities

}


                    //-----------------//
                    //     Outputs     //
                    //-----------------//

// Save in a file
// --------------
void Black_hole::sauve(FILE* fich) const {

    lapconf_rs.sauve(fich) ;
    shift_rs.sauve(fich) ;
    confo.sauve(fich) ;

    fwrite(&kerrschild, sizeof(bool), 1, fich) ;
    fwrite(&mass_bh, sizeof(double), 1, fich) ;

}

// Printing
// --------
ostream& operator<<(ostream& ost, const Black_hole& bh) {

    bh >> ost ;
    return ost ;

}

ostream& Black_hole::operator>>(ostream& ost) const {

  using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    int nt = mg->get_nt(1) ;

    ost << endl ;
    if (kerrschild) {
        ost << "Kerr-Schild background" << endl ;
	ost << "----------------------" << endl ;
    }
    else {
        ost << "Conformally flat background" << endl ;
	ost << "---------------------------" << endl ;
    }

    ost << "lapconf on the AH :  "
	<< lapconf.val_grid_point(1,0,nt-1,0) << endl ;
    ost << "lapse on the AH :    "
	<< lapse.val_grid_point(1,0,nt-1,0) << endl ;
    ost << "shift(1) on the AH : "
	<< shift(1).val_grid_point(1,0,nt-1,0) << endl ;
    ost << "shift(2) on the AH : "
	<< shift(2).val_grid_point(1,0,nt-1,0) << endl ;
    ost << "shift(3) on the AH : "
	<< shift(3).val_grid_point(1,0,nt-1,0) << endl ;
    ost << "confo on the AH :    "
	<< confo.val_grid_point(1,0,nt-1,0) << endl ;
    ost << "Gravitational mass : "
	<< mass_bh / msol << " M_sol" << endl ;
    ost << "Irreducible mass :   "
	<< mass_irr() / msol << " M_sol" << endl ;
    ost << "ADM mass :           "
	<< mass_adm() / msol << " M_sol" << endl ;
    ost << "Komar mass :         "
	<< mass_kom() / msol << " M_sol" << endl ;

    double irr_gm, adm_gm, kom_gm ;
    irr_gm = mass_irr() / mass_bh - 1. ;
    adm_gm = mass_adm() / mass_bh - 1. ;
    kom_gm = mass_kom() / mass_bh - 1. ;
    ost << "Diff. (Mirr-Mg)/Mg : " << irr_gm << endl ;
    ost << "Diff. (Madm-Mg)/Mg : " << adm_gm << endl ;
    ost << "Diff. (Mkom-Mg)/Mg : " << kom_gm << endl ;

    return ost ;

}
}
