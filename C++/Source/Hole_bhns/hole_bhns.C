/*
 *  Methods of class Hole_bhns
 *
 *    (see file hole_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2006 Keisuke Taniguchi
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

char hole_bhns_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:24:16  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <math.h>

// Lorene headers
#include "hole_bhns.h"
#include "unites.h"


                    //----------------------//
                    //     Constructors     //
                    //----------------------//

// Standard constructor
// --------------------
Hole_bhns::Hole_bhns(Map& mp_i, bool kerrschild_i, bool bc_nd_i,
		     bool bc_fs_i, bool irrot_i, double massbh)
      : Black_hole(mp_i, kerrschild_i, massbh),
	bc_lapse_nd(bc_nd_i),
	bc_lapse_fs(bc_fs_i),
	irrotational(irrot_i),
	lapse_auto_rs(mp_i),
	lapse_auto_bh(mp_i),
	lapse_auto(mp_i),
	lapse_comp(mp_i),
	lapse_tot(mp_i),
	d_lapse_auto_rs(mp_i, COV, mp_i.get_bvect_cart()),
	d_lapse_auto_bh(mp_i, COV, mp_i.get_bvect_cart()),
	d_lapse_auto(mp_i, COV, mp_i.get_bvect_cart()),
	d_lapse_comp(mp_i, COV, mp_i.get_bvect_cart()),
	shift_auto_rs(mp_i, CON, mp_i.get_bvect_cart()),
	shift_auto_bh(mp_i, CON, mp_i.get_bvect_cart()),
	shift_auto(mp_i, CON, mp_i.get_bvect_cart()),
	shift_comp(mp_i, CON, mp_i.get_bvect_cart()),
	shift_tot(mp_i, CON, mp_i.get_bvect_cart()),
	d_shift_auto_rs(mp_i, 2, CON, mp_i.get_bvect_cart()),
	d_shift_auto_bh(mp_i, 2, CON, mp_i.get_bvect_cart()),
	d_shift_auto(mp_i, 2, CON, mp_i.get_bvect_cart()),
	d_shift_comp(mp_i, 2, CON, mp_i.get_bvect_cart()),
	confo_auto_rs(mp_i),
	confo_auto_bh(mp_i),
	confo_auto(mp_i),
	confo_comp(mp_i),
	confo_tot(mp_i),
	d_confo_auto_rs(mp_i, COV, mp_i.get_bvect_cart()),
	d_confo_auto_bh(mp_i, COV, mp_i.get_bvect_cart()),
	d_confo_auto(mp_i, COV, mp_i.get_bvect_cart()),
	d_confo_comp(mp_i, COV, mp_i.get_bvect_cart()),
	taij_tot_rs(mp_i, CON, mp_i.get_bvect_cart()),
	taij_tot_rot(mp_i, CON, mp_i.get_bvect_cart()),
	taij_tot_bh(mp_i, CON, mp_i.get_bvect_cart()),
	taij_tot(mp_i, CON, mp_i.get_bvect_cart()),
	taij_auto_rs(mp_i, CON, mp_i.get_bvect_cart()),
	taij_auto(mp_i, CON, mp_i.get_bvect_cart()),
	taij_comp(mp_i, CON, mp_i.get_bvect_cart()),
	taij_quad_tot_rs(mp_i),
	taij_quad_tot_rot(mp_i),
	taij_quad_tot_bh(mp_i),
	taij_quad_tot(mp_i),
	taij_quad_auto(mp_i),
	taij_quad_comp(mp_i)
{

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;

    // The metric quantities are initialized to the flat one or zero
    lapse_auto_rs = 0. ;
    lapse_auto_rs.std_spectral_base() ;
    lapse_auto_bh = 1. ;
    lapse_auto_bh.std_spectral_base() ;
    lapse_auto = 1. ;
    lapse_auto.std_spectral_base() ;
    lapse_comp = 0. ;
    lapse_comp.std_spectral_base() ;
    lapse_tot = 1. ;
    lapse_tot.std_spectral_base() ;

    d_lapse_auto_rs.set_etat_zero() ;
    d_lapse_auto_bh.set_etat_zero() ;
    d_lapse_auto.set_etat_zero() ;
    d_lapse_comp.set_etat_zero() ;

    shift_auto_rs.set_etat_zero() ;
    shift_auto_bh.set_etat_zero() ;
    shift_auto.set_etat_zero() ;
    shift_comp.set_etat_zero() ;
    shift_tot.set_etat_zero() ;

    d_shift_auto_rs.set_etat_zero() ;
    d_shift_auto_bh.set_etat_zero() ;
    d_shift_auto.set_etat_zero() ;
    d_shift_comp.set_etat_zero() ;

    confo_auto_rs = 0. ;
    confo_auto_rs.std_spectral_base() ;
    confo_auto_bh = 1. ;
    confo_auto_bh.std_spectral_base() ;
    confo_auto = 1. ;
    confo_auto.std_spectral_base() ;
    confo_comp = 0. ;
    confo_comp.std_spectral_base() ;
    confo_tot = 1. ;
    confo_tot.std_spectral_base() ;

    d_confo_auto_rs.set_etat_zero() ;
    d_confo_auto_bh.set_etat_zero() ;
    d_confo_auto.set_etat_zero() ;
    d_confo_comp.set_etat_zero() ;

    taij_tot_rs.set_etat_zero() ;
    taij_tot_rot.set_etat_zero() ;
    taij_tot_bh.set_etat_zero() ;
    taij_tot.set_etat_zero() ;
    taij_auto_rs.set_etat_zero() ;
    taij_auto.set_etat_zero() ;
    taij_comp.set_etat_zero() ;

    taij_quad_tot_rs = 0. ;
    taij_quad_tot_rs.std_spectral_base() ;
    taij_quad_tot_rot = 0. ;
    taij_quad_tot_rot.std_spectral_base() ;
    taij_quad_tot_bh = 0. ;
    taij_quad_tot_bh.std_spectral_base() ;
    taij_quad_tot = 0. ;
    taij_quad_tot.std_spectral_base() ;
    taij_quad_auto = 0. ;
    taij_quad_auto.std_spectral_base() ;
    taij_quad_comp = 0. ;
    taij_quad_comp.std_spectral_base() ;

}

// Copy constructor
// ----------------
Hole_bhns::Hole_bhns(const Hole_bhns& hole)
      : Black_hole(hole),
	bc_lapse_nd(hole.bc_lapse_nd),
	bc_lapse_fs(hole.bc_lapse_fs),
	irrotational(hole.irrotational),
	lapse_auto_rs(hole.lapse_auto_rs),
	lapse_auto_bh(hole.lapse_auto_bh),
	lapse_auto(hole.lapse_auto),
	lapse_comp(hole.lapse_comp),
	lapse_tot(hole.lapse_tot),
	d_lapse_auto_rs(hole.d_lapse_auto_rs),
	d_lapse_auto_bh(hole.d_lapse_auto_bh),
	d_lapse_auto(hole.d_lapse_auto),
	d_lapse_comp(hole.d_lapse_comp),
	shift_auto_rs(hole.shift_auto_rs),
	shift_auto_bh(hole.shift_auto_bh),
	shift_auto(hole.shift_auto),
	shift_comp(hole.shift_comp),
	shift_tot(hole.shift_tot),
	d_shift_auto_rs(hole.d_shift_auto_rs),
	d_shift_auto_bh(hole.d_shift_auto_bh),
	d_shift_auto(hole.d_shift_auto),
	d_shift_comp(hole.d_shift_comp),
	confo_auto_rs(hole.confo_auto_rs),
	confo_auto_bh(hole.confo_auto_bh),
	confo_auto(hole.confo_auto),
	confo_comp(hole.confo_comp),
	confo_tot(hole.confo_tot),
	d_confo_auto_rs(hole.d_confo_auto_rs),
	d_confo_auto_bh(hole.d_confo_auto_bh),
	d_confo_auto(hole.d_confo_auto),
	d_confo_comp(hole.d_confo_comp),
	taij_tot_rs(hole.taij_tot_rs),
	taij_tot_rot(hole.taij_tot_rot),
	taij_tot_bh(hole.taij_tot_bh),
	taij_tot(hole.taij_tot),
	taij_auto_rs(hole.taij_auto_rs),
	taij_auto(hole.taij_auto),
	taij_comp(hole.taij_comp),
	taij_quad_tot_rs(hole.taij_quad_tot_rs),
	taij_quad_tot_rot(hole.taij_quad_tot_rot),
	taij_quad_tot_bh(hole.taij_quad_tot_bh),
	taij_quad_tot(hole.taij_quad_tot),
	taij_quad_auto(hole.taij_quad_auto),
	taij_quad_comp(hole.taij_quad_comp)
{
    set_der_0x0() ;
}

// Constructor from a file
// -----------------------
Hole_bhns::Hole_bhns(Map& mp_i, FILE* fich)
      : Black_hole(mp_i, fich),
	lapse_auto_rs(mp_i, *(mp_i.get_mg()), fich),
	lapse_auto_bh(mp_i),
	lapse_auto(mp_i),
	lapse_comp(mp_i),
	lapse_tot(mp_i),
	d_lapse_auto_rs(mp_i, COV, mp_i.get_bvect_cart()),
	d_lapse_auto_bh(mp_i, COV, mp_i.get_bvect_cart()),
	d_lapse_auto(mp_i, COV, mp_i.get_bvect_cart()),
	d_lapse_comp(mp_i, COV, mp_i.get_bvect_cart()),
	shift_auto_rs(mp_i, mp_i.get_bvect_cart(), fich),
	shift_auto_bh(mp_i, CON, mp_i.get_bvect_cart()),
	shift_auto(mp_i, CON, mp_i.get_bvect_cart()),
	shift_comp(mp_i, CON, mp_i.get_bvect_cart()),
	shift_tot(mp_i, CON, mp_i.get_bvect_cart()),
	d_shift_auto_rs(mp_i, 2, CON, mp_i.get_bvect_cart()),
	d_shift_auto_bh(mp_i, 2, CON, mp_i.get_bvect_cart()),
	d_shift_auto(mp_i, 2, CON, mp_i.get_bvect_cart()),
	d_shift_comp(mp_i, 2, CON, mp_i.get_bvect_cart()),
	confo_auto_rs(mp_i, *(mp_i.get_mg()), fich),
	confo_auto_bh(mp_i),
	confo_auto(mp_i),
	confo_comp(mp_i),
	confo_tot(mp_i),
	d_confo_auto_rs(mp_i, COV, mp_i.get_bvect_cart()),
	d_confo_auto_bh(mp_i, COV, mp_i.get_bvect_cart()),
	d_confo_auto(mp_i, COV, mp_i.get_bvect_cart()),
	d_confo_comp(mp_i, COV, mp_i.get_bvect_cart()),
	taij_tot_rs(mp_i, CON, mp_i.get_bvect_cart()),
	taij_tot_rot(mp_i, CON, mp_i.get_bvect_cart()),
	taij_tot_bh(mp_i, CON, mp_i.get_bvect_cart()),
	taij_tot(mp_i, CON, mp_i.get_bvect_cart()),
	taij_auto_rs(mp_i, CON, mp_i.get_bvect_cart()),
	taij_auto(mp_i, CON, mp_i.get_bvect_cart()),
	taij_comp(mp_i, CON, mp_i.get_bvect_cart()),
	taij_quad_tot_rs(mp_i),
	taij_quad_tot_rot(mp_i),
	taij_quad_tot_bh(mp_i),
	taij_quad_tot(mp_i),
	taij_quad_auto(mp_i),
	taij_quad_comp(mp_i)
{

    fread(&bc_lapse_nd, sizeof(bool), 1, fich) ;
    fread(&bc_lapse_fs, sizeof(bool), 1, fich) ;
    fread(&irrotational, sizeof(bool), 1, fich) ;

    // All other quantities are initialized to zero
    // --------------------------------------------

    lapse_auto_bh = 1. ;
    lapse_auto_bh.std_spectral_base() ;
    lapse_auto = 1. ;
    lapse_auto.std_spectral_base() ;
    lapse_comp = 0. ;
    lapse_comp.std_spectral_base() ;
    lapse_tot = 1. ;
    lapse_tot.std_spectral_base() ;

    d_lapse_auto_rs.set_etat_zero() ;
    d_lapse_auto_bh.set_etat_zero() ;
    d_lapse_auto.set_etat_zero() ;
    d_lapse_comp.set_etat_zero() ;

    shift_auto_bh.set_etat_zero() ;
    shift_auto.set_etat_zero() ;
    shift_comp.set_etat_zero() ;
    shift_tot.set_etat_zero() ;
    d_shift_auto_rs.set_etat_zero() ;
    d_shift_auto_bh.set_etat_zero() ;
    d_shift_auto.set_etat_zero() ;
    d_shift_comp.set_etat_zero() ;

    confo_auto_bh = 1. ;
    confo_auto_bh.std_spectral_base() ;
    confo_auto = 1. ;
    confo_auto.std_spectral_base() ;
    confo_comp = 0. ;
    confo_comp.std_spectral_base() ;
    confo_tot = 1. ;
    confo_tot.std_spectral_base() ;

    d_confo_auto_rs.set_etat_zero() ;
    d_confo_auto_bh.set_etat_zero() ;
    d_confo_auto.set_etat_zero() ;
    d_confo_comp.set_etat_zero() ;

    taij_tot_rs.set_etat_zero() ;
    taij_tot_rot.set_etat_zero() ;
    taij_tot_bh.set_etat_zero() ;
    taij_tot.set_etat_zero() ;
    taij_auto_rs.set_etat_zero() ;
    taij_auto.set_etat_zero() ;
    taij_comp.set_etat_zero() ;
    taij_quad_tot_rs = 0. ;
    taij_quad_tot_rs.std_spectral_base() ;
    taij_quad_tot_rot = 0. ;
    taij_quad_tot_rot.std_spectral_base() ;
    taij_quad_tot_bh = 0. ;
    taij_quad_tot_bh.std_spectral_base() ;
    taij_quad_tot = 0. ;
    taij_quad_tot.std_spectral_base() ;
    taij_quad_auto = 0. ;
    taij_quad_auto.std_spectral_base() ;
    taij_quad_comp = 0. ;
    taij_quad_comp.std_spectral_base() ;

    // Pointers of derived quantities initialized to zero
    // --------------------------------------------------
    set_der_0x0() ;

}


                    //--------------------//
                    //     Destructor     //
                    //--------------------//

Hole_bhns::~Hole_bhns()
{

    del_deriv() ;

}


                    //------------------------------------------//
                    //     Management of derived quantities     //
                    //------------------------------------------//

void Hole_bhns::del_deriv() const {

    Black_hole::del_deriv() ;

    if (p_mass_irr_bhns != 0x0) delete p_mass_irr_bhns ;

    set_der_0x0() ;

}

void Hole_bhns::set_der_0x0() const {

    Black_hole::set_der_0x0() ;

    p_mass_irr_bhns = 0x0 ;

}


                    //--------------------//
                    //     Assignment     //
                    //--------------------//

// Assignment to another Hole_bhns
// -------------------------------
void Hole_bhns::operator=(const Hole_bhns& hole) {

    // Assignment of quantities common to the derived classes of Black_hole
    Black_hole::operator=(hole) ;

    // Assignment of proper quantities of class Hole_bhns
    bc_lapse_nd = hole.bc_lapse_nd ;
    bc_lapse_fs = hole.bc_lapse_fs ;
    irrotational = hole.irrotational ;
    lapse_auto_rs = hole.lapse_auto_rs ;
    lapse_auto_bh = hole.lapse_auto_bh ;
    lapse_auto = hole.lapse_auto ;
    lapse_comp = hole.lapse_comp ;
    lapse_tot = hole.lapse_tot ;
    d_lapse_auto_rs = hole.d_lapse_auto_rs ;
    d_lapse_auto_bh = hole.d_lapse_auto_bh ;
    d_lapse_auto = hole.d_lapse_auto ;
    d_lapse_comp = hole.d_lapse_comp ;
    shift_auto_rs = hole.shift_auto_rs ;
    shift_auto_bh = hole.shift_auto_bh ;
    shift_auto = hole.shift_auto ;
    shift_comp = hole.shift_comp ;
    shift_tot = hole.shift_tot ;
    d_shift_auto_rs = hole.d_shift_auto_rs ;
    d_shift_auto_bh = hole.d_shift_auto_bh ;
    d_shift_auto = hole.d_shift_auto ;
    d_shift_comp = hole.d_shift_comp ;
    confo_auto_rs = hole.confo_auto_rs ;
    confo_auto_bh = hole.confo_auto_bh ;
    confo_auto = hole.confo_auto ;
    confo_comp = hole.confo_comp ;
    confo_tot = hole.confo_tot ;
    d_confo_auto_rs = hole.d_confo_auto_rs ;
    d_confo_auto_bh = hole.d_confo_auto_bh ;
    d_confo_auto = hole.d_confo_auto ;
    d_confo_comp = hole.d_confo_comp ;
    taij_tot_rs = hole.taij_tot_rs ;
    taij_tot_rot = hole.taij_tot_rot ;
    taij_tot_bh = hole.taij_tot_bh ;
    taij_tot = hole.taij_tot ;
    taij_auto_rs = hole.taij_auto_rs ;
    taij_auto = hole.taij_auto ;
    taij_comp = hole.taij_comp ;
    taij_quad_tot_rs = hole.taij_quad_tot_rs ;
    taij_quad_tot_rot = hole.taij_quad_tot_rot ;
    taij_quad_tot_bh = hole.taij_quad_tot_bh ;
    taij_quad_tot = hole.taij_quad_tot ;
    taij_quad_auto = hole.taij_quad_auto ;
    taij_quad_comp = hole.taij_quad_comp ;

    del_deriv() ;

}


Scalar& Hole_bhns::set_lapse_auto_bh() {

    del_deriv() ;
    return lapse_auto_bh ;

}

Scalar& Hole_bhns::set_lapse_auto_rs() {

    del_deriv() ;
    return lapse_auto_rs ;

}

Scalar& Hole_bhns::set_lapse_auto() {

    del_deriv() ;
    return lapse_auto ;

}

Scalar& Hole_bhns::set_lapse_comp() {

    del_deriv() ;
    return lapse_comp ;

}

Scalar& Hole_bhns::set_lapse_tot() {

    del_deriv() ;
    return lapse_tot ;

}

Vector& Hole_bhns::set_shift_auto_bh() {

    del_deriv() ;
    return shift_auto_bh ;

}

Vector& Hole_bhns::set_shift_auto_rs() {

    del_deriv() ;
    return shift_auto_rs ;

}

Vector& Hole_bhns::set_shift_auto() {

    del_deriv() ;
    return shift_auto ;

}

Vector& Hole_bhns::set_shift_comp() {

    del_deriv() ;
    return shift_comp ;

}

Vector& Hole_bhns::set_shift_tot() {

    del_deriv() ;
    return shift_tot ;

}

Scalar& Hole_bhns::set_confo_auto_rs() {

    del_deriv() ;
    return confo_auto_rs ;

}

Scalar& Hole_bhns::set_confo_auto_bh() {

    del_deriv() ;
    return confo_auto_bh ;

}

Scalar& Hole_bhns::set_confo_auto() {

    del_deriv() ;
    return confo_auto ;

}

Scalar& Hole_bhns::set_confo_comp() {

    del_deriv() ;
    return confo_comp ;

}

Scalar& Hole_bhns::set_confo_tot() {

    del_deriv() ;
    return confo_tot ;

}


                    //-----------------//
                    //     Outputs     //
                    //-----------------//

// Save in a file
// --------------
void Hole_bhns::sauve(FILE* fich) const {

    Black_hole::sauve(fich) ;

    lapse_auto_rs.sauve(fich) ;
    shift_auto_rs.sauve(fich) ;
    confo_auto_rs.sauve(fich) ;

    fwrite(&bc_lapse_nd, sizeof(bool), 1, fich) ;
    fwrite(&bc_lapse_fs, sizeof(bool), 1, fich) ;
    fwrite(&irrotational, sizeof(bool), 1, fich) ;

}

// Printing
// --------
ostream& Hole_bhns::operator>>(ostream& ost) const {

    using namespace Unites ;

    //    Black_hole::operator>>(ost) ;

    ost << endl ;
    ost << "Black hole in a BHNS binary" << endl ;
    ost << "---------------------------" << endl ;

    int nt = mp.get_mg()->get_nt(1) ;

    ost << "Irreducible mass of BH :         "
	<< mass_irr_bhns() / msol << " [Mo]" << endl ;
    ost << "Mass in the background :         "
	<< mass_bh / msol << " [Mo]" << endl ;
    ost << "Radius of the apparent horison : "
	<< rad_ah() / km << " [km]" << endl ;
    ost << "Lapse function on the AH :       "
	<< lapse_tot.val_grid_point(1,0,nt-1,0) << endl ;
    ost << "Conformal factor on the AH :     "
	<< confo_tot.val_grid_point(1,0,nt-1,0) << endl ;
    ost << "shift(1) on the AH :             "
	<< shift_tot(1).val_grid_point(1,0,nt-1,0) << endl ;
    ost << "shift(2) on the AH :             "
	<< shift_tot(2).val_grid_point(1,0,nt-1,0) << endl ;
    ost << "shift(3) on the AH :             "
	<< shift_tot(3).val_grid_point(1,0,nt-1,0) << endl ;

    return ost ;

}

                    //--------------------------------//
                    //     Computational routines     //
                    //--------------------------------//

void Hole_bhns::relax_bhns(const Hole_bhns& hole_prev,
			   double relax_met, int mer, int fmer_met) {

    double relax_met_jm1 = 1. - relax_met ;

    if ( (mer != 0) && (mer % fmer_met == 0)) {

        lapse_auto_rs = relax_met * lapse_auto_rs
	    + relax_met_jm1 * hole_prev.lapse_auto_rs ;

	shift_auto_rs = relax_met * shift_auto_rs
	    + relax_met_jm1 * hole_prev.shift_auto_rs ;

	confo_auto_rs = relax_met * confo_auto_rs
	    + relax_met_jm1 * hole_prev.confo_auto_rs ;

    }

    del_deriv() ;

}

