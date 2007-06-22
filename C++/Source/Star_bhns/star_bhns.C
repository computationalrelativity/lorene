/*
 *  Methods of class Star_bhns
 *
 *    (see file star_bhns.h for documentation).
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

char star_bhns_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:30:10  k_taniguchi
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
#include "etoile.h"
#include "star.h"
#include "star_bhns.h"
#include "eos.h"
#include "unites.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ;

                    //---------------------//
                    //     Constructor     //
                    //---------------------//

// Standard constructor
// --------------------
Star_bhns::Star_bhns(Map& mp_i, int nzet_i, const Eos& eos_i, bool irrot_i)
      : Star(mp_i, nzet_i, eos_i),
	mp_aff(mp_i),
	irrotational(irrot_i),
	psi0(mp_i),
	d_psi(mp_i, COV, mp_i.get_bvect_cart()),
	wit_w(mp_i, CON, mp_i.get_bvect_cart()),
	loggam(mp_i),
	bsn(mp_i, CON, mp_i.get_bvect_cart()),
	gam(mp_i),
	gam0(mp_i),
	pot_centri(mp_i),
        lapse_auto(mp_i),
	lapse_comp(mp_i),
	lapse_tot(mp_i),
	d_lapse_auto(mp_i, COV, mp_i.get_bvect_cart()),
	d_lapse_comp(mp_i, COV, mp_i.get_bvect_cart()),
	shift_auto(mp_i, CON, mp_i.get_bvect_cart()),
	shift_comp(mp_i, CON, mp_i.get_bvect_cart()),
	shift_tot(mp_i, CON, mp_i.get_bvect_cart()),
	d_shift_auto(mp_i, 2, CON, mp_i.get_bvect_cart()),
	d_shift_comp(mp_i, 2, CON, mp_i.get_bvect_cart()),
	confo_auto(mp_i),
	confo_comp(mp_i),
	confo_tot(mp_i),
	d_confo_auto(mp_i, COV, mp_i.get_bvect_cart()),
	d_confo_comp(mp_i, COV, mp_i.get_bvect_cart()),
	psi4(mp_i),
	taij_auto(mp_i, CON, mp_i.get_bvect_cart()),
	taij_quad_auto(mp_i),
	flat(mp_i, mp_i.get_bvect_cart()),
	ssjm1_lapse(mp_i),
	ssjm1_confo(mp_i),
	ssjm1_khi(mp_i),
	ssjm1_wshift(mp_i, CON, mp_i.get_bvect_cart())
{

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;

    // Quantities defined on a spherical triad in star.C are put on
    // a cartesian one
    u_euler.change_triad(mp_i.get_bvect_cart()) ;

    // All quantities are initialized to zero
    psi0 = 0. ;
    psi0.std_spectral_base() ;
    d_psi.set_etat_zero() ;
    wit_w.set_etat_zero() ;
    loggam = 0. ;
    loggam.std_spectral_base() ;
    bsn.set_etat_zero() ;
    gam = 0. ;
    gam.std_spectral_base() ;
    gam0 = 0. ;
    gam0.std_spectral_base() ;
    pot_centri = 0. ;
    pot_centri.std_spectral_base() ;

    lapse_auto = 1. ;
    lapse_auto.std_spectral_base() ;
    lapse_comp = 0. ;
    lapse_comp.std_spectral_base() ;
    lapse_tot = 1. ;
    lapse_tot.std_spectral_base() ;
    d_lapse_auto.set_etat_zero() ;
    d_lapse_comp.set_etat_zero() ;
    shift_auto.set_etat_zero() ;
    shift_comp.set_etat_zero() ;
    shift_tot.set_etat_zero() ;
    d_shift_auto.set_etat_zero() ;
    d_shift_comp.set_etat_zero() ;
    confo_auto = 1. ;
    confo_auto.std_spectral_base() ;
    confo_comp = 0. ;
    confo_comp.std_spectral_base() ;
    confo_tot = 1. ;
    confo_tot.std_spectral_base() ;
    d_confo_auto.set_etat_zero() ;
    d_confo_comp.set_etat_zero() ;
    psi4 = 1. ;
    psi4.std_spectral_base() ;

    taij_auto.set_etat_zero() ;
    taij_quad_auto = 0. ;
    taij_quad_auto.std_spectral_base() ;

    ssjm1_lapse = 0. ;
    ssjm1_lapse.std_spectral_base() ;
    ssjm1_confo = 0. ;
    ssjm1_confo.std_spectral_base() ;
    ssjm1_khi = 0. ;
    ssjm1_khi.std_spectral_base() ;
    ssjm1_wshift.set_etat_zero() ;

}

// Copy constructor
// ----------------
Star_bhns::Star_bhns(const Star_bhns& star)
      : Star(star),
	mp_aff(star.mp_aff),
	irrotational(star.irrotational),
	psi0(star.psi0),
	d_psi(star.d_psi),
	wit_w(star.wit_w),
	loggam(star.loggam),
	bsn(star.bsn),
	gam(star.gam),
	gam0(star.gam0),
	pot_centri(star.pot_centri),
	lapse_auto(star.lapse_auto),
	lapse_comp(star.lapse_comp),
	lapse_tot(star.lapse_tot),
	d_lapse_auto(star.d_lapse_auto),
	d_lapse_comp(star.d_lapse_comp),
	shift_auto(star.shift_auto),
	shift_comp(star.shift_comp),
	shift_tot(star.shift_tot),
	d_shift_auto(star.d_shift_auto),
	d_shift_comp(star.d_shift_comp),
	confo_auto(star.confo_auto),
	confo_comp(star.confo_comp),
	confo_tot(star.confo_tot),
	d_confo_auto(star.d_confo_auto),
	d_confo_comp(star.d_confo_comp),
	psi4(star.psi4),
	taij_auto(star.taij_auto),
	taij_quad_auto(star.taij_quad_auto),
	flat(star.flat),
	ssjm1_lapse(star.ssjm1_lapse),
	ssjm1_confo(star.ssjm1_confo),
	ssjm1_khi(star.ssjm1_khi),
	ssjm1_wshift(star.ssjm1_wshift)
{
    set_der_0x0() ;
}

// Constructor from a file
// -----------------------
Star_bhns::Star_bhns(Map& mp_i, const Eos& eos_i, FILE* fich)
      : Star(mp_i, eos_i, fich),
	mp_aff(mp_i),
	psi0(mp_i),
	d_psi(mp_i, COV, mp_i.get_bvect_cart()),
	wit_w(mp_i, CON, mp_i.get_bvect_cart()),
	loggam(mp_i),
	bsn(mp_i, CON, mp_i.get_bvect_cart()),
	gam(mp_i),
	gam0(mp_i),
	pot_centri(mp_i),
	lapse_auto(mp_i, *(mp_i.get_mg()), fich),
	lapse_comp(mp_i),
	lapse_tot(mp_i),
	d_lapse_auto(mp_i, COV, mp_i.get_bvect_cart()),
	d_lapse_comp(mp_i, COV, mp_i.get_bvect_cart()),
	shift_auto(mp_i, mp_i.get_bvect_cart(), fich),
	shift_comp(mp_i, CON, mp_i.get_bvect_cart()),
	shift_tot(mp_i, CON, mp_i.get_bvect_cart()),
	d_shift_auto(mp_i, 2, CON, mp_i.get_bvect_cart()),
	d_shift_comp(mp_i, 2, CON, mp_i.get_bvect_cart()),
	confo_auto(mp_i, *(mp_i.get_mg()), fich),
	confo_comp(mp_i),
	confo_tot(mp_i),
	d_confo_auto(mp_i, COV, mp_i.get_bvect_cart()),
	d_confo_comp(mp_i, COV, mp_i.get_bvect_cart()),
	psi4(mp_i),
	taij_auto(mp_i, CON, mp_i.get_bvect_cart()),
	taij_quad_auto(mp_i),
	flat(mp_i, mp_i.get_bvect_cart()),
	ssjm1_lapse(mp_i, *(mp_i.get_mg()), fich),
	ssjm1_confo(mp_i, *(mp_i.get_mg()), fich),
	ssjm1_khi(mp_i, *(mp_i.get_mg()), fich),
	ssjm1_wshift(mp_i, mp_i.get_bvect_cart(), fich)
{

    // Star parameter
    fread(&irrotational, sizeof(bool), 1, fich) ;

    // Read of the saved fields
    // ------------------------

    if (irrotational) {
        Scalar gam_euler_file(mp, *(mp.get_mg()), fich) ;
	gam_euler = gam_euler_file ;

	Scalar psi0_file(mp, *(mp.get_mg()), fich) ;
	psi0 = psi0_file ;
    }

    // Quantities defined on a spherical triad in star.C are put on
    // a cartesian one
    u_euler.change_triad(mp_i.get_bvect_cart()) ;

    // All other quantities are initialized to zero
    // --------------------------------------------

    d_psi.set_etat_zero() ;
    wit_w.set_etat_zero() ;
    loggam = 0. ;
    loggam.std_spectral_base() ;
    bsn.set_etat_zero() ;
    gam = 0. ;
    gam.std_spectral_base() ;
    gam0 = 0. ;
    gam0.std_spectral_base() ;
    pot_centri = 0. ;
    pot_centri.std_spectral_base() ;

    lapse_comp = 0. ;
    lapse_comp.std_spectral_base() ;
    lapse_tot = 1. ;
    lapse_tot.std_spectral_base() ;
    d_lapse_auto.set_etat_zero() ;
    d_lapse_comp.set_etat_zero() ;
    shift_comp.set_etat_zero() ;
    shift_tot.set_etat_zero() ;
    d_shift_auto.set_etat_zero() ;
    d_shift_comp.set_etat_zero() ;
    confo_comp = 0. ;
    confo_comp.std_spectral_base() ;
    confo_tot = 1. ;
    confo_tot.std_spectral_base() ;
    d_confo_auto.set_etat_zero() ;
    d_confo_comp.set_etat_zero() ;
    psi4 = 1. ;
    psi4.std_spectral_base() ;
    taij_auto.set_etat_zero() ;
    taij_quad_auto = 0. ;
    taij_quad_auto.std_spectral_base() ;

    // Pointers of derived quantities initialized to zero
    // --------------------------------------------------
    set_der_0x0() ;

}


                    //--------------------//
                    //     Destructor     //
                    //--------------------//

Star_bhns::~Star_bhns()
{

    del_deriv() ;

}


                    //------------------------------------------//
                    //     Management of derived quantities     //
                    //------------------------------------------//

void Star_bhns::del_deriv() const {

    Star::del_deriv() ;

    if (p_mass_b_bhns != 0x0) delete p_mass_b_bhns ;
    if (p_mass_g_bhns != 0x0) delete p_mass_g_bhns ;

    set_der_0x0() ;

}

void Star_bhns::set_der_0x0() const {

    Star::set_der_0x0() ;

    p_mass_b_bhns = 0x0 ;
    p_mass_g_bhns = 0x0 ;

}


                    //--------------------//
                    //     Assignment     //
                    //--------------------//

// Assignment to another Star_bhns
// -------------------------------
void Star_bhns::operator=(const Star_bhns& star) {

    // Assignment of quantities common to the derived classes of Star
    Star::operator=(star) ;

    // Assignment of proper quantities of class Star_bhns
    mp_aff = star.mp_aff ;
    irrotational = star.irrotational ;
    psi0 = star.psi0 ;
    d_psi = star.d_psi ;
    wit_w = star.wit_w ;
    loggam = star.loggam ;
    bsn = star.bsn ;
    gam = star.gam ;
    gam0 = star.gam0 ;
    pot_centri = star.pot_centri ;
    lapse_auto = star.lapse_auto ;
    lapse_comp = star.lapse_comp ;
    lapse_tot = star.lapse_tot ;
    d_lapse_auto = star.d_lapse_auto ;
    d_lapse_comp = star.d_lapse_comp ;
    shift_auto = star.shift_auto ;
    shift_comp = star.shift_comp ;
    shift_tot = star.shift_tot ;
    d_shift_auto = star.d_shift_auto ;
    d_shift_comp = star.d_shift_comp ;
    confo_auto = star.confo_auto ;
    confo_comp = star.confo_comp ;
    confo_tot = star.confo_tot ;
    d_confo_auto = star.d_confo_auto ;
    d_confo_comp = star.d_confo_comp ;
    psi4 = star.psi4 ;
    taij_auto = star.taij_auto ;
    taij_quad_auto = star.taij_quad_auto ;
    flat = star.flat ;
    ssjm1_lapse = star.ssjm1_lapse ;
    ssjm1_confo = star.ssjm1_confo ;
    ssjm1_khi = star.ssjm1_khi ;
    ssjm1_wshift = star.ssjm1_wshift ;

    del_deriv() ;

}

Scalar& Star_bhns::set_pot_centri() {

    del_deriv() ;
    return pot_centri ;

}

Scalar& Star_bhns::set_lapse_auto() {

    del_deriv() ;
    return lapse_auto ;

}

Scalar& Star_bhns::set_lapse_comp() {

    del_deriv() ;
    return lapse_comp ;

}

Vector& Star_bhns::set_shift_auto() {

    del_deriv() ;
    return shift_auto ;

}

Vector& Star_bhns::set_shift_comp() {

    del_deriv() ;
    return shift_comp ;

}

Scalar& Star_bhns::set_confo_auto() {

    del_deriv() ;
    return confo_auto ;

}

Scalar& Star_bhns::set_confo_comp() {

    del_deriv() ;
    return confo_comp ;

}


                    //-----------------//
                    //     Outputs     //
                    //-----------------//

// Save in a file
// --------------
void Star_bhns::sauve(FILE* fich) const {

    Star::sauve(fich) ;

    lapse_auto.sauve(fich) ;
    shift_auto.sauve(fich) ;
    confo_auto.sauve(fich) ;

    ssjm1_lapse.sauve(fich) ;
    ssjm1_confo.sauve(fich) ;
    ssjm1_khi.sauve(fich) ;
    ssjm1_wshift.sauve(fich) ;

    fwrite(&irrotational, sizeof(bool), 1, fich) ;

    if (irrotational) {
        gam_euler.sauve(fich) ; // required to construct d_psi from psi0
	psi0.sauve(fich) ;
    }

}

// Printing
// --------

ostream& Star_bhns::operator>>(ostream& ost) const {

    using namespace Unites ;

    //    Star::operator>>(ost) ;

    ost << endl ;
    ost << "Neutron star in a BHNS binary" << endl ;
    ost << "-----------------------------" << endl ;

    ost << "Coordinate radius R_eq_tow :           "
	<< ray_eq_pi() / km << " [km]" << endl ;
    ost << "Coordinate radius R_eq_opp :           "
	<< ray_eq() / km << " [km]" << endl ;
    ost << "Coordinate radius R_eq_orb :           "
	<< ray_eq_pis2() / km << " [km]" << endl ;
    ost << "Coordinate radius R_pole :             "
	<< ray_pole() / km << " [km]" << endl ;
    ost << "Central enthalph H_ent :               "
	<< ent.val_grid_point(0,0,0,0) << endl ;
    ost << "Lapse function at the center of NS :   "
	<< lapse_tot.val_grid_point(0,0,0,0) << endl ;
    ost << "Conformal factor at the center of NS : "
	<< confo_tot.val_grid_point(0,0,0,0) << endl ;
    ost << "shift(1) at the center of NS :         "
	<< shift_tot(1).val_grid_point(0,0,0,0) << endl ;
    ost << "shift(2) at the center of NS :         "
	<< shift_tot(2).val_grid_point(0,0,0,0) << endl ;
    ost << "shift(3) at the center of NS :         "
	<< shift_tot(3).val_grid_point(0,0,0,0) << endl ;

    return ost ;

}


                    //--------------------------------//
                    //     Computational routines     //
                    //--------------------------------//

void Star_bhns::fait_d_psi_bhns() {

    if (!irrotational) {
        d_psi.set_etat_nondef() ;
	return ;
    }

    // Specific relativistic enthalpy          ---> hhh
    // ------------------------------

    Scalar hhh = exp(ent) ;  // = 1 at the Newtonian limit

    // Computation of W^i = h Gamma_n B^i / N
    // --------------------------------------

    Vector www = hhh * gam_euler * bsn * psi4 ;

    // Constant value of W^i at the center of the star
    // -----------------------------------------------

    Vector v_orb(mp, COV, mp.get_bvect_cart()) ;

    for (int i=1; i<=3; i++) {
        v_orb.set(i) = www(i).val_grid_point(0,0,0,0) ;
    }

    // Gradient of psi
    // ---------------

    Vector d_psi0(mp, COV, mp.get_bvect_cart()) ;
    d_psi0.set_etat_qcq() ;
    for (int i=1; i<=3; i++)
        d_psi0.set(i) = psi0.deriv(i) ;

    d_psi0.std_spectral_base() ;

    d_psi = d_psi0 + v_orb ;

    for (int i=1; i<=3; i++) {
        if (d_psi(i).get_etat() == ETATZERO)
	    d_psi.set(i).annule_hard() ;
    }

    // C^1 continuation of d_psi outside the star
    // (to ensure a smooth enthalpy field accross the stellar surface)
    // ---------------------------------------------------------------

    int nzm1 = mp.get_mg()->get_nzone() - 1 ;
    d_psi.annule(nzet, nzm1) ;
    for (int i=1; i<=3; i++) {
        Cmp d_psi_i (d_psi(i)) ;
	d_psi_i.va.set_base( d_psi0(i).get_spectral_va().base ) ;
	d_psi_i = raccord_c1(d_psi_i, nzet) ;
	d_psi.set(i) = d_psi_i ;
    }

}


void Star_bhns::relax_bhns(const Star_bhns& star_prev, double relax_ent,
			   double relax_met, int mer, int fmer_met) {

    double relax_ent_jm1 = 1. - relax_ent ;
    double relax_met_jm1 = 1. - relax_met ;

    ent = relax_ent * ent + relax_ent_jm1 * star_prev.ent ;

    if ( (mer != 0) && (mer % fmer_met == 0)) {

        lapse_auto = relax_met * lapse_auto
	    + relax_met_jm1 * star_prev.lapse_auto ;

	shift_auto = relax_met * shift_auto
	    + relax_met_jm1 * star_prev.shift_auto ;

	confo_auto = relax_met * confo_auto
	    + relax_met_jm1 * star_prev.confo_auto ;

    }

    del_deriv() ;

    equation_of_state() ;

}
