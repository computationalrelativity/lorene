/*
 *  Methods for the class Et_bin_nsbh
 *
 *    (see file et_bin_nsbh.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Keisuke Taniguchi
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

char et_bin_nsbh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2003/12/16 05:30:18  k_taniguchi
 * Some changes in the constructor from file and the saving to file.
 *
 * Revision 1.4  2003/11/28 04:30:22  k_taniguchi
 * Add the output operator.
 *
 * Revision 1.3  2003/10/24 12:30:58  k_taniguchi
 * Add some derivative terms
 *
 * Revision 1.2  2003/10/22 08:12:46  k_taniguchi
 * Supress some terms in Et_bin_nsbh::sauve.
 *
 * Revision 1.1  2003/10/21 11:50:38  k_taniguchi
 * Methods for class Et_bin_nsbh
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "et_bin_nsbh.h"
#include "etoile.h"
#include "eos.h"
#include "utilitaires.h"
#include "param.h"

                            //--------------//
                            // Constructors //
                            //--------------//

// Standard constructor
// --------------------
Et_bin_nsbh::Et_bin_nsbh(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
			 bool irrot, const Base_vect& ref_triad_i)
    : Etoile_bin(mp_i, nzet_i, relat, eos_i, irrot, ref_triad_i),
      n_auto(mp_i),
      n_comp(mp_i),
      d_n_auto(mp_i, 1, COV, ref_triad_i),
      d_n_comp(mp_i, 1, COV, ref_triad_i),
      confpsi(mp_i),
      confpsi_auto(mp_i),
      confpsi_comp(mp_i),
      d_confpsi_auto(mp_i, 1, COV, ref_triad_i),
      d_confpsi_comp(mp_i, 1, COV, ref_triad_i),
      taij_auto(mp_i, 2, CON, ref_triad_i),
      taij_comp(mp_i, 2, CON, ref_triad_i),
      taij_tot(mp_i, 2, CON, ref_triad_i),
      tkij_auto(mp_i, 2, CON, ref_triad_i),
      tkij_tot(mp_i, 2, CON, ref_triad_i),
      ssjm1_lapse(mp_i),
      ssjm1_confpsi(mp_i) {

    // Pointers of derived quantities initialized to zero :
    set_der_0x0() ;

    // The metric is initialized to the flat one :
    n_auto = 0.5 ;
    n_auto.set_std_base() ;
    n_comp = 0.5 ;
    n_comp.set_std_base() ;
    d_n_auto = 0 ;
    d_n_comp = 0 ;
    confpsi = 1 ;
    confpsi.set_std_base() ;
    confpsi_auto = 0.5 ;
    confpsi_auto.set_std_base() ;
    confpsi_comp = 0.5 ;
    confpsi_comp.set_std_base() ;
    d_confpsi_auto = 0 ;
    d_confpsi_comp = 0 ;

    taij_auto.set_etat_zero() ;
    taij_comp.set_etat_zero() ;
    taij_tot.set_etat_zero() ;
    tkij_auto.set_etat_zero() ;
    tkij_tot.set_etat_zero() ;

    ssjm1_lapse.set_etat_qcq() ;
    ssjm1_lapse = 0.5 ;
    ssjm1_confpsi.set_etat_qcq() ;
    ssjm1_confpsi = 0.5 ;

}

// Copy constructor
// ----------------
Et_bin_nsbh::Et_bin_nsbh(const Et_bin_nsbh& et)
    : Etoile_bin(et),
      n_auto(et.n_auto),
      n_comp(et.n_comp),
      d_n_auto(et.d_n_auto),
      d_n_comp(et.d_n_comp),
      confpsi(et.confpsi),
      confpsi_auto(et.confpsi_auto),
      confpsi_comp(et.confpsi_comp),
      d_confpsi_auto(et.d_confpsi_auto),
      d_confpsi_comp(et.d_confpsi_comp),
      taij_auto(et.taij_auto),
      taij_comp(et.taij_comp),
      taij_tot(et.taij_tot),
      tkij_auto(et.tkij_auto),
      tkij_tot(et.tkij_tot),
      ssjm1_lapse(et.ssjm1_lapse),
      ssjm1_confpsi(et.ssjm1_confpsi) {

    set_der_0x0() ;

}

// Constructor from a file
// -----------------------
Et_bin_nsbh::Et_bin_nsbh(Map& mp_i, const Eos& eos_i,
			 const Base_vect& ref_triad_i, FILE* fich)
    : Etoile_bin(mp_i, eos_i, ref_triad_i, fich),
      n_auto(mp_i),
      n_comp(mp_i),
      d_n_auto(mp_i, 1, COV, ref_triad_i),
      d_n_comp(mp_i, 1, COV, ref_triad_i),
      confpsi(mp_i),
      confpsi_auto(mp_i),
      confpsi_comp(mp_i),
      d_confpsi_auto(mp_i, 1, COV, ref_triad_i),
      d_confpsi_comp(mp_i, 1, COV, ref_triad_i),
      taij_auto(mp_i, 2, CON, ref_triad_i),
      taij_comp(mp_i, 2, CON, ref_triad_i),
      taij_tot(mp_i, 2, CON, ref_triad_i),
      tkij_auto(mp_i, 2, CON, ref_triad_i),
      tkij_tot(mp_i, 2, CON, ref_triad_i),
      ssjm1_lapse(mp_i),
      ssjm1_confpsi(mp_i) {

    // Read of the saved fields :
    Cmp ssjm1_lapse_file(mp_i, *(mp_i.get_mg()), fich) ;
    ssjm1_lapse = ssjm1_lapse_file ;

    Cmp ssjm1_confpsi_file(mp_i, *(mp_i.get_mg()), fich) ;
    ssjm1_confpsi = ssjm1_confpsi_file ;

    // Construct from data in Etoile_bin
    n_auto = exp(logn_auto) ;
    n_auto.set_std_base() ;
    confpsi_auto = exp(0.5*(beta_auto-logn_auto)) ;
    confpsi_auto.set_std_base() ;

    // All other fields are initialized to zero or some constants :
    // ----------------------------------------------------------
    n_comp = 0.5 ;
    d_n_auto = 0 ;
    d_n_comp = 0 ;
    confpsi_comp = 0.5 ;
    d_confpsi_auto = 0 ;
    d_confpsi_comp = 0 ;
    taij_auto.set_etat_zero() ;
    taij_comp.set_etat_zero() ;
    taij_tot.set_etat_zero() ;
    tkij_auto.set_etat_zero() ;
    tkij_tot.set_etat_zero() ;

    // Pointers of derived quantities initialized to zero
    // --------------------------------------------------
    set_der_0x0() ;

}

                            //------------//
                            // Destructor //
                            //------------//

Et_bin_nsbh::~Et_bin_nsbh(){

    del_deriv() ;

}

                            //--------------//
                            //  Assignment  //
                            //--------------//

// Assignment to another Et_bin_nsbh
// ---------------------------------
void Et_bin_nsbh::operator=(const Et_bin_nsbh& et) {

    // Assignment of quantities common to the derived classes of Etoile_bin
    Etoile_bin::operator=(et) ;

    n_auto = et.n_auto ;
    n_comp = et.n_comp ;
    d_n_auto = et.d_n_auto ;
    d_n_comp = et.d_n_comp ;
    confpsi = et.confpsi ;
    confpsi_auto = et.confpsi_auto ;
    confpsi_comp = et.confpsi_comp ;
    d_confpsi_auto = et.d_confpsi_auto ;
    d_confpsi_comp = et.d_confpsi_comp ;
    taij_auto = et.taij_auto ;
    taij_comp = et.taij_comp ;
    taij_tot = et.taij_tot ;
    tkij_auto = et.tkij_auto ;
    tkij_tot = et.tkij_tot ;
    ssjm1_lapse = et.ssjm1_lapse ;
    ssjm1_confpsi = et.ssjm1_confpsi ;

    del_deriv() ;  // Deletes all derived quantities
}

Tenseur& Et_bin_nsbh::set_n_comp() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return n_comp ;

}

Tenseur& Et_bin_nsbh::set_confpsi_comp() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return confpsi_comp ;

}

                            //--------------//
                            //    Outputs   //
                            //--------------//

// Save in a file
// --------------
void Et_bin_nsbh::sauve(FILE* fich) const {

    Etoile_bin::sauve(fich) ;

    ssjm1_lapse.sauve(fich) ;
    ssjm1_confpsi.sauve(fich) ;
}


// Printing
// --------

ostream& Et_bin_nsbh::operator>>(ostream& ost) const {

    #include "unites.h"
    // To avoid some compilation warnings
    if (&ost == 0x0) {
      cout << qpig << msol << f_unit << mevpfm3 << endl ;
    }

    Etoile::operator>>(ost) ;

    ost << endl ;
    ost << "Neutron star in a binary system" << endl ;
    ost << "-------------------------------" << endl ;

    if (irrotational) {
      ost << "irrotational configuration" << endl ;
    }
    else {
      ost << "corotating configuration" << endl ;
    }

    ost << "Absolute abscidia of the stellar center: " <<
      mp.get_ori_x() / km << " km" << endl ;

    ost << "Absolute abscidia of the barycenter of the baryon density : " <<
      xa_barycenter() / km << " km" << endl ;

    double r_0 = 0.5 * ( ray_eq() + ray_eq_pi() ) ;
    double d_ns = fabs( mp.get_ori_x() ) + ray_eq_pi() - r_0 ;
    double d_tilde = 2 * d_ns / r_0 ;

    ost << "d_tilde : " << d_tilde << endl ;

    ost << "Orientation with respect to the absolute frame : " <<
      mp.get_rot_phi() << " rad" << endl ;

    ost << "Central value of gam_euler : "
	<< gam_euler()(0, 0, 0, 0)  << endl ;

    ost << "Central u_euler (U^X, U^Y, U^Z) [c] : "
	<< u_euler(0)(0, 0, 0, 0) << "  "
	<< u_euler(1)(0, 0, 0, 0) << "  "
	<< u_euler(2)(0, 0, 0, 0) << endl ;

    if (irrotational) {
      ost << "Central d_psi (X, Y, Z) [c] :         "
	  << d_psi(0)(0, 0, 0, 0) << "  "
	  << d_psi(1)(0, 0, 0, 0) << "  "
	  << d_psi(2)(0, 0, 0, 0) << endl ;

      ost << "Central vel. / co-orb. (W^X, W^Y, W^Z) [c] : "
	  << wit_w(0)(0, 0, 0, 0) << "  "
	  << wit_w(1)(0, 0, 0, 0) << "  "
	  << wit_w(2)(0, 0, 0, 0) << endl ;

      ost << "Max vel. / co-orb. (W^X, W^Y, W^Z) [c] : "
	  << max(max(wit_w(0))) << "  "
	  << max(max(wit_w(1))) << "  "
	  << max(max(wit_w(2))) << endl ;

      ost << "Min vel. / co-orb. (W^X, W^Y, W^Z) [c] : "
	  << min(min(wit_w(0))) << "  "
	  << min(min(wit_w(1))) << "  "
	  << min(min(wit_w(2))) << endl ;

      double r_surf = mp.val_r(0,1.,M_PI/4,M_PI/4) ;

      ost << "Velocity at (r_surf,pi/4,pi/4) / co-orb. [c] : "
	  << wit_w(0).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	  << wit_w(1).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	  << wit_w(2).val_point(r_surf,M_PI/4,M_PI/4) << endl ;

      ost << "Central value of loggam : "
	  << loggam()(0, 0, 0, 0)  << endl ;
    }

    ost << "Central value of lapse(N) auto :             "
	<< n_auto()(0, 0, 0, 0) << endl ;

    ost << "Central value of confpsi auto :              "
	<< confpsi_auto()(0, 0, 0, 0) << endl ;

    ost << "Central value of shift (N^X, N^Y, N^Z) [c] : "
	<< shift(0)(0, 0, 0, 0) << "  "
	<< shift(1)(0, 0, 0, 0) << "  "
	<< shift(2)(0, 0, 0, 0) << endl ;

    ost << "  ... shift_auto part of it [c] :            "
	<< shift_auto(0)(0, 0, 0, 0) << "  "
	<< shift_auto(1)(0, 0, 0, 0) << "  "
	<< shift_auto(2)(0, 0, 0, 0) << endl ;

    ost << "  ... shift_comp part of it [c] :            "
	<< shift_comp(0)(0, 0, 0, 0) << "  "
	<< shift_comp(1)(0, 0, 0, 0) << "  "
	<< shift_comp(2)(0, 0, 0, 0) << endl ;

    ost << "  ... w_shift (NB: components in the star Cartesian frame) [c] :  "
	<< endl
	<< w_shift(0)(0, 0, 0, 0) << "  "
	<< w_shift(1)(0, 0, 0, 0) << "  "
	<< w_shift(2)(0, 0, 0, 0) << endl ;

    ost << "Central value of khi_shift [km c] : "
        << khi_shift()(0, 0, 0, 0) / km << endl ;

    ost << endl << "Central value of (B^X, B^Y, B^Z)/N [c] : "
	<< bsn(0)(0, 0, 0, 0) << "  "
	<< bsn(1)(0, 0, 0, 0) << "  "
	<< bsn(2)(0, 0, 0, 0) << endl ;

    ost << endl <<
      "Central (d/dX,d/dY,d/dZ)(logn_auto) [km^{-1}] : "
	<< d_n_auto(0)(0, 0, 0, 0) * km << "  "
	<< d_n_auto(1)(0, 0, 0, 0) * km << "  "
	<< d_n_auto(2)(0, 0, 0, 0) * km << endl ;

    ost << endl <<
      "Central (d/dX,d/dY,d/dZ)(beta_auto) [km^{-1}] : "
	<< d_confpsi_auto(0)(0, 0, 0, 0) * km << "  "
	<< d_confpsi_auto(1)(0, 0, 0, 0) * km << "  "
	<< d_confpsi_auto(2)(0, 0, 0, 0) * km << endl ;

    return ost ;
}
