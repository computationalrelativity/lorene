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
      confpsi(mp_i),
      confpsi_auto(mp_i),
      confpsi_comp(mp_i),
      ssjm1_n_auto(mp_i),
      ssjm1_confpsi(mp_i) {

    // The metric is initialized to the flat one :
    n_auto = 0.5 ;
    n_auto.set_std_base() ;
    n_comp = 0.5 ;
    n_comp.set_std_base() ;
    confpsi = 1 ;
    confpsi.set_std_base() ;
    confpsi_auto = 0.5 ;
    confpsi_auto.set_std_base() ;
    confpsi_comp = 0.5 ;
    confpsi_comp.set_std_base() ;

    ssjm1_n_auto.set_etat_qcq() ;
    ssjm1_n_auto = 0.5 ;
    ssjm1_confpsi.set_etat_qcq() ;
    ssjm1_confpsi = 0.5 ;

}

// Copy constructor
// ----------------
Et_bin_nsbh::Et_bin_nsbh(const Et_bin_nsbh& et)
    : Etoile_bin(et),
      n_auto(et.n_auto),
      n_comp(et.n_comp),
      confpsi(et.confpsi),
      confpsi_auto(et.confpsi_auto),
      confpsi_comp(et.confpsi_comp),
      ssjm1_n_auto(et.ssjm1_n_auto),
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
      confpsi(mp_i),
      confpsi_auto(mp_i),
      confpsi_comp(mp_i),
      ssjm1_n_auto(mp_i),
      ssjm1_confpsi(mp_i) {

    // Read of the saved fields :
    Tenseur ent_file(mp_i, fich) ;
    ent = ent_file ;

    Tenseur n_auto_file(mp_i, fich) ;
    n_auto = n_auto_file ;

    Tenseur confpsi_auto_file(mp_i, fich) ;
    confpsi_auto = confpsi_auto_file ;

    Cmp ssjm1_n_auto_file(mp_i, *(mp_i.get_mg()), fich) ;
    ssjm1_n_auto = ssjm1_n_auto_file ;

    Cmp ssjm1_confpsi_file(mp_i, *(mp_i.get_mg()), fich) ;
    ssjm1_confpsi = ssjm1_confpsi_file ;

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
    confpsi = et.confpsi ;
    confpsi_auto = et.confpsi_auto ;
    confpsi_comp = et.confpsi_comp ;
    ssjm1_n_auto = et.ssjm1_n_auto ;
    ssjm1_confpsi = et.ssjm1_confpsi ;

}

                            //--------------//
                            //    Outputs   //
                            //--------------//

// Save in a file
// --------------
void Et_bin_nsbh::sauve(FILE* fich) const {

    Etoile::sauve(fich) ;

    fwrite(&irrotational, sizeof(bool), 1, fich) ;

    if (irrotational) {
	gam_euler.sauve(fich) ; // required to construct d_psi from psi0
	psi0.sauve(fich) ;
    }

    w_shift.sauve(fich) ;
    khi_shift.sauve(fich) ;

    ssjm1_n_auto.sauve(fich) ;
    ssjm1_confpsi.sauve(fich) ;
    ssjm1_khi.sauve(fich) ;
    ssjm1_wshift.sauve(fich) ;

}
