/*
 *  Methods Et_bin_bhns_extr::update_metric_der_comp_extr_ks
 *  and Et_bin_bhns_extr::update_metric_der_comp_extr_cf
 *
 *    (see file et_bin_bhns_extr.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004-2005 Keisuke Taniguchi
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

char et_bin_bhns_extr_um_der_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/02/28 23:15:52  k_taniguchi
 * Modification to include the case of the conformally flat background metric
 *
 * Revision 1.1  2004/11/30 20:51:09  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include<math.h>

// Lorene headers
#include "et_bin_bhns_extr.h"
#include "etoile.h"
#include "coord.h"
#include "unites.h"

void Et_bin_bhns_extr::update_metric_der_comp_extr(const double& mass,
						   const double& sepa) {

  using namespace Unites ;

    if (kerrschild) {

        // Computation of quantities coming from the companion (K-S BH)
        // ------------------------------------------------------------

        const Coord& xx = mp.x ;
	const Coord& yy = mp.y ;
	const Coord& zz = mp.z ;

	Tenseur r_bh(mp) ;
	r_bh.set_etat_qcq() ;
	r_bh.set() = pow( (xx+sepa)*(xx+sepa) + yy*yy + zz*zz, 0.5) ;
	r_bh.set_std_base() ;

	Tenseur xx_cov(mp, 1, COV, ref_triad) ;
	xx_cov.set_etat_qcq() ;
	xx_cov.set(0) = xx + sepa ;
	xx_cov.set(1) = yy ;
	xx_cov.set(2) = zz ;
	xx_cov.set_std_base() ;

	Tenseur xsr_cov(mp, 1, COV, ref_triad) ;
	xsr_cov = xx_cov / r_bh ;
	xsr_cov.set_std_base() ;

	Tenseur xx_con(mp, 1, CON, ref_triad) ;
	xx_con.set_etat_qcq() ;
	xx_con.set(0) = xx + sepa ;
	xx_con.set(1) = yy ;
	xx_con.set(2) = zz ;
	xx_con.set_std_base() ;

	Tenseur xsr_con(mp, 1, CON, ref_triad) ;
	xsr_con = xx_con / r_bh ;
	xsr_con.set_std_base() ;

	Tenseur msr(mp) ;
	msr = ggrav * mass / r_bh ;
	msr.set_std_base() ;

	// Computation of d_logn_comp
	// --------------------------

	Tenseur lapse_bh2(mp) ;  // lapse_bh * lapse_bh
	lapse_bh2 = 1. / (1.+2.*msr) ;
	lapse_bh2.set_std_base() ;

	d_logn_comp.set_etat_qcq() ;

	d_logn_comp.set(0) = lapse_bh2()%msr()%xsr_cov(0)/r_bh() ;
	d_logn_comp.set(1) = lapse_bh2()%msr()%xsr_cov(1)/r_bh() ;
	d_logn_comp.set(2) = lapse_bh2()%msr()%xsr_cov(2)/r_bh() ;

	d_logn_comp.set_std_base() ;
	d_logn_comp.change_triad(ref_triad) ;

	// Computation of d_beta_comp (Just the same that d_logn_comp)
	// --------------------------

	d_beta_comp = d_logn_comp ;
	d_beta_comp.change_triad(ref_triad) ;

	// Computation of tkij_comp
	// ------------------------
	/** The definition of tkij_comp:
	 *  A_{tot}^2 A_{BH}^{ij} =-0.5/nnn*(DN_{BH}+DN_{BH}-2/3gammaDN_{BH})
	 */

	// Components of shift_comp with respect to the Cartesian triad
	//  (d/dx, d/dy, d/dz) of the mapping :

	// Computation of A^2 A^{ij}
	// -------------------------
	tkij_comp.set_etat_qcq() ;

	for (int i=0; i<3; i++) {
	    for (int j=i; j<3; j++) {
	        tkij_comp.set(i, j) = - 3.*xsr_con(i) % xsr_con(j)
		  - 2*msr() % xsr_con(i) % xsr_con(j) ;
	    }
	    tkij_comp.set(i, i) += 1. + 2.*msr() ;
	}

	tkij_comp = (double(2)/double(3)) * pow(lapse_bh2, 3.) % msr
	  % (2.*tkij_comp +3.*msr % tkij_comp) / nnn / r_bh ;

	tkij_comp.set_triad(ref_triad) ;
	tkij_comp.set_std_base() ;

	if (relativistic) {
	    // Computation of akcar_comp
	    // -------------------------

	    Tenseur lapse_bh8(mp) ;
	    lapse_bh8 = 1. / pow(1.+2.*msr, 4.) ;
	    lapse_bh8.set_std_base() ;

	    akcar_comp.set_etat_qcq() ;
	    akcar_comp.set() = 0 ;

	    akcar_comp.set() = 8.*lapse_bh8()
	      * pow(2.*msr()+3.*msr()*msr(), 2.) / 3.
	      / nnn() / nnn() / r_bh() / r_bh() ;

	    akcar_comp.set_std_base() ;
	    akcar_comp = a_car % akcar_comp ;

	}

	// The derived quantities are obsolete
	// -----------------------------------

	Etoile_bin::del_deriv() ;

    }
    else {

        // Computation of quantities coming from the companion (CF Sch. BH)
        // ----------------------------------------------------------------

        const Coord& xx = mp.x ;
	const Coord& yy = mp.y ;
	const Coord& zz = mp.z ;

	Tenseur r_bh(mp) ;
	r_bh.set_etat_qcq() ;
	r_bh.set() = pow( (xx+sepa)*(xx+sepa) + yy*yy + zz*zz, 0.5) ;
	r_bh.set_std_base() ;

	Tenseur xx_cov(mp, 1, COV, ref_triad) ;
	xx_cov.set_etat_qcq() ;
	xx_cov.set(0) = xx + sepa ;
	xx_cov.set(1) = yy ;
	xx_cov.set(2) = zz ;
	xx_cov.set_std_base() ;

	Tenseur msr(mp) ;
	msr = ggrav * mass / r_bh ;
	msr.set_std_base() ;

	// Computation of d_logn_comp
	// --------------------------

	Tenseur tmp(mp) ;
	tmp = 1. / (1. - 0.25*msr*msr) ;
	tmp.set_std_base() ;

	d_logn_comp.set_etat_qcq() ;

	d_logn_comp.set(0) = tmp()%msr()%xx_cov(0)/r_bh()/r_bh() ;
	d_logn_comp.set(1) = tmp()%msr()%xx_cov(1)/r_bh()/r_bh() ;
	d_logn_comp.set(2) = tmp()%msr()%xx_cov(2)/r_bh()/r_bh() ;

	d_logn_comp.set_std_base() ;
	d_logn_comp.change_triad(ref_triad) ;

	// Computation of d_beta_comp
	// --------------------------

	d_beta_comp.set_etat_qcq() ;

	d_beta_comp.set(0) = 0.5*tmp()%msr()%msr()%xx_cov(0)/r_bh()/r_bh() ;
	d_beta_comp.set(1) = 0.5*tmp()%msr()%msr()%xx_cov(1)/r_bh()/r_bh() ;
	d_beta_comp.set(2) = 0.5*tmp()%msr()%msr()%xx_cov(2)/r_bh()/r_bh() ;

	d_beta_comp.set_std_base() ;
	d_beta_comp.change_triad(ref_triad) ;

	// Computation of tkij_comp
	// ------------------------

	// Computation of A^2 A^{ij}
	// -------------------------
	tkij_comp.set_etat_qcq() ;

	for (int i=0; i<3; i++) {
	    for (int j=i; j<3; j++) {
	        tkij_comp.set(i, j) = 0. ;
	    }
	}

	tkij_comp.set_triad(ref_triad) ;
	tkij_comp.set_std_base() ;

	if (relativistic) {
	    // Computation of akcar_comp
	    // -------------------------

	    akcar_comp.set_etat_qcq() ;
	    akcar_comp.set() = 0 ;
	    akcar_comp.set_std_base() ;

	}

	// The derived quantities are obsolete
	// -----------------------------------

	Etoile_bin::del_deriv() ;

    }

}
