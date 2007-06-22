/*
 *  Method of class Hole_bhns to compute the extrinsic curvature tensor
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

char hole_bhns_extr_curv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:24:56  k_taniguchi
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
#include "utilitaires.h"
#include "unites.h"
//#include "graphique.h"

void Hole_bhns::extr_curv_bhns(double omega_orb, double x_rot, double y_rot) {

    //----------------------------------
    // Total extrinsic curvature tensor
    //----------------------------------

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    // Coordinates
    // -----------

    double mass = ggrav * mass_bh ;

    Scalar rr(mp) ;
    rr = mp.r ;
    rr.std_spectral_base() ;
    Scalar st(mp) ;
    st = mp.sint ;
    st.std_spectral_base() ;
    Scalar ct(mp) ;
    ct = mp.cost ;
    ct.std_spectral_base() ;
    Scalar sp(mp) ;
    sp = mp.sinp ;
    sp.std_spectral_base() ;
    Scalar cp(mp) ;
    cp = mp.cosp ;
    cp.std_spectral_base() ;

    Vector ll(mp, CON, mp.get_bvect_cart()) ;
    ll.set_etat_qcq() ;
    ll.set(1) = st * cp ;
    ll.set(2) = st * sp ;
    ll.set(3) = ct ;
    ll.std_spectral_base() ;

    Scalar divshift(mp) ;
    divshift = shift_auto_rs(1).deriv(1) + shift_auto_rs(2).deriv(2)
      + shift_auto_rs(3).deriv(3) + d_shift_comp(1,1)
      + d_shift_comp(2,2) + d_shift_comp(3,3) ;
    divshift.std_spectral_base() ;

    if (kerrschild) {

	Scalar orb_rot_x(mp) ;
	orb_rot_x = omega_orb * (mp.get_ori_x() - x_rot) ;
	orb_rot_x.std_spectral_base() ;

	Scalar orb_rot_y(mp) ;
	orb_rot_y = omega_orb * (mp.get_ori_y() - y_rot) ;
	orb_rot_y.std_spectral_base() ;

	Vector uv(mp, CON, mp.get_bvect_cart()) ; // unit vector
	uv.set_etat_qcq() ;
	uv.set(1) = - orb_rot_y ;
	uv.set(2) = orb_rot_x ;
	uv.set(3) = 0. ;
	uv.std_spectral_base() ;

	// Computation of \tilde{A}^{ij}
	// -----------------------------

	Sym_tensor flat_taij(mp, CON, mp.get_bvect_cart()) ;
	flat_taij.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_taij.set(i,j) = shift_auto_rs(j).deriv(i)
		  + shift_auto_rs(i).deriv(j) + d_shift_comp(i,j)
		  + d_shift_comp(j,i)
		  - 2. * divshift * flat.con()(i,j) / 3. ;
	    }
	}

	flat_taij.std_spectral_base() ;

	Sym_tensor curv_taij(mp, CON, mp.get_bvect_cart()) ;
	curv_taij.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        curv_taij.set(i,j) = -2. * lapse_auto_bh * lapse_auto_bh * mass
		  * (ll(i) * (shift_auto_rs(j).dsdr()
			      + ll(1)*d_shift_comp(1,j)
			      + ll(2)*d_shift_comp(2,j)
			      + ll(3)*d_shift_comp(3,j))
		     + ll(j) * (shift_auto_rs(i).dsdr()
				+ ll(1)*d_shift_comp(1,i)
				+ ll(2)*d_shift_comp(2,i)
				+ ll(3)*d_shift_comp(3,i))
		     - 2. * ll(i) * ll(j) * divshift / 3.) / rr ;
	    }
	}

	curv_taij.std_spectral_base() ;

	Sym_tensor resi_taij(mp, CON, mp.get_bvect_cart()) ;
	resi_taij.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        resi_taij.set(i,j) = 2. * lapse_auto_bh * lapse_auto_bh * mass
		  * ( ll(i) * (shift_auto_rs(j) + shift_comp(j))
		      + ll(j) * (shift_auto_rs(i) + shift_comp(i))
		      + ( flat.con()(i,j)
			  - lapse_auto_bh*lapse_auto_bh
			  *(9.+14.*mass/rr)*ll(i)*ll(j) )
		      * ( ll(1) * (shift_auto_rs(1) + shift_comp(1))
			  + ll(2) * (shift_auto_rs(2) + shift_comp(2))
			  + ll(3) * (shift_auto_rs(3) + shift_comp(3)) ) / 3. )
		  / rr / rr ;
	    }
	}

	resi_taij.std_spectral_base() ;
	resi_taij.inc_dzpuis(2) ;

	taij_tot_rs = 0.5 * pow(confo_tot, 6.)
	  * (flat_taij + curv_taij + resi_taij) / lapse_tot ;

	taij_tot_rs.std_spectral_base() ;
	taij_tot_rs.annule_domain(0) ;

	taij_tot_rot.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_tot_rot.set(i,j) = pow(confo_tot,6.)
		  * lapse_auto_bh * lapse_auto_bh * mass
		  * ( ll(i) * uv(j) + ll(j) * uv(i)
		      + ( flat.con()(i,j)
			  - lapse_auto_bh*lapse_auto_bh
			  *(9.+14.*mass/rr)*ll(i)*ll(j) )
		      * (ll(2)*orb_rot_x - ll(1)*orb_rot_y) / 3. )
		  / lapse_tot / rr / rr ;
	    }
	}

	taij_tot_rot.std_spectral_base() ;
	taij_tot_rot.annule_domain(0) ;
	taij_tot_rot.inc_dzpuis(2) ;

	taij_tot_bh.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_tot_bh.set(i,j) = 2. * pow(confo_tot,6.)
		  * pow(lapse_auto_bh,6.) * mass * (2.+3.*mass/rr)
		  * ( (1.+2.*mass/rr)*flat.con()(i,j)
		      - (3.+2.*mass/rr) * ll(i) * ll(j) )
		  / 3. / lapse_tot / rr / rr ;
	    }
	}

	taij_tot_bh.std_spectral_base() ;
	taij_tot_bh.annule_domain(0) ;
	taij_tot_bh.inc_dzpuis(2) ;

	taij_tot = taij_tot_rs + taij_tot_rot + taij_tot_bh ;

	taij_tot.std_spectral_base() ;
	taij_tot.annule_domain(0) ;

	// Computation of \tilde{A}^{ij} \tilde{A}_{ij}
	// --------------------------------------------

	Sym_tensor flat_dshift(mp, COV, mp.get_bvect_cart()) ;
	flat_dshift.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_dshift.set(i,j) =
		  flat.cov()(j,1)*(shift_auto_rs(1).deriv(i)
				   + d_shift_comp(i,1))
		  + flat.cov()(j,2)*(shift_auto_rs(2).deriv(i)
				     + d_shift_comp(i,2))
		  + flat.cov()(j,3)*(shift_auto_rs(3).deriv(i)
				     + d_shift_comp(i,3))
		  + flat.cov()(i,1)*(shift_auto_rs(1).deriv(j)
				     + d_shift_comp(j,1))
		  + flat.cov()(i,2)*(shift_auto_rs(2).deriv(j)
				     + d_shift_comp(j,2))
		  + flat.cov()(i,3)*(shift_auto_rs(3).deriv(j)
				     + d_shift_comp(j,3))
		  - 2. * divshift * flat.cov()(i,j) / 3. ;
	    }
	}

	flat_dshift.std_spectral_base() ;

	Sym_tensor curv_dshift(mp, COV, mp.get_bvect_cart()) ;
	curv_dshift.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        curv_dshift.set(i,j) = 2. * mass
		  *( ll(j) * ( ll(1)*(shift_auto_rs(1).deriv(i)
				      + d_shift_comp(i,1))
			       + ll(2)*(shift_auto_rs(2).deriv(i)
					+ d_shift_comp(i,2))
			       + ll(3)*(shift_auto_rs(3).deriv(i)
					+ d_shift_comp(i,3)))
		     + ll(i) * ( ll(1)*(shift_auto_rs(1).deriv(j)
					+ d_shift_comp(j,1))
				 + ll(2)*(shift_auto_rs(2).deriv(j)
					  + d_shift_comp(j,2))
				 + ll(3)*(shift_auto_rs(3).deriv(j)
					  + d_shift_comp(j,3)))
		     - 2. * divshift * ll(i) * ll(j) / 3. ) / rr ;
	    }
	}

	curv_dshift.std_spectral_base() ;

	Sym_tensor tmp1(mp, COV, mp.get_bvect_cart()) ;
	tmp1.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        tmp1.set(i,j) = 2. * mass
		  *(ll(j)*( (flat.cov()(i,1)+2.*mass*ll(i)*ll(1)/rr)
			    * (shift_auto_rs(1) + shift_comp(1))
			    + (flat.cov()(i,2)+2.*mass*ll(i)*ll(2)/rr)
			    * (shift_auto_rs(2) + shift_comp(2))
			    + (flat.cov()(i,3)+2.*mass*ll(i)*ll(3)/rr)
			    * (shift_auto_rs(3) + shift_comp(3))
			    )
		    + ll(i)*( (flat.cov()(j,1)+2.*mass*ll(j)*ll(1)/rr)
			      * (shift_auto_rs(1) + shift_comp(1))
			      + (flat.cov()(j,2)+2.*mass*ll(j)*ll(2)/rr)
			      * (shift_auto_rs(2) + shift_comp(2))
			      + (flat.cov()(j,3)+2.*mass*ll(j)*ll(3)/rr)
			      * (shift_auto_rs(3) + shift_comp(3)) )
		    ) / rr / rr ;
	    }
	}
	tmp1.std_spectral_base() ;
	tmp1.inc_dzpuis(2) ;

	Sym_tensor tmp2(mp, COV, mp.get_bvect_cart()) ;
	tmp2.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        tmp2.set(i,j) = 2. * mass * lapse_auto_bh * lapse_auto_bh
		  * ( ll(1) * (shift_auto_rs(1) + shift_comp(1))
		      + ll(2) * (shift_auto_rs(2) + shift_comp(2))
		      + ll(3) * (shift_auto_rs(3) + shift_comp(3)) )
		  * (flat.cov()(i,j)
		     - (9.+28.*mass/rr+24.*mass*mass/rr/rr)*ll(i)*ll(j))
		  / 3. / rr / rr ;
	    }
	}
	tmp2.std_spectral_base() ;
	tmp2.inc_dzpuis(2) ;

	Sym_tensor taij_down_rs(mp, COV, mp.get_bvect_cart()) ;
	taij_down_rs.set_etat_qcq() ;

	taij_down_rs = 0.5 * pow(confo_tot,6.)
	  * (flat_dshift + curv_dshift + tmp1 + tmp2) / lapse_tot ;

	taij_down_rs.std_spectral_base() ;
	taij_down_rs.annule_domain(0) ;

	Sym_tensor taij_down_rot(mp, COV, mp.get_bvect_cart()) ;
	taij_down_rot.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	      taij_down_rot.set(i,j) = mass * pow(confo_tot,6.)
		* ( ll(j)*uv(i) + ll(i)*uv(j)
		    + lapse_auto_bh * lapse_auto_bh
		    * (ll(2)*orb_rot_x - ll(1)*orb_rot_y)
		    * ( flat.cov()(i,j) - (9.+16.*mass/rr)*ll(i)*ll(j) ) / 3.
		    ) / lapse_tot / rr / rr ;
	    }
	}
	taij_down_rot.std_spectral_base() ;
	taij_down_rot.annule_domain(0) ;
	taij_down_rot.inc_dzpuis(2) ;

	Sym_tensor taij_down_bh(mp, COV, mp.get_bvect_cart()) ;
	taij_down_bh.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	      taij_down_bh.set(i,j) = 2. * mass * pow(confo_tot,6.)
		* pow(lapse_auto_bh,4.) * (2.+3.*mass/rr)
		* (flat.cov()(i,j) - (3.+4.*mass/rr) * ll(i) * ll(j))
		/ 3. / lapse_auto / rr / rr ;
	    }
	}
	taij_down_bh.std_spectral_base() ;
	taij_down_bh.annule_domain(0) ;
	taij_down_bh.inc_dzpuis(2) ;

	Scalar taij_quad_rstot(mp) ;
	taij_quad_rstot = 0. ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad_rstot += taij_down_rs(i,j) * taij_tot(i,j) ;
	    }
	}
	taij_quad_rstot.std_spectral_base() ;

	Scalar taij_quad_rsrotbh(mp) ;
	taij_quad_rsrotbh = 0. ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad_rsrotbh += taij_tot_rs(i,j)
		  * (taij_down_rot(i,j) + taij_down_bh(i,j)) ;
	    }
	}
	taij_quad_rsrotbh.std_spectral_base() ;

	taij_quad_tot_rs = taij_quad_rstot + taij_quad_rsrotbh ;
	taij_quad_tot_rs.std_spectral_base() ;

	taij_quad_tot_rot = 8.*mass*mass*pow(confo_tot,12.)
	  * pow(lapse_auto_bh,6.) * (2.+3.*mass/rr)
	  * (ll(2)*orb_rot_x - ll(1)*orb_rot_y)
	  / 3. / lapse_tot / lapse_tot / pow(rr,4.)
	  + 2.*mass*mass*pow(confo_tot,12.)*pow(lapse_auto_bh,4.)
	  * (3.*(1.+2.*mass/rr)*(orb_rot_x*orb_rot_x+orb_rot_y*orb_rot_y)
	     -2.*(1.+3.*mass/rr)*(ll(2)*orb_rot_x-ll(1)*orb_rot_y)
	     *(ll(2)*orb_rot_x-ll(1)*orb_rot_y))
	  / 3. / lapse_tot / lapse_tot / pow(rr,4.) ;

	taij_quad_tot_rot.std_spectral_base() ;
	taij_quad_tot_rot.inc_dzpuis(4) ;

	taij_quad_tot_bh = 8.*mass*mass*pow(confo_tot,12.)
	  * pow(lapse_auto_bh,8.) * (2.+3.*mass/rr) * (2.+3.*mass/rr)
	  / 3. / lapse_tot / lapse_tot / pow(rr,4.) ;

	taij_quad_tot_bh.std_spectral_base() ;
	taij_quad_tot_bh.inc_dzpuis(4) ;

	taij_quad_tot = taij_quad_tot_rs + taij_quad_tot_rot
	  + taij_quad_tot_bh ;
	taij_quad_tot.std_spectral_base() ;

    }
    else { // Isotropic coordinates with the maximal slicing

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------
        double cc ;

	if (bc_lapse_nd) {  // Neumann boundary condition
	    if (bc_lapse_fs) {  // First condition
	        // d\alpha/dr = 0
	        // --------------
	        cc = 2. ;
	    }
	    else {  // Second condition
	        // d\alpha/dr = \alpha/(2 rah)
	        // ---------------------------
	        cc = 0.5 * (sqrt(17.) - 1.) ;
	    }
	}
	else {  // Dirichlet boundary condition
	    if (bc_lapse_fs) {  // First condition
	        // \alpha = 1/2
	        // ------------
	        cc = 2. ;
	    }
	    else {  // Second condition
	        // \alpha = 1/sqrt(2.)
	        // -------------------
	        cc = 2. * sqrt(2.) ;
	    }
	}

	Scalar r_are(mp) ;
	r_are = r_coord(bc_lapse_nd, bc_lapse_fs) ;
	r_are.std_spectral_base() ;

        // Computation of \tilde{A}^{ij}
        // -----------------------------

	Sym_tensor flat_taij(mp, CON, mp.get_bvect_cart()) ;
	flat_taij.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_taij.set(i,j) = shift_auto_rs(j).deriv(i)
		  + shift_auto_rs(i).deriv(j) + d_shift_comp(i,j)
		  + d_shift_comp(j,i)
		  - 2. * divshift * flat.con()(i,j) / 3. ;
	    }
	}

	flat_taij.std_spectral_base() ;

	taij_tot_rs = 0.5 * pow(confo_tot, 6.) * flat_taij / lapse_tot ;

	taij_tot_rs.std_spectral_base() ;
	taij_tot_rs.annule_domain(0) ;

	taij_tot_bh.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_tot_bh.set(i,j) = pow(confo_tot,6.)*mass*mass*cc
		  * sqrt(1. - 2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
		  * (flat.con()(i,j) - 3.*ll(i)*ll(j)) / lapse_tot
		  / pow(r_are*rr,3.) ;
	    }
	}

	taij_tot_bh.std_spectral_base() ;
	taij_tot_bh.annule_domain(0) ;
	taij_tot_bh.inc_dzpuis(2) ;

	taij_tot = taij_tot_rs + taij_tot_bh ;

	taij_tot.std_spectral_base() ;
	taij_tot.annule_domain(0) ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
		taij_tot_rot.set(i,j) = 0. ;
	    }
	}
	taij_tot_rot.std_spectral_base() ;

	// Computation of \tilde{A}_{BH}^{ij}
	// ----------------------------------

	Scalar divshift_auto(mp) ;
	divshift_auto = shift_auto_rs(1).deriv(1)
	  + shift_auto_rs(2).deriv(2) + shift_auto_rs(3).deriv(3) ;
	divshift_auto.std_spectral_base() ;

	Sym_tensor flat_taij_auto_rs(mp, CON, mp.get_bvect_cart()) ;
	flat_taij_auto_rs.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_taij_auto_rs.set(i,j) = shift_auto_rs(j).deriv(i)
		  + shift_auto_rs(i).deriv(j)
		  - 2. * divshift_auto * flat.con()(i,j) / 3. ;
	    }
	}

	flat_taij_auto_rs.std_spectral_base() ;

	taij_auto_rs = 0.5 * pow(confo_tot, 6.) * flat_taij_auto_rs
	  / lapse_tot ;

	taij_auto_rs.std_spectral_base() ;
	taij_auto_rs.annule_domain(0) ;

	taij_auto = taij_auto_rs + taij_tot_bh ;

	taij_auto.std_spectral_base() ;
	taij_auto.annule_domain(0) ;


	// Computation of \tilde{A}_{NS}^{ij}
	// ----------------------------------

	Scalar divshift_comp(mp) ;
	divshift_comp = d_shift_comp(1,1) + d_shift_comp(2,2)
	  + d_shift_comp(3,3) ;
	divshift_comp.std_spectral_base() ;

	Sym_tensor flat_taij_comp(mp, CON, mp.get_bvect_cart()) ;
	flat_taij_comp.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_taij_comp.set(i,j) = d_shift_comp(i,j)
		  + d_shift_comp(j,i)
		  - 2. * divshift_comp * flat.con()(i,j) / 3. ;
	    }
	}

	flat_taij_comp.std_spectral_base() ;

	taij_comp = 0.5 * pow(confo_comp+0.5, 6.) * flat_taij_comp
	  / (lapse_comp+0.5) ;

	taij_comp.std_spectral_base() ;
	taij_comp.annule_domain(0) ;


	// Computation of \tilde{A}^{ij} \tilde{A}_{ij}
	// --------------------------------------------

	Sym_tensor flat_dshift(mp, COV, mp.get_bvect_cart()) ;
	flat_dshift.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_dshift.set(i,j) =
		  flat.cov()(j,1)*(shift_auto_rs(1).deriv(i)
				   + d_shift_comp(i,1))
		  + flat.cov()(j,2)*(shift_auto_rs(2).deriv(i)
				     + d_shift_comp(i,2))
		  + flat.cov()(j,3)*(shift_auto_rs(3).deriv(i)
				     + d_shift_comp(i,3))
		  + flat.cov()(i,1)*(shift_auto_rs(1).deriv(j)
				     + d_shift_comp(j,1))
		  + flat.cov()(i,2)*(shift_auto_rs(2).deriv(j)
				     + d_shift_comp(j,2))
		  + flat.cov()(i,3)*(shift_auto_rs(3).deriv(j)
				     + d_shift_comp(j,3))
		  - 2. * divshift * flat.cov()(i,j) / 3. ;
	    }
	}

	Sym_tensor taij_down_rs(mp, COV, mp.get_bvect_cart()) ;
	taij_down_rs.set_etat_qcq() ;

	taij_down_rs = 0.5 * pow(confo_tot, 6.) * flat_dshift / lapse_tot ;

	taij_down_rs.std_spectral_base() ;
	taij_down_rs.annule_domain(0) ;

	Sym_tensor taij_down_bh(mp, COV, mp.get_bvect_cart()) ;
	taij_down_bh.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_down_bh.set(i,j) = pow(confo_tot,6.)*mass*mass*cc
		  * sqrt(1. - 2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
		  * (flat.cov()(i,j) - 3.*ll(i)*ll(j)) / lapse_tot
		  / pow(r_are*rr,3.) ;
	    }
	}
	taij_down_bh.std_spectral_base() ;
	taij_down_bh.annule_domain(0) ;
	taij_down_bh.inc_dzpuis(2) ;

	Scalar taij_quad_rstot(mp) ;
	taij_quad_rstot = 0. ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad_rstot += taij_down_rs(i,j) * taij_tot(i,j) ;
	    }
	}
	taij_quad_rstot.std_spectral_base() ;

	Scalar taij_quad_rsbh(mp) ;
	taij_quad_rsbh = 0. ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad_rsbh += taij_tot_rs(i,j) * taij_down_bh(i,j) ;
	    }
	}
	taij_quad_rsbh.std_spectral_base() ;

	taij_quad_tot_rs = taij_quad_rstot + taij_quad_rsbh ;
	taij_quad_tot_rs.std_spectral_base() ;

	taij_quad_tot_rot = 0. ;
	taij_quad_tot_rot.std_spectral_base() ;

	taij_quad_tot_bh = 6.*pow(confo_tot,12.)*pow(mass*mass*cc,2.)
	  * (1. - 2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
	  / lapse_tot / lapse_tot / pow(r_are*rr, 6.) ;
	taij_quad_tot_bh.std_spectral_base() ;
	taij_quad_tot_bh.annule_domain(0) ;
	taij_quad_tot_bh.inc_dzpuis(4) ;

	taij_quad_tot = taij_quad_tot_rs + taij_quad_tot_bh ;
	taij_quad_tot.std_spectral_base() ;


	// -------------------------
	Scalar taij_quad_auto1(mp) ;
	taij_quad_auto1 = 0. ;
	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad_auto1 += taij_auto_rs(i,j)
		  * (taij_down_rs(i,j)
		     + pow(confo_tot/(confo_comp+0.5),6.)*(lapse_comp+0.5)
		     * taij_comp(i,j) / lapse_tot) ;
	    }
	}
	taij_quad_auto1.std_spectral_base() ;

	Scalar taij_quad_auto2(mp) ;
	taij_quad_auto2 = 0. ;
	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad_auto2 += taij_tot_bh(i,j) * taij_down_rs(i,j) ;
	    }
	}
	taij_quad_auto2.std_spectral_base() ;

	taij_quad_auto = taij_quad_auto1 + 2.*taij_quad_auto2 ;
	taij_quad_auto.std_spectral_base() ;

	// Computation of \tilde{A}_{NS}^{ij} \tilde{A}^{NS}_{ij}
	// ------------------------------------------------------

	taij_quad_comp = 0. ;
	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	      taij_quad_comp += taij_comp(i,j) * taij_comp(i,j) ;
	    }
	}
	taij_quad_comp.std_spectral_base() ;

    }

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;

}
