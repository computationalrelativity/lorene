/*
 *  Method of class Black_hole to compute the extrinsic curvature tensor
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

char blackhole_extr_curv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/05/15 19:27:14  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:19:32  k_taniguchi
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
#include "blackhole.h"
#include "utilitaires.h"
#include "unites.h"

void Black_hole::extr_curv_bh() {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (kerrschild) {

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


	// Computation of \tilde{A}^{ij}
	// -----------------------------

	Scalar divshift(mp) ;
	divshift = shift_rs(1).dsdx() + shift_rs(2).dsdy()
	  + shift_rs(3).dsdz() ;
	divshift.std_spectral_base() ;

	Sym_tensor flat_taij(mp, CON, mp.get_bvect_cart()) ;
	flat_taij.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_taij.set(i,j) = shift_rs(j).deriv(i)
		  + shift_rs(i).deriv(j)
		  - 2. * divshift * flat.con()(i,j) / 3. ;
	    }
	}

	flat_taij.std_spectral_base() ;

	Scalar lapse_bh(mp) ;
	lapse_bh = 1. / sqrt(1. + 2. * mass / rr) ;
	lapse_bh.std_spectral_base() ;
	lapse_bh.annule_domain(0) ;
	lapse_bh.raccord(1) ;

	Sym_tensor curv_taij(mp, CON, mp.get_bvect_cart()) ;
	curv_taij.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        curv_taij.set(i,j) = -2. * lapse_bh * lapse_bh * mass
		  * (ll(i)*(shift_rs(j).dsdr()) + ll(j)*(shift_rs(i).dsdr())
		     - 2. * ll(i) * ll(j) * divshift / 3.) / rr ;
	    }
	}

	curv_taij.std_spectral_base() ;

	Sym_tensor resi_taij(mp, CON, mp.get_bvect_cart()) ;
	resi_taij.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        resi_taij.set(i,j) = 2. * lapse_bh * lapse_bh * mass
		  * ( ll(i) * shift_rs(j) + ll(j) * shift_rs(i)
		      + ( flat.con()(i,j)
			  - lapse_bh*lapse_bh*(9.+14.*mass/rr)*ll(i)*ll(j) )
		      * ( ll(1)*shift_rs(1)+ll(2)*shift_rs(2)
			  +ll(3)*shift_rs(3) )/ 3. )
		  / rr / rr ;
	    }
	}

	resi_taij.std_spectral_base() ;
	resi_taij.inc_dzpuis(2) ;

	taij_rs = 0.5 * pow(confo, 7.)
	  * (flat_taij + curv_taij + resi_taij) / lapconf ;

	taij_rs.std_spectral_base() ;
	taij_rs.annule_domain(0) ;

	Sym_tensor taij_bh(mp, CON, mp.get_bvect_cart()) ;
	taij_bh.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_bh.set(i,j) = 2.*pow(lapse_bh,6.)*mass*(2.+3.*mass/rr)
		  *( (1.+2.*mass/rr) * flat.con()(i,j)
		     - (3.+2.*mass/rr) * ll(i) * ll(j) )
		  *pow(confo, 7.)/lapconf/3./rr/rr ;
	    }
	}

	taij_bh.std_spectral_base() ;
	taij_bh.inc_dzpuis(2) ;
	taij_bh.annule_domain(0) ;

	taij = taij_rs + taij_bh ;
	taij.std_spectral_base() ;
	taij.annule_domain(0) ;

	/*
	Sym_tensor taij_ks_con(mp, CON, mp.get_bvect_cart()) ;
	taij_ks_con.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_ks_con.set(i,j) = 2.*pow(lap_bh2,2.5)*mass
		  * (2.+3.*mass/rr)
		  * ( (1.+2.*mass/rr)*flat.con()(i,j)
		      - (3.+2.*mass/rr)*ll(i)*ll(j) ) / 3. / rr / rr ;
	    }
	}
	taij_ks_con.std_spectral_base() ;
	taij_ks_con.annule_domain(0) ;
	taij_ks_con.inc_dzpuis(2) ;

	cout << "taij(1,1) - taij_ks_con(1,1) : " << endl ;
	cout << taij(1,1) - taij_ks_con(1,1) << endl ;
	arrete() ;

	cout << "taij_ks_con(1,1) : " << endl ;
	cout << taij_ks_con(1,1) << endl ;
	arrete() ;

	cout << "taij(1,1) : " << endl ;
	cout << taij(1,1) << endl ;
	arrete() ;
	*/

	// Computation of \tilde{A}^{ij} \tilde{A}_{ij}
	// --------------------------------------------

	Sym_tensor flat_dshift(mp, COV, mp.get_bvect_cart()) ;
	flat_dshift.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_dshift.set(i,j) = flat.cov()(j,1)*(shift_rs(1).deriv(i))
		  + flat.cov()(j,2)*(shift_rs(2).deriv(i))
		  + flat.cov()(j,3)*(shift_rs(3).deriv(i))
		  + flat.cov()(i,1)*(shift_rs(1).deriv(j))
		  + flat.cov()(i,2)*(shift_rs(2).deriv(j))
		  + flat.cov()(i,3)*(shift_rs(3).deriv(j))
		  - 2. * divshift * flat.cov()(i,j) / 3. ;
	    }
	}

	flat_dshift.std_spectral_base() ;

	Sym_tensor curv_dshift(mp, COV, mp.get_bvect_cart()) ;
	curv_dshift.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        curv_dshift.set(i,j) = 2. * mass
		  *( ll(j) *( ll(1)*(shift_rs(1).deriv(i))
			      + ll(2)*(shift_rs(2).deriv(i))
			      + ll(3)*(shift_rs(3).deriv(i)) )
		     + ll(i) *( ll(1)*(shift_rs(1).deriv(j))
				+ ll(2)*(shift_rs(2).deriv(j))
				+ ll(3)*(shift_rs(3).deriv(j)) )
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
			    *shift_rs(1)
			    + (flat.cov()(i,2)+2.*mass*ll(i)*ll(2)/rr)
			    *shift_rs(2)
			    + (flat.cov()(i,3)+2.*mass*ll(i)*ll(3)/rr)
			    *shift_rs(3)
			    )
		    + ll(i)*( (flat.cov()(j,1)+2.*mass*ll(j)*ll(1)/rr)
			      *shift_rs(1)
			      + (flat.cov()(j,2)+2.*mass*ll(j)*ll(2)/rr)
			      *shift_rs(2)
			      + (flat.cov()(j,3)+2.*mass*ll(j)*ll(3)/rr)
			      *shift_rs(3) )
		    ) / rr / rr ;
	    }
	}
	tmp1.std_spectral_base() ;
	tmp1.inc_dzpuis(2) ;

	Sym_tensor tmp2(mp, COV, mp.get_bvect_cart()) ;
	tmp2.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        tmp2.set(i,j) = 2. * mass * lapse_bh * lapse_bh
		  * (ll(1)*shift_rs(1)+ll(2)*shift_rs(2)+ll(3)*shift_rs(3))
		  * (flat.cov()(i,j)
		     - (9.+28.*mass/rr+24.*mass*mass/rr/rr)*ll(i)*ll(j))
		  / 3. / rr / rr ;
	    }
	}
	tmp2.std_spectral_base() ;
	tmp2.inc_dzpuis(2) ;

	Sym_tensor taij_rs_down(mp, COV, mp.get_bvect_cart()) ;
	taij_rs_down.set_etat_qcq() ;

	taij_rs_down = 0.5 * pow(confo, 7.)
	  * (flat_dshift + curv_dshift + tmp1 + tmp2) / lapconf ;

	taij_rs_down.std_spectral_base() ;
	taij_rs_down.annule_domain(0) ;

	Sym_tensor taij_bh_down(mp, COV, mp.get_bvect_cart()) ;
	taij_bh_down.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	      taij_bh_down.set(i,j) = 2.*pow(lapse_bh,4.)*mass*(2.+3.*mass/rr)
		*pow(confo,7.)*(flat.cov()(i,j)-(3.+4.*mass/rr)*ll(i)*ll(j))
		/lapconf/3./rr/rr ;
	    }
	}

	taij_bh_down.std_spectral_base() ;
	taij_bh_down.inc_dzpuis(2) ;
	taij_bh_down.annule_domain(0) ;

	/*
	Sym_tensor taij_ks_cov(mp, COV, mp.get_bvect_cart()) ;
	taij_ks_cov.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	      taij_ks_cov.set(i,j) = 2.*pow(lap_bh2,1.5)*mass * (2.+3.*mass/rr)
		* ( flat.cov()(i,j)
		    - (3.+4.*mass/rr)*ll(i)*ll(j) ) / 3. / rr / rr ;
	    }
	}
	taij_ks_cov.std_spectral_base() ;
	taij_ks_cov.annule_domain(0) ;
	taij_ks_cov.inc_dzpuis(2) ;

	cout << "taij_down(1,1) - taij_ks_cov(1,1) : " << endl ;
	cout << taij_down(1,1) - taij_ks_cov(1,1) << endl ;
	arrete() ;

	cout << "taij_ks_cov(1,1) : " << endl ;
	cout << taij_ks_cov(1,1) << endl ;
	arrete() ;

	cout << "taij_down(1,1) : " << endl ;
	cout << taij_down(1,1) << endl ;
	arrete() ;
	*/

	Scalar taij_quad_rsrs(mp) ;
	taij_quad_rsrs = 0. ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad_rsrs += taij_rs_down(i,j) * taij_rs(i,j) ;
	    }
	}
	taij_quad_rsrs.std_spectral_base() ;

	Scalar taij_quad_rsbh1(mp) ;
	taij_quad_rsbh1 = 0. ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad_rsbh1 += taij_rs_down(i,j) * taij_bh(i,j) ;
	    }
	}
	taij_quad_rsbh1.std_spectral_base() ;

	Scalar taij_quad_rsbh2(mp) ;
	taij_quad_rsbh2 = 0. ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad_rsbh2 += taij_bh_down(i,j) * taij_rs(i,j) ;
	    }
	}
	taij_quad_rsbh2.std_spectral_base() ;

	taij_quad_rs = taij_quad_rsrs + taij_quad_rsbh1 + taij_quad_rsbh2 ;
	taij_quad_rs.std_spectral_base() ;

	Scalar taij_quad_bh(mp) ;
	taij_quad_bh = 8.*pow(lapse_bh,10.)*mass*mass*(2.+3.*mass/rr)
	  *(2.+3.*mass/rr)*pow(confo,12.)/3./pow(rr,4.)/lapconf/lapconf ;
	taij_quad_bh.std_spectral_base() ;
	taij_quad_bh.inc_dzpuis(4) ;

	taij_quad = taij_quad_rs + taij_quad_bh ;

	taij_quad.std_spectral_base() ;

	/*
	Scalar taij_quad_ks(mp) ;
	taij_quad_ks = 8. * pow(lap_bh2,3.) * mass * mass * (2.+3.*mass/rr)
	  * (2.+3.*mass/rr) / 3. / pow(rr, 4.) ;
	taij_quad_ks.std_spectral_base() ;
	taij_quad_ks.annule_domain(0) ;
	taij_quad_ks.inc_dzpuis(4) ;

	cout << "taij_quad - taij_quad_ks : " << endl ;
	cout << taij_quad - taij_quad_ks << endl ;
	arrete() ;
	*/
    }
    else {  // Isotropic coordinates with the maximal slicing

        // Computation of \tilde{A}^{ij}
	// -----------------------------

	Scalar divshift(mp) ;
	divshift = shift(1).dsdx() + shift(2).dsdy() + shift(3).dsdz() ;
	divshift.std_spectral_base() ;

	Sym_tensor flat_taij(mp, CON, mp.get_bvect_cart()) ;
	flat_taij.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_taij.set(i,j) = shift(j).deriv(i) + shift(i).deriv(j)
		  - 2. * divshift * flat.con()(i,j) / 3. ;
	    }
	}

	flat_taij.std_spectral_base() ;

	taij = 0.5 * pow(confo, 7.) * flat_taij / lapconf ;

	taij.std_spectral_base() ;
	taij.annule_domain(0) ;


	// Computation of \tilde{A}^{ij} \tilde{A}_{ij}
	// --------------------------------------------

	Sym_tensor flat_dshift(mp, COV, mp.get_bvect_cart()) ;
	flat_dshift.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        flat_dshift.set(i,j) = flat.cov()(j,1)*(shift(1).deriv(i))
		  + flat.cov()(j,2)*(shift(2).deriv(i))
		  + flat.cov()(j,3)*(shift(3).deriv(i))
		  + flat.cov()(i,1)*(shift(1).deriv(j))
		  + flat.cov()(i,2)*(shift(2).deriv(j))
		  + flat.cov()(i,3)*(shift(3).deriv(j))
		  - 2. * divshift * flat.cov()(i,j) / 3. ;
	    }
	}

	Sym_tensor taij_down(mp, COV, mp.get_bvect_cart()) ;
	taij_down.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	      taij_down.set(i,j) = 0.5 * pow(confo, 7.) * flat_dshift(i,j)
		  / lapconf ;
	    }
	}

	taij_down.std_spectral_base() ;
	taij_down.annule_domain(0) ;

	taij_quad = 0. ;

	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        taij_quad += taij_down(i,j) * taij(i,j) ;
	    }
	}
	taij_quad.std_spectral_base() ;

    }

}
