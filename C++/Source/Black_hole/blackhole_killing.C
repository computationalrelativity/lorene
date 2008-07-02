/*
 *  Methods of class Black_hole to compute Killing vectors
 *
 *    (see file blackhole.h for documentation).
 *
 */

/*
 *   Copyright (c) 2007 Keisuke Taniguchi
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

char blackhole_killing_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/07/02 20:45:07  k_taniguchi
 * A bug removed.
 *
 * Revision 1.1  2008/05/15 19:33:12  k_taniguchi
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
#include "unites.h"
#include "utilitaires.h"

               //---------------------------------------------//
               //          Killing vectors on the AH          //
               //---------------------------------------------//

Vector Black_hole::killing_vect_bh(const Tbl& xi_i, const double& phi_i,
				   const double& theta_i, const int& nrk_phi,
				   const int& nrk_theta) const {

    using namespace Unites ;

    assert(xi_i.get_ndim() == 1) ;
    assert(xi_i.get_dim(0) == 3) ;

    const Mg3d* mg = mp.get_mg() ;
    int nr = mg->get_nr(1) ;
    int nt = mg->get_nt(1) ;
    int np = mg->get_np(1) ;

    // Vector which is returned to the main code
    // Spherical basis, covariant
    Vector killing(mp, COV, mp.get_bvect_spher()) ;

    if (kerrschild) {

      cout << "Not yet prepared!!!" << endl ;
      abort() ;

    }
    else {  // Isotropic coordinates

        // Solution of the Killing vector on the equator
        // ---------------------------------------------

        double dp = 2. * M_PI / double(np) ;  // Angular step

	// Killing vector on the equator
	// np+1 is for the check of xi(phi=0)= xi(phi=2pi)
	Tbl xi_t(np+1) ;
	xi_t.set_etat_qcq() ;
	Tbl xi_p(np+1) ;
	xi_p.set_etat_qcq() ;
	Tbl xi_l(np+1) ;
	xi_l.set_etat_qcq() ;

	xi_t.set(0) = xi_i(0) ;
	xi_p.set(0) = xi_i(1) ;
	xi_l.set(0) = xi_i(2) ;

	Tbl xi(3) ;
	xi.set_etat_qcq() ;

	Tbl xi_ini(3) ;
	xi_ini.set_etat_qcq() ;
	xi_ini.set(0) = xi_i(0) ;
	xi_ini.set(1) = xi_i(1) ;
	xi_ini.set(2) = xi_i(2) ;

	double pp_0 = phi_i ;  // azimuthal angle phi

	for (int i=1; i<np+1; i++) {

	    xi = runge_kutta_phi_bh(xi_ini, pp_0, nrk_phi) ;

	    xi_t.set(i) = xi(0) ;
	    xi_p.set(i) = xi(1) ;
	    xi_l.set(i) = xi(2) ;

	    // New data for the next step
	    //  -------------------------
	    pp_0 += dp ;  // New angle
	    xi_ini = xi ;

	}

	/*
	for (int i=0; i<np+1; i++) {

	  cout << "xi_t    xi_p    xi_l" << endl ;
	  cout << xi_t(i) << "  " << xi_p(i) << "  " << xi_l(i) << endl ;

	}
	arrete() ;
	*/

	// Normalization of the Killing vector
	// -----------------------------------

	// Initialization of the Killing vector to the phi direction
	Scalar xi_phi(mp) ;
	xi_phi = 0.5 ;  // If we use "1." for the initialization value,
	                // the state flag becomes "ETATUN" which does not
	                // work for "set_grid_point".

	for (int k=0; k<np; k++) {
	  xi_phi.set_grid_point(0, k, nt-1, nr-1) = xi_p(k) ;
	  xi_phi.set_grid_point(1, k, nt-1, 0) = xi_p(k) ;
	}
	xi_phi.std_spectral_base() ;
	/*
	for (int l=0; l<2; l++) {
	  for (int k=0; k<np; k++) {
	    for (int j=0; j<nt; j++) {
	      for (int i=0; i<nr; i++) {
		cout << "(l,k,j,i)=" << l << "," << k << "," << j
		     << "," << i << ": "
		     << xi_phi.val_grid_point(l,k,j,i) << endl ;
	      }
	      arrete() ;
	    }
	    arrete() ;
	  }
	  arrete() ;
	}
	*/

	// Normalization of the Killing vector
	Scalar rr(mp) ;
	rr = mp.r ;
	rr.std_spectral_base() ;

	Scalar st(mp) ;
	st = mp.sint ;
	st.std_spectral_base() ;

	Scalar source_phi(mp) ;
	source_phi = pow(confo, 2.) * rr * st / xi_phi ;
	source_phi.std_spectral_base() ;

	double rah = rad_ah() ;

	int nn = 1000 ;
	double hh = 2. * M_PI / double(nn) ;
	double integ = 0. ;

	int mm ;
	double t1, t2, t3, t4, t5 ;

	// Boole's Rule (Newton-Cotes Integral) for integration
	// ----------------------------------------------------

	assert(nn%4 == 0) ;
	mm = nn/4 ;

	for (int i=0; i<mm; i++) {

	    t1 = hh * double(4*i) ;
	    t2 = hh * double(4*i+1) ;
	    t3 = hh * double(4*i+2) ;
	    t4 = hh * double(4*i+3) ;
	    t5 = hh * double(4*i+4) ;

	    integ += (hh/45.) * (14.*source_phi.val_point(rah,M_PI/2.,t1)
				 + 64.*source_phi.val_point(rah,M_PI/2.,t2)
				 + 24.*source_phi.val_point(rah,M_PI/2.,t3)
				 + 64.*source_phi.val_point(rah,M_PI/2.,t4)
				 + 14.*source_phi.val_point(rah,M_PI/2.,t5)
				 ) ;

	}

	cout << "Black_hole:: t_f = " << integ << endl ;
	double ratio = 0.5 * integ / M_PI ;

	cout << "Black_hole:: t_f / 2M_PI = " << ratio << endl ;

	for (int k=0; k<np; k++) {
	  xi_p.set(k) = xi_phi.val_grid_point(1, k, nt-1, 0) * ratio ;
	}

	/*
	for (int k=0; k<np; k++) {
	  cout << "Normalized xi_p" << "(" << k << ") :" << xi_p(k) << endl ;
	}
	*/

        // Solution of the Killing vector to the pole angle
        // ------------------------------------------------

	double dt = 0.5 * M_PI / double(nt-1) ;  // Angular step

	// Killing vector to the polar angle
	Tbl xi_th(nt, np) ;
	xi_th.set_etat_qcq() ;
	Tbl xi_ph(nt, np) ;
	xi_ph.set_etat_qcq() ;
	Tbl xi_ll(nt, np) ;
	xi_ll.set_etat_qcq() ;

	// Values on the equator
	for (int i=0; i<np; i++) {

	    xi_th.set(nt-1, i) = xi_t(i) ;
	    xi_ph.set(nt-1, i) = xi_p(i) ;
	    xi_ll.set(nt-1, i) = xi_l(i) ;

	}

	for (int i=0; i<np; i++) {

	    // Value at theta=pi/2, phi=phi_i
	    xi_ini.set(0) = xi_t(i) ;
	    xi_ini.set(1) = xi_p(i) ;
	    xi_ini.set(2) = xi_l(i) ;

	    double pp = double(i) * dp ;
	    double tt_0 = theta_i ;  // polar angle theta

	    for (int j=1; j<nt; j++) {

	        xi = runge_kutta_theta_bh(xi_ini, tt_0, pp, nrk_theta) ;

		xi_th.set(nt-1-j, i) = xi(0) ;
		xi_ph.set(nt-1-j, i) = xi(1) ;
		xi_ll.set(nt-1-j, i) = xi(2) ;

		// New data for the nxt step
		// -------------------------
		tt_0 -= dt ;  // New angle
		xi_ini = xi ;

	    }  // End of the loop to the polar direction

	}  // End of the loop to the azimuhtal direction


	// Construction of the Killing vector
	// ----------------------------------

	// Definition of the vector is at the top of this code
	killing.set_etat_qcq() ;
	killing.set(1) = 0. ;  // r component
	killing.set(2) = 0.5 ;  // initialization of the theta component
	killing.set(3) = 0.5 ;  // initialization of the phi component

	for (int l=0; l<2; l++) {
	  for (int i=0; i<nr; i++) {
	    for (int j=0; j<nt; j++) {
	      for (int k=0; k<np; k++) {
		(killing.set(2)).set_grid_point(l, k, j, i) = xi_th(j, k) ;
		(killing.set(3)).set_grid_point(l, k, j, i) = xi_ph(j, k) ;
	      }
	    }
	  }
	}
	killing.std_spectral_base() ;

	// Check the normalization
	// -----------------------

	double check_norm = 0. ;
	source_phi = pow(confo, 2.) * rr * st / killing(3) ;
	source_phi.std_spectral_base() ;

	for (int i=0; i<mm; i++) {

	    t1 = hh * double(4*i) ;
	    t2 = hh * double(4*i+1) ;
	    t3 = hh * double(4*i+2) ;
	    t4 = hh * double(4*i+3) ;
	    t5 = hh * double(4*i+4) ;

	    check_norm += (hh/45.) *
	      ( 14.*source_phi.val_point(rah,M_PI/4.,t1)
		+ 64.*source_phi.val_point(rah,M_PI/4.,t2)
		+ 24.*source_phi.val_point(rah,M_PI/4.,t3)
		+ 64.*source_phi.val_point(rah,M_PI/4.,t4)
		+ 14.*source_phi.val_point(rah,M_PI/4.,t5) ) ;

	}

	cout << "Black_hole:: t_f for M_PI/4 = " << check_norm / M_PI
	     << " M_PI" << endl ;

    }  // End of the loop for isotropic coordinates

    return killing ;

}
