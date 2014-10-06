/*
 *  Methods of class Hole_bhns to compute Killing vectors
 *
 *    (see file hole_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2006-2007 Keisuke Taniguchi
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

char hole_bhns_killing_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:10  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/07/02 21:01:11  k_taniguchi
 * A bug removed.
 *
 * Revision 1.1  2008/05/15 19:09:53  k_taniguchi
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
#include "hole_bhns.h"
#include "unites.h"
#include "utilitaires.h"

               //---------------------------------------------//
               //          Killing vectors on the AH          //
               //---------------------------------------------//

Vector Hole_bhns::killing_vect(const Tbl& xi_i, const double& phi_i,
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

	    xi = runge_kutta_phi(xi_ini, pp_0, nrk_phi) ;

	    xi_t.set(i) = xi(0) ;
	    xi_p.set(i) = xi(1) ;
	    xi_l.set(i) = xi(2) ;

	    // New data for the next step
	    //  -------------------------
	    pp_0 += dp ;  // New angle
	    xi_ini = xi ;

	}

	cout << "Hole_bhns::killing_vect :" << endl ;
	cout.precision(16) ;
	cout << "   xi_p(0) :  " << xi_p(0) << endl ;
	cout << "   xi_p(np) : " << xi_p(np) << endl ;

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

	// Normalization of the Killing vector
	Scalar rr(mp) ;
	rr = mp.r ;
	rr.std_spectral_base() ;

	Scalar st(mp) ;
	st = mp.sint ;
	st.std_spectral_base() ;

	Scalar source_phi(mp) ;
	source_phi = pow(confo_tot, 2.) * rr * st / xi_phi ;
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

	cout << "Hole_bhns:: t_f = " << integ << endl ;
	double ratio = 0.5 * integ / M_PI ;

	cout << "Hole_bhns:: t_f / 2M_PI = " << ratio << endl ;

	for (int k=0; k<np; k++) {
	  xi_p.set(k) = xi_phi.val_grid_point(1, k, nt-1, 0) * ratio ;
	}


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

	        xi = runge_kutta_theta(xi_ini, tt_0, pp, nrk_theta) ;

		xi_th.set(nt-1-j, i) = xi(0) ;
		xi_ph.set(nt-1-j, i) = xi(1) ;
		xi_ll.set(nt-1-j, i) = xi(2) ;

		// New data for the nxt step
		// -------------------------
		tt_0 -= dt ;  // New angle
		xi_ini = xi ;

	    }  // End of the loop to the polar direction

	}  // End of the loop to the azimuhtal direction

	cout << "Hole_bhns::killing_vect :" << endl ;
	cout.precision(16) ;
	cout << "   xi_ph(nt-1,0) :  " << xi_ph(nt-1,0) << endl ;
	cout << "   xi_ph(0,0) :     " << xi_ph(0,0) << endl ;
	cout << "   xi_ph(0,np/4) :  " << xi_ph(0,np/4) << endl ;
	cout << "   xi_ph(0,np/2) :  " << xi_ph(0,np/2) << endl ;
	cout << "   xi_ph(0,3np/4) : " << xi_ph(0,3*np/4) << endl ;


	// Construction of the Killing vector
	// ----------------------------------

	// Definition of the vector is at the top of this code
	killing.set_etat_qcq() ;
	killing.set(1) = 0.5 ;  // initialization of L instead of
	                        //  the r component
	killing.set(2) = 0.5 ;  // initialization of the theta component
	killing.set(3) = 0.5 ;  // initialization of the phi component

	for (int l=0; l<2; l++) {
	  for (int i=0; i<nr; i++) {
	    for (int j=0; j<nt; j++) {
	      for (int k=0; k<np; k++) {
		(killing.set(1)).set_grid_point(l, k, j, i) = xi_ll(j, k) ;
		(killing.set(2)).set_grid_point(l, k, j, i) = xi_th(j, k) ;
		(killing.set(3)).set_grid_point(l, k, j, i) = xi_ph(j, k) ;
	      }
	    }
	  }
	}
	killing.std_spectral_base() ;

	// Check the normalization
	// -----------------------

	double check_norm1 = 0. ;
	double check_norm2 = 0. ;
	source_phi = pow(confo_tot, 2.) * rr * st / killing(3) ;
	source_phi.std_spectral_base() ;

	for (int i=0; i<mm; i++) {

	    t1 = hh * double(4*i) ;
	    t2 = hh * double(4*i+1) ;
	    t3 = hh * double(4*i+2) ;
	    t4 = hh * double(4*i+3) ;
	    t5 = hh * double(4*i+4) ;

	    check_norm1 += (hh/45.) *
	      ( 14.*source_phi.val_point(rah,M_PI/4.,t1)
		+ 64.*source_phi.val_point(rah,M_PI/4.,t2)
		+ 24.*source_phi.val_point(rah,M_PI/4.,t3)
		+ 64.*source_phi.val_point(rah,M_PI/4.,t4)
		+ 14.*source_phi.val_point(rah,M_PI/4.,t5) ) ;

	    check_norm2 += (hh/45.) *
	      ( 14.*source_phi.val_point(rah,M_PI/8.,t1)
		+ 64.*source_phi.val_point(rah,M_PI/8.,t2)
		+ 24.*source_phi.val_point(rah,M_PI/8.,t3)
		+ 64.*source_phi.val_point(rah,M_PI/8.,t4)
		+ 14.*source_phi.val_point(rah,M_PI/8.,t5) ) ;

	}

	cout.precision(16) ;
	cout << "Hole_bhns:: t_f for M_PI/4 = " << check_norm1 / M_PI
	     << " M_PI" << endl ;
	cout << "Hole_bhns:: t_f for M_PI/8 = " << check_norm2 / M_PI
	     << " M_PI" << endl ;

	// Checking the accuracy of the Killing vector
	// -------------------------------------------

	// xi^i D_i L = 0
	Scalar dldt(mp) ;
	dldt = killing(1).dsdt() ;

	Scalar dldp(mp) ;
	dldp = killing(1).stdsdp() ;

	Scalar xidl(mp) ;
	xidl = killing(2) * dldt + killing(3) * dldp ;
	xidl.std_spectral_base() ;

	double xidl_error1 = 0. ;
	double xidl_error2 = 0. ;
	double xidl_norm = 0. ;

	for (int j=0; j<nt; j++) {
	  for (int k=0; k<np/2; k++) {
	    xidl_error1 += xidl.val_grid_point(1, k, j, 0) ;
	  }
	}

	for (int j=0; j<nt; j++) {
	  for (int k=np/2; k<np; k++) {
	    xidl_error2 += xidl.val_grid_point(1, k, j, 0) ;
	  }
	}

	for (int j=0; j<nt; j++) {
	  for (int k=0; k<np; k++) {
	    xidl_norm += fabs(xidl.val_grid_point(1, k, j, 0)) ;
	  }
	}

	cout.precision(6) ;
	cout << "Relative error in the 1st condition : "
	     << xidl_error1 / xidl_norm / double(nt) / double(np) * 2.
	     << "  "
	     << xidl_error2 / xidl_norm / double(nt) / double(np) * 2.
	     << "  "
	     << (xidl_error1+xidl_error2)/xidl_norm/double(nt)/double(np)
	     << endl ;

	// D_i xi^i = 0
	Scalar xitst(mp) ;
	xitst = pow(confo_tot, 2.) * rr * st * killing(2) ;
	xitst.std_spectral_base() ;

	Scalar dxitstdt(mp) ;
	dxitstdt = xitst.dsdt() ;

	Scalar xip(mp) ;
	xip = pow(confo_tot, 2.) * rr * killing(3) ;
	xip.std_spectral_base() ;

	Scalar dxipdp(mp) ;
	dxipdp = xip.stdsdp() ;

	Scalar dxi(mp) ;
	dxi = dxitstdt + st * dxipdp ;
	dxi.std_spectral_base() ;

	double dxi_error1 = 0. ;
	double dxi_error2 = 0. ;
	double dxi_norm = 0. ;

	for (int j=0; j<nt; j++) {
	  for (int k=0; k<np/2; k++) {
	    dxi_error1 += dxi.val_grid_point(1, k, j, 0) ;
	  }
	}

	for (int j=0; j<nt; j++) {
	  for (int k=np/2; k<np; k++) {
	    dxi_error2 += dxi.val_grid_point(1, k, j, 0) ;
	  }
	}

	for (int j=0; j<nt; j++) {
	  for (int k=0; k<np; k++) {
	    dxi_norm += fabs(dxi.val_grid_point(1, k, j, 0)) ;
	  }
	}

	cout << "Relative error in the 2nd condition : "
	     << dxi_error1 / dxi_norm / double(nt) / double(np) * 2.
	     << "  "
	     << dxi_error2 / dxi_norm / double(nt) / double(np) * 2.
	     << "  "
	     << (dxi_error1+dxi_error2)/dxi_norm/double(nt)/double(np)
	     << endl ;

	// e^{ij} D_i \xi_j = 2L
	Scalar xipst(mp) ;
	xipst = pow(confo_tot, 2.) * rr * st * killing(3) ;
	xipst.std_spectral_base() ;

	Scalar dxipstdt(mp) ;
	dxipstdt = xipst.dsdt() ;

	Scalar xit(mp) ;
	xit = pow(confo_tot, 2.) * rr * killing(2) ;
	xit.std_spectral_base() ;

	Scalar dxitdp(mp) ;
	dxitdp = xit.stdsdp() ;

	Scalar dxi2l(mp) ;
	dxi2l = dxipstdt - st * dxitdp
	  - 2. * killing(1) * pow(confo_tot, 4.) * rr * rr * st ;
	dxi2l.std_spectral_base() ;

	double dxi2l_error1 = 0. ;
	double dxi2l_error2 = 0. ;
	double dxi2l_norm = 0. ;

	for (int j=0; j<nt; j++) {
	  for (int k=0; k<np/2; k++) {
	    dxi2l_error1 += dxi2l.val_grid_point(1, k, j, 0) ;
	  }
	}

	for (int j=0; j<nt; j++) {
	  for (int k=np/2; k<np; k++) {
	    dxi2l_error2 += dxi2l.val_grid_point(1, k, j, 0) ;
	  }
	}

	for (int j=0; j<nt; j++) {
	  for (int k=0; k<np; k++) {
	    dxi2l_norm += fabs(dxi2l.val_grid_point(1, k, j, 0)) ;
	  }
	}

	cout << "Relative error in the 3rd condition : "
	     << dxi2l_error1 / dxi2l_norm / double(nt) / double(np) * 2.
	     << "  "
	     << dxi2l_error2 / dxi2l_norm / double(nt) / double(np) * 2.
	     << "  "
	     << (dxi2l_error1+dxi2l_error2)/dxi2l_norm/double(nt)/double(np)
	     << endl ;

	cout.precision(4) ;

    }  // End of the loop for isotropic coordinates

    return killing ;

}
