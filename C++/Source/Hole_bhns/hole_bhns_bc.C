/*
 *  Methods of class Hole_bhns to compute the inner boundary condition
 *  at the excised surface
 *
 *    (see file hile_bhns.h for documentation).
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

char hole_bhns_bc_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:10  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/05/15 19:04:10  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:23:56  k_taniguchi
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
#include "valeur.h"
#include "grilles.h"
#include "unites.h"

                    //----------------------------------//
                    //     Inner boundary condition     //
                    //----------------------------------//

const Valeur Hole_bhns::bc_lapconf() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    Scalar tmp(mp) ;

    //    double cc ; // C/M^2

    if (bc_lapconf_nd) {

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

	Scalar rr(mp) ;
	rr = mp.r ;
	rr.std_spectral_base() ;

        if (bc_lapconf_fs) {  // dlapconf/dr = 0

	    if (kerrschild) {
	        cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
		abort() ;
	    }
	    else {  // Isotropic coordinates with the maximal slicing
	        tmp = - d_lapconf_comp(1) % st % cp
		  - d_lapconf_comp(2) % st % sp - d_lapconf_comp(3) % ct ;
	    }

	}
	else {  // dlapconf/dr = 0.5*lapconf/rr

	    Scalar tmp1(mp) ;
	    tmp1 = 0.5 * (lapconf_auto_rs + lapconf_comp) / rr ;
	    tmp1.std_spectral_base() ;
	    tmp1.inc_dzpuis(2) ;  // dzpuis : 0 -> 2

	    if (kerrschild) {  // dlapconf/dr = 0.5*lapconf/rr
	        cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
		abort() ;
	    }
	    else {  // Isotropic coordinates with the maximal slicing
	            // dlapconf/dr = 0.5*lapconf/rr

	        tmp = - d_lapconf_comp(1) % st % cp
		  - d_lapconf_comp(2) % st % sp - d_lapconf_comp(3) % ct
		  + tmp1 ;
	    }

	}
    }
    else {

        if (bc_lapconf_fs) {  // The poisson solver in LORENE assumes
	                      // the asymptotic behavior of the function -> 0

	    if (kerrschild) {
	        cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
		abort() ;
		// lapconf_auto -> 0.5 <-> lapconf_auto_rs -> -0.5
	    }
	    else {  // Isotropic coordinates with the maximal slicing
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
		//	        tmp = -lapconf_comp + 0.5 ;  // lapconf = 0.5

	    }

	}
	else {

	    if (kerrschild) {
	        cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
		abort() ;
	    }
	    else {  // Isotropic coordinates with the maximal slicing
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
		//	        tmp = -lapconf_comp + 0.5 ;

	    }

	}
    }

    bc = 1. ;
    for (int j=0; j<nt; j++) {
        for (int k=0; k<np; k++) {
	    bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
	}
    }

    bc.std_base_scal() ;
    return bc ;

}

const Valeur Hole_bhns::bc_shift_x(double ome_orb, double y_rot) const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

    Base_val** bases = mp.get_mg()->std_base_vect_cart() ;

    Scalar rr(mp) ;
    rr = mp.r ;
    rr.std_spectral_base() ;
    Scalar st(mp) ;
    st = mp.sint ;
    st.std_spectral_base() ;
    Scalar cp(mp) ;
    cp = mp.cosp ;
    cp.std_spectral_base() ;
    Scalar yy(mp) ;
    yy = mp.y ;
    yy.std_spectral_base() ;

    double mass = ggrav * mass_bh ;
    double ori_y_bh = mp.get_ori_y() ;

    Scalar tmp(mp) ;

    if (kerrschild) {

	// Exact solution of an isolated BH is extracted

        cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	abort() ;

    }
    else {  // Isotropic coordinates with the maximal slicing

        double cc ;

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------

	if (bc_lapconf_nd) {  // Neumann boundary condition
	    if (bc_lapconf_fs) {  // First condition
	        // d(\alpha \psi)/dr = 0
	        // ---------------------
	        cc = 2. * (sqrt(13.) - 1.) / 3. ;
	    }
	    else {  // Second condition
	        // d(\alpha \psi)/dr = (\alpha \psi)/(2 rah)
	        // -----------------------------------------
	        cc = 4. / 3. ;
	    }
	}
	else {  // Dirichlet boundary condition
	    if (bc_lapconf_fs) {  // First condition
	        // (\alpha \psi) = 1/2
	        // -------------------
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
	    }
	    else {  // Second condition
	        // (\alpha \psi) = 1/sqrt(2.) \psi_KS
	        // ----------------------------------
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
		//	        cc = 2. * sqrt(2.) ;
	    }
	}

        Scalar r_are(mp) ;
	r_are = r_coord(bc_lapconf_nd, bc_lapconf_fs) ;
	r_are.std_spectral_base() ;

	/*
        tmp = ((lapse_tot / confo_tot / confo_tot)
	       - mass*mass*cc/rr/rr/pow(r_are,3.)) * st * cp
	  - shift_comp(1)
	  + (ome_orb - omega_spin) * yy + ome_orb * (ori_y_bh - y_rot) ;
	*/
        tmp = ((lapconf_tot / pow(confo_tot,3.)) - (0.25*cc/r_are)) * st * cp
	  - shift_comp(1)
	  + (ome_orb - omega_spin) * yy + ome_orb * (ori_y_bh - y_rot) ;

    }

    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    bc = 1. ;
    for (int j=0; j<nt; j++) {
        for (int k=0; k<np; k++) {
	    bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
	}
    }

    bc.base = *bases[0] ;

    for (int i=0; i<3; i++)
        delete bases[i] ;

    delete [] bases ;

    return bc ;

}

const Valeur Hole_bhns::bc_shift_y(double ome_orb, double x_rot) const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

    Base_val** bases = mp.get_mg()->std_base_vect_cart() ;

    Scalar rr(mp) ;
    rr = mp.r ;
    rr.std_spectral_base() ;
    Scalar st(mp) ;
    st = mp.sint ;
    st.std_spectral_base() ;
    Scalar sp(mp) ;
    sp = mp.sinp ;
    sp.std_spectral_base() ;
    Scalar xx(mp) ;
    xx = mp.x ;
    xx.std_spectral_base() ;

    double mass = ggrav * mass_bh ;
    double ori_x_bh = mp.get_ori_x() ;

    Scalar tmp(mp) ;

    if (kerrschild) {

	// Exact solution of an isolated BH is extracted

        cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	abort() ;

    }
    else {  // Isotropic coordinates with the maximal slicing

        double cc ;

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------

	if (bc_lapconf_nd) {  // Neumann boundary condition
	    if (bc_lapconf_fs) {  // First condition
	        // d(\alpha \psi)/dr = 0
	        // ---------------------
	        cc = 2. * (sqrt(13.) - 1.) / 3. ;
	    }
	    else {  // Second condition
	        // d(\alpha \psi)/dr = (\alpha \psi)/(2 rah)
	        // -----------------------------------------
	        cc = 4. / 3. ;
	    }
	}
	else {  // Dirichlet boundary condition
	    if (bc_lapconf_fs) {  // First condition
	        // (\alpha \psi) = 1/2
	        // -------------------
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
	    }
	    else {  // Second condition
	        // (\alpha \psi) = 1/sqrt(2.) \psi_KS
	        // ----------------------------------
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
		//	        cc = 2. * sqrt(2.) ;
	    }
	}

        Scalar r_are(mp) ;
	r_are = r_coord(bc_lapconf_nd, bc_lapconf_fs) ;
	r_are.std_spectral_base() ;

	/*
        tmp = ((lapse_tot / confo_tot / confo_tot)
	       - mass*mass*cc/rr/rr/pow(r_are,3.)) * st * sp
	  - shift_comp(2)
	  - (ome_orb - omega_spin) * xx - ome_orb * (ori_x_bh - x_rot) ;
	*/
        tmp = ((lapconf_tot / pow(confo_tot,3.)) - (0.25*cc/r_are)) * st * sp
	  - shift_comp(2)
	  - (ome_orb - omega_spin) * xx - ome_orb * (ori_x_bh - x_rot) ;

    }

    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    bc = 1. ;
    for (int j=0; j<nt; j++) {
        for (int k=0; k<np; k++) {
	    bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
	}
    }

    bc.base = *bases[1] ;

    for (int i=0; i<3; i++)
        delete bases[i] ;

    delete [] bases ;

    return bc ;

}

const Valeur Hole_bhns::bc_shift_z() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

    Base_val** bases = mp.get_mg()->std_base_vect_cart() ;

    Scalar rr(mp) ;
    rr = mp.r ;
    rr.std_spectral_base() ;
    Scalar ct(mp) ;
    ct = mp.cost ;
    ct.std_spectral_base() ;

    double mass = ggrav * mass_bh ;

    Scalar tmp(mp) ;

    if (kerrschild) {

	// Exact solution of an isolated BH is extracted

        cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	abort() ;

    }
    else {  // Isotropic coordinates with the maximal slicing

        double cc ;

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------

	if (bc_lapconf_nd) {  // Neumann boundary condition
	    if (bc_lapconf_fs) {  // First condition
	        // d(\alpha \psi)/dr = 0
	        // ---------------------
	        cc = 2. * (sqrt(13.) - 1.) / 3. ;
	    }
	    else {  // Second condition
	        // d(\alpha \psi)/dr = (\alpha \psi)/(2 rah)
	        // -----------------------------------------
	        cc = 4. / 3. ;
	    }
	}
	else {  // Dirichlet boundary condition
	    if (bc_lapconf_fs) {  // First condition
	        // (\alpha \psi) = 1/2
	        // -------------------
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
	    }
	    else {  // Second condition
	        // (\alpha \psi) = 1/sqrt(2.) \psi_KS
	        // ----------------------------------
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
		//	        cc = 2. * sqrt(2.) ;
	    }
	}

        Scalar r_are(mp) ;
	r_are = r_coord(bc_lapconf_nd, bc_lapconf_fs) ;
	r_are.std_spectral_base() ;

	/*
        tmp = ((lapse_tot / confo_tot / confo_tot)
	       - mass*mass*cc/rr/rr/pow(r_are,3.)) * ct - shift_comp(3) ;
	*/
        tmp = ((lapconf_tot / pow(confo_tot,3.)) - (0.25*cc/r_are)) * ct
	  - shift_comp(3) ;

    }

    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    bc = 1. ;
    for (int j=0; j<nt; j++) {
        for (int k=0; k<np; k++) {
	    bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
	}
    }

    bc.base = *bases[2] ;

    for (int i=0; i<3; i++)
      delete bases[i] ;

    delete [] bases ;

    return bc ;

}

const Valeur Hole_bhns::bc_confo(double ome_orb, double x_rot,
				 double y_rot) const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

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
    ll.set(1) = st % cp ;
    ll.set(2) = st % sp ;
    ll.set(3) = ct ;
    ll.std_spectral_base() ;

    Scalar divshift(mp) ;  // dzpuis = 2
    divshift = shift_auto_rs(1).deriv(1) + shift_auto_rs(2).deriv(2)
      + shift_auto_rs(3).deriv(3) + d_shift_comp(1,1) + d_shift_comp(2,2)
      + d_shift_comp(3,3) ;
    divshift.std_spectral_base() ;

    Scalar llshift(mp) ;   // dzpuis = 0
    llshift = ll(1) % (shift_auto_rs(1) + shift_comp(1))
      + ll(2) % (shift_auto_rs(2) + shift_comp(2))
      + ll(3) % (shift_auto_rs(3) + shift_comp(3)) ;
    llshift.std_spectral_base() ;

    Scalar llshift_auto_rs(mp) ;   // dzpuis = 0
    llshift_auto_rs = ll(1)%shift_auto_rs(1) + ll(2)%shift_auto_rs(2)
      + ll(3)%shift_auto_rs(3) ;
    llshift_auto_rs.std_spectral_base() ;

    Scalar lldllsh = llshift_auto_rs.dsdr()
      + ll(1) * ( ll(1)%d_shift_comp(1,1) + ll(2)%d_shift_comp(1,2)
		  + ll(3)%d_shift_comp(1,3) )
      + ll(2) * ( ll(1)%d_shift_comp(2,1) + ll(2)%d_shift_comp(2,2)
		  + ll(3)%d_shift_comp(2,3) )
      + ll(3) * ( ll(1)%d_shift_comp(3,1) + ll(2)%d_shift_comp(3,2)
		  + ll(3)%d_shift_comp(3,3) ) ; // dzpuis = 2
    lldllsh.std_spectral_base() ;

    Scalar tmp2 = divshift ;
    Scalar tmp3 = -3.*lldllsh ;

    tmp2.dec_dzpuis(2) ;
    tmp3.dec_dzpuis(2) ;

    Scalar tmp(mp) ;

    double mass = ggrav * mass_bh ;

    if (kerrschild) {

        cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	abort() ;

    }
    else {  // Isotropic coordinates with the maximal slicing

        double cc ;

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------

	if (bc_lapconf_nd) {  // Neumann boundary condition
	    if (bc_lapconf_fs) {  // First condition
	        // d(\alpha \psi)/dr = 0
	        // ---------------------
	        cc = 2. * (sqrt(13.) - 1.) / 3. ;
	    }
	    else {  // Second condition
	        // d(\alpha \psi)/dr = (\alpha \psi)/(2 rah)
	        // -----------------------------------------
	        cc = 4. / 3. ;
	    }
	}
	else {  // Dirichlet boundary condition
	    if (bc_lapconf_fs) {  // First condition
	        // (\alpha \psi) = 1/2
	        // -------------------
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
	    }
	    else {  // Second condition
	        // (\alpha \psi) = 1/sqrt(2.) \psi_KS
	        // ----------------------------------
	        cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
		abort() ;
		//	        cc = 2. * sqrt(2.) ;
	    }
	}

        Scalar r_are(mp) ;
	r_are = r_coord(bc_lapconf_nd, bc_lapconf_fs) ;
	r_are.std_spectral_base() ;

        Scalar tmp1 = - 0.5 * (confo_auto_rs + confo_comp) / rr ;
        Scalar tmp7 = - ll(1)%d_confo_comp(1) - ll(2)%d_confo_comp(2)
	  - ll(3)%d_confo_comp(3) ;
	tmp7.std_spectral_base() ;
	tmp7.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	/*
	Scalar tmp8 = 0.5 * sqrt(1. - 2.*mass/r_are/rr
				 + cc*cc*pow(mass/r_are/rr,4.))
	  * (pow(confo_tot,3.)*mass*mass*cc/lapse_tot/pow(r_are*rr,3.)
	     - sqrt(r_are) / rr) ;
	*/
	Scalar tmp8 = 0.125*cc*(0.25*cc*pow(confo_tot,4.)/r_are/lapconf_tot
				- sqrt(r_are)) / rr ;
	tmp8.std_spectral_base() ;

        tmp = tmp7 + tmp1
	  + pow(confo_tot,4.) * (tmp2 + tmp3) / 12. / lapconf_tot
	  + tmp8 ;

    }

    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    bc = 1. ;
    for (int j=0; j<nt; j++) {
        for (int k=0; k<np; k++) {
	    bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
	}
    }

    bc.std_base_scal() ;
    return bc ;

}
