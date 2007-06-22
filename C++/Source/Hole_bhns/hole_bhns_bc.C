/*
 *  Methods of class Hole_bhns to compute the inner boundary condition
 *  at the excised surface
 *
 *    (see file hile_bhns.h for documentation).
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

char hole_bhns_bc_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
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
#include <math.h>

// Lorene headers
#include "hole_bhns.h"
#include "valeur.h"
#include "grilles.h"
#include "unites.h"

                    //----------------------------------//
                    //     Inner boundary condition     //
                    //----------------------------------//

const Valeur Hole_bhns::bc_lapse() const {

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

    if (bc_lapse_nd) {

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

        if (bc_lapse_fs) {  // dlapse/dr = 0

	    if (kerrschild) {
	        tmp = - d_lapse_comp(1) * st * cp
		  - d_lapse_comp(2) * st * sp - d_lapse_comp(3) * ct
		  - 0.125*sqrt(2.) / rr ;
	    }
	    else {  // Isotropic coordinates with the maximal slicing
	        tmp = - d_lapse_comp(1) * st * cp
		  - d_lapse_comp(2) * st * sp - d_lapse_comp(3) * ct ;
	    }

	}
	else {  // dlapse/dr = 0.25*lapse/rr

	    Scalar tmp1(mp) ;
	    tmp1 = 0.25 * (lapse_auto_rs + lapse_comp) / rr ;
	    tmp1.std_spectral_base() ;
	    tmp1.inc_dzpuis(2) ;  // dzpuis : 0 -> 2

	    if (kerrschild) {  // dlapse/dr = 0.25*lapse/rr
	        tmp = - d_lapse_comp(1) * st * cp
		  - d_lapse_comp(2) * st * sp - d_lapse_comp(3) * ct
		  + tmp1 ;
	    }
	    else {  // Isotropic coordinates with the maximal slicing
	            // dlapse/dr = 0.5*lapse/rr

	        tmp = - d_lapse_comp(1) * st * cp
		  - d_lapse_comp(2) * st * sp - d_lapse_comp(3) * ct
		  + 2. * tmp1 ;
	    }

	}
    }
    else {

        if (bc_lapse_fs) {  // The poisson solver in LORENE assumes
	                    // the asymptotic behavior of the function -> 0

	    if (kerrschild) {
	        tmp = -lapse_comp + 1. - 1./sqrt(2.) ;
		// lapse_auto -> 0.5 <-> lapse_auto_rs -> -0.5
	    }
	    else {  // Isotropic coordinates with the maximal slicing
	        tmp = -lapse_comp + 0.5 ;  // lapse = 0.5

		/*
	        cc = 2. ;
	        tmp = -lapse_comp + 1. - 0.25*cc ;  // lapse = 0.5
		*/
	    }

	}
	else {

	    if (kerrschild) {
	        tmp = -lapse_comp + 0.5 ;
		// lapse_auto -> 1/sqrt(2)
	    }
	    else {  // Isotropic coordinates with the maximal slicing
	        tmp = -lapse_comp + 0.5 ;

		/*
	        cc = 2. * sqrt(2.) ;
	        tmp = -lapse_comp + 1./sqrt(2.) + 0.5 - 0.25*cc ;
		// lapse = 1/sqrt(2)
		*/
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

const Valeur Hole_bhns::bc_shift_x(double ome_orb, double ome_spin,
				   double y_rot) const {

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

        tmp = lapse_auto_bh * st * cp
	  * (lapse_auto_rs + lapse_comp
	     + lapse_auto_bh * (1.-2.*mass*confo_tot*confo_tot/rr))
	  / confo_tot / confo_tot
	  - shift_comp(1)
	  + (ome_orb - ome_spin) * yy + ome_orb * (ori_y_bh - y_rot) ;

	/*
	tmp = ((lapse_auto_rs+lapse_comp)/sqrt(2.)/confo_tot/confo_tot
	       +0.5*(1./confo_tot/confo_tot - 1.)) * st * cp
	  - shift_comp(1)
	  + (ome_orb - ome_spin) * yy + ome_orb * (ori_y_bh - y_rot) ;
	*/
    }
    else {  // Isotropic coordinates with the maximal slicing

        double cc ;

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------

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

	/*
        tmp = ((lapse_tot / confo_tot / confo_tot)
	       - mass*mass*cc/rr/rr/pow(r_are,3.)) * st * cp
	  - shift_comp(1)
	  + (ome_orb - ome_spin) * yy + ome_orb * (ori_y_bh - y_rot) ;
	*/
        tmp = ((lapse_tot / confo_tot / confo_tot)
	       - (0.25*cc/r_are)) * st * cp
	  - shift_comp(1)
	  + (ome_orb - ome_spin) * yy + ome_orb * (ori_y_bh - y_rot) ;

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

const Valeur Hole_bhns::bc_shift_y(double ome_orb, double ome_spin,
				   double x_rot) const {

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

        tmp = lapse_auto_bh * st * sp
	  * (lapse_auto_rs + lapse_comp
	     + lapse_auto_bh * (1.-2.*mass*confo_tot*confo_tot/rr))
	  / confo_tot / confo_tot
	  - shift_comp(2)
	  - (ome_orb - ome_spin) * xx - ome_orb * (ori_x_bh - x_rot) ;

	/*
	tmp = ((lapse_auto_rs+lapse_comp)/sqrt(2.)/confo_tot/confo_tot
	       +0.5*(1./confo_tot/confo_tot - 1.)) * st * sp
	  - shift_comp(2)
	  - (ome_orb - ome_spin) * xx - ome_orb * (ori_x_bh - x_rot) ;
	*/
    }
    else {  // Isotropic coordinates with the maximal slicing

        double cc ;

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------

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

	/*
        tmp = ((lapse_tot / confo_tot / confo_tot)
	       - mass*mass*cc/rr/rr/pow(r_are,3.)) * st * sp
	  - shift_comp(2)
	  - (ome_orb - ome_spin) * xx - ome_orb * (ori_x_bh - x_rot) ;
	*/
        tmp = ((lapse_tot / confo_tot / confo_tot)
	       - (0.25*cc/r_are)) * st * sp
	  - shift_comp(2)
	  - (ome_orb - ome_spin) * xx - ome_orb * (ori_x_bh - x_rot) ;

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

        tmp = lapse_auto_bh * ct
	  * (lapse_auto_rs + lapse_comp
	     + lapse_auto_bh * (1.-2.*mass*confo_tot*confo_tot/rr))
	  / confo_tot / confo_tot
	  - shift_comp(3) ;

	/*
	tmp = ((lapse_auto_rs+lapse_comp)/sqrt(2.)/confo_tot/confo_tot
	       +0.5*(1./confo_tot/confo_tot - 1.)) * ct
	  - shift_comp(3) ;
	*/
    }
    else {  // Isotropic coordinates with the maximal slicing

        double cc ;

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------

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

	/*
        tmp = ((lapse_tot / confo_tot / confo_tot)
	       - mass*mass*cc/rr/rr/pow(r_are,3.)) * ct - shift_comp(3) ;
	*/
        tmp = ((lapse_tot / confo_tot / confo_tot)
	       - (0.25*cc/r_are)) * ct - shift_comp(3) ;

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
    ll.set(1) = st * cp ;
    ll.set(2) = st * sp ;
    ll.set(3) = ct ;
    ll.std_spectral_base() ;

    Scalar divshift(mp) ;  // dzpuis = 2
    divshift = shift_auto_rs(1).deriv(1) + shift_auto_rs(2).deriv(2)
      + shift_auto_rs(3).deriv(3) + d_shift_comp(1,1) + d_shift_comp(2,2)
      + d_shift_comp(3,3) ;
    divshift.std_spectral_base() ;

    Scalar llshift(mp) ;   // dzpuis = 0
    llshift = ll(1) * (shift_auto_rs(1) + shift_comp(1))
      + ll(2) * (shift_auto_rs(2) + shift_comp(2))
      + ll(3) * (shift_auto_rs(3) + shift_comp(3)) ;
    llshift.std_spectral_base() ;

    Scalar llshift_auto_rs(mp) ;   // dzpuis = 0
    llshift_auto_rs = ll(1)*shift_auto_rs(1) + ll(2)*shift_auto_rs(2)
      + ll(3)*shift_auto_rs(3) ;
    llshift_auto_rs.std_spectral_base() ;

    Scalar lldllsh = llshift_auto_rs.dsdr()
      + ll(1) * ( ll(1)*d_shift_comp(1,1) + ll(2)*d_shift_comp(1,2)
		  + ll(3)*d_shift_comp(1,3) )
      + ll(2) * ( ll(1)*d_shift_comp(2,1) + ll(2)*d_shift_comp(2,2)
		  + ll(3)*d_shift_comp(2,3) )
      + ll(3) * ( ll(1)*d_shift_comp(3,1) + ll(2)*d_shift_comp(3,2)
		  + ll(3)*d_shift_comp(3,3) ) ; // dzpuis = 2
    lldllsh.std_spectral_base() ;

    Scalar tmp2 = divshift ;
    Scalar tmp3 = -3.*lldllsh ;

    tmp2.dec_dzpuis(2) ;
    tmp3.dec_dzpuis(2) ;

    Scalar tmp(mp) ;

    double mass = ggrav * mass_bh ;

    if (kerrschild) {

	Scalar tmp1 = 0.5 * confo_tot
	  * (lapse_auto_bh * confo_tot * confo_tot / lapse_tot - 1.) / rr ;

	Scalar tmp4 = 2. * mass * lapse_auto_bh * lapse_auto_bh
	  * llshift / rr / rr ;
	// dzpuis = 0

	Scalar tmp5 = 4. * mass * pow(lapse_auto_bh,3.) * (1.+3.*mass/rr)
	  * (lapse_auto_rs + lapse_comp) / rr / rr ;
	// dzpuis = 0

	Scalar tmp6 = mass * lapse_auto_bh * pow(confo_tot,3.) * ome_orb
	  * ( ll(2) * (mp.get_ori_x() - x_rot)
	      - ll(1) * (mp.get_ori_y() - y_rot) )
	  / 6. / lapse_tot / rr / rr ;
	// dzpuis = 0

	Scalar tmp7 = - ll(1)*d_confo_comp(1) - ll(2)*d_confo_comp(2)
	  - ll(3)*d_confo_comp(3) ;
	tmp7.std_spectral_base() ;
	tmp7.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	tmp = tmp7 + tmp1 + tmp6
	  + pow(confo_tot,3.) * (tmp2 + tmp3 + tmp4 + tmp5)
	  / 12. / lapse_tot / lapse_auto_bh ;

    }
    else {  // Isotropic coordinates with the maximal slicing

        double cc ;

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------

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

        Scalar tmp1 = - 0.5 * (confo_auto_rs + confo_comp) / rr ;
        Scalar tmp7 = - ll(1)*d_confo_comp(1) - ll(2)*d_confo_comp(2)
	  - ll(3)*d_confo_comp(3) ;
	tmp7.std_spectral_base() ;
	tmp7.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	/*
	Scalar tmp8 = 0.5 * sqrt(1. - 2.*mass/r_are/rr
				 + cc*cc*pow(mass/r_are/rr,4.))
	  * (pow(confo_tot,3.)*mass*mass*cc/lapse_tot/pow(r_are*rr,3.)
	     - sqrt(r_are) / rr) ;
	*/
	Scalar tmp8 = 0.125*cc*(0.25*cc*pow(confo_tot,3.)/r_are/lapse_tot
				- sqrt(r_are)) / rr ;
	tmp8.std_spectral_base() ;

        tmp = tmp7 + tmp1
	  + pow(confo_tot,3.) * (tmp2 + tmp3) / 12. / lapse_tot
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
