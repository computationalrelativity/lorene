/*
 *  Method of class Et_bin_bhns_extr to compute the extrinsic curvature tensor
 *  in the Kerr-Schild background metric with extreme mass ratio
 *
 *    (see file et_bin_bhns_extr.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Keisuke Taniguchi
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

char et_bin_bhns_extr_excurv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/11/30 20:49:13  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "et_bin_bhns_extr.h"
#include "etoile.h"
#include "coord.h"
#include "unites.h"

void Et_bin_bhns_extr::extrinsic_curv_extr(const double& mass,
					   const double& sepa) {

  using namespace Unites ;

    /** In the codes related BH-NS binary systems with extreme mass ratio,
     *  "tkij_auto" means
     *  A_{tot}^2 A_{NS}^{ij} =-0.5/nnn*(DN_{NS}+DN_{NS}-2/3gammaDN_{NS})
     *  and "akcar_auto"
     *  A_{tot}^2 A^{NS}_{ij} A_{NS}^{ij} =A_{tot}^2 gamma gamma
     *                                     tkij_auto tkij_auto
     */

    // Components of shift_auto with respect to the Cartesian triad
    //  (d/dx, d/dy, d/dz) of the mapping :

    Tenseur shift_auto_local = shift_auto ;
    shift_auto_local.change_triad( mp.get_bvect_cart() ) ;

    // Gradient (partial derivatives with respect to the Cartesian
    //           coordinates of the mapping)
    // dn(i, j) = D_i N^j

    Tenseur dn = shift_auto_local.gradient() ;

    // Return to the absolute reference frame
    dn.change_triad(ref_triad) ;

    // Trace of D_i N^j = divergence of N^j :
    Tenseur divn = contract(dn, 0, 1) ;

    // Computation of quantities coming from the companion (Kerr-Schild BH)
    // --------------------------------------------------------------------

    const Coord& xx = mp.x ;
    const Coord& yy = mp.y ;
    const Coord& zz = mp.z ;

    Tenseur r_bh(mp) ;
    r_bh.set_etat_qcq() ;
    r_bh.set() = pow( (xx+sepa)*(xx+sepa) + yy*yy + zz*zz, 0.5) ;
    r_bh.set_std_base() ;

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

    Tenseur lapse_bh2(mp) ;  // lapse_bh * lapse_bh
    lapse_bh2 = 1. / (1.+2.*msr) ;
    lapse_bh2.set_std_base() ;

    // Computation of some auxiliary functions
    // ---------------------------------------

    shift_auto_local.change_triad( ref_triad ) ;

    Tenseur tmp1(mp, 2, CON, ref_triad) ;
    Tenseur tmp2(mp, 2, CON, ref_triad) ;
    Tenseur tmp3(mp, 2, CON, ref_triad) ;
    tmp1.set_etat_qcq() ;
    tmp2.set_etat_qcq() ;
    tmp3.set_etat_qcq() ;

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
	    tmp1.set(i, j) = -2.*lapse_bh2()%msr()%xsr_con(i)%xsr_con(j) ;

	    tmp2.set(i, j) = -3.*lapse_bh2()%xsr_con(i)%xsr_con(j)
	      -4.*lapse_bh2()*msr()%xsr_con(i)%xsr_con(j) ;

	    tmp3.set(i, j) = xsr_con(i)%shift_auto_local(j) ;
	}
    }

    tmp1.set_std_base() ;
    tmp2.set_std_base() ;
    tmp3.set_std_base() ;

    Tenseur tmp4(mp) ;
    tmp4.set_etat_qcq() ;
    tmp4.set() = 0 ;
    tmp4.set_std_base() ;

    for (int i=0; i<3; i++)
        tmp4.set() += xsr_con(i) % shift_auto_local(i) ;

    tmp4.set_std_base() ;

    // Computation of contraction
    // --------------------------

    Tenseur tmp1dn = contract(tmp1, 1, dn, 0) ;

    // Computation of A^2 A^{ij}
    // -------------------------
    tkij_auto.set_etat_qcq() ;

    for (int i=0; i<3; i++) {
        for (int j=i; j<3; j++) {
	    tkij_auto.set(i, j) = dn(i, j) + dn(j, i)
	      + tmp1dn(i, j) + tmp1dn(j, i)
	      + 2.*lapse_bh2()%msr()/r_bh()%( tmp3(i, j) + tmp3(j, i)
					      + tmp4() % tmp2(i, j) )
	      - double(2)/double(3) * tmp1(i, j)
	      * (divn() - lapse_bh2() % msr() / r_bh() % tmp4()) ;
	}
	tkij_auto.set(i, i) -= double(2) /double(3)
	  * (divn() - lapse_bh2() % msr() / r_bh() % tmp4()) ;
    }

    tkij_auto = - 0.5 * tkij_auto / nnn ;

    tkij_auto.set_std_base() ;

    // Computation of A^2 A_{ij} A^{ij}
    // --------------------------------

    Tenseur xx_cov(mp, 1, COV, ref_triad) ;
    xx_cov.set_etat_qcq() ;
    xx_cov.set(0) = xx + sepa ;
    xx_cov.set(1) = yy ;
    xx_cov.set(2) = zz ;
    xx_cov.set_std_base() ;

    Tenseur xsr_cov(mp, 1, COV, ref_triad) ;
    xsr_cov = xx_cov / r_bh ;
    xsr_cov.set_std_base() ;

    Tenseur tmp5(mp, 2, COV, ref_triad) ;
    Tenseur tmp6(mp, 2, CON, ref_triad) ;
    tmp5.set_etat_qcq() ;
    tmp6.set_etat_qcq() ;

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
	    tmp6.set(i, j) = 0 ;
	}
    }
    tmp6.set_std_base() ;

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
	    tmp5.set(i, j) = 2.*msr()%xsr_cov(i)%xsr_cov(j) ;
	}
    }

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
	    for (int k=0; k<3; k++) {
	        tmp6.set(i, j) += tkij_auto(i,k) % tkij_auto(j,k) ;
	    }
	}
    }

    tmp5.set_std_base() ;
    tmp6.set_std_base() ;

    Tenseur tmp7(mp) ;
    tmp7.set_etat_qcq() ;
    tmp7.set() = 0 ;
    tmp7.set_std_base() ;

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
	    tmp7.set() += tmp5(i,j) % tmp6(i,j) ;
	}
    }
    tmp7.set_std_base() ;

    Tenseur tmp8(mp) ;
    tmp8.set_etat_qcq() ;
    tmp8.set() = 0 ;
    tmp8.set_std_base() ;

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
	    tmp8.set() += tmp5(i,j) % tkij_auto(i,j) ;
	}
    }
    tmp8.set_std_base() ;

    akcar_auto.set_etat_qcq() ;
    akcar_auto.set() = 2.*tmp7() + tmp8() % tmp8() ;
    akcar_auto.set_std_base() ;

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
	    akcar_auto.set() += tkij_auto(i, j) % tkij_auto(i, j) ;
	}
    }

    akcar_auto.set_std_base() ;

    akcar_auto = a_car % akcar_auto ;

}
