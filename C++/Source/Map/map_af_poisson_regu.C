/*
 * Method of the class Map_af for the resolution of the scalar Poisson
 *  equation by using regularized source.
 */

/*
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char map_af_poisson_regu_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/12/19 16:21:43  j_novak
 * Shadow hunt
 *
 * Revision 1.2  2003/10/03 15:58:48  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.13  2000/10/25  16:05:11  keisuke
 * Remake for the arbitrary regularization degree (k_div).
 *
 * Revision 2.12  2000/10/06  15:31:41  keisuke
 * Suppress assertions for nilm_cos and nilm_sin.
 * Add warnings for nilm_cos and nilm_sin.
 *
 * Revision 2.11  2000/09/25  15:02:48  keisuke
 * modify the derivative duu_div.
 *
 * Revision 2.10  2000/09/11  13:50:57  keisuke
 * Set the basis of duu_div to the spherical one.
 *
 * Revision 2.9  2000/09/11  10:15:06  keisuke
 * Change the method to reconstruct source_regu.
 *
 * Revision 2.8  2000/09/09  14:51:20  keisuke
 * Suppress uu_regu.set_dzpuis(0).
 *
 * Revision 2.7  2000/09/07  15:29:50  keisuke
 * Add a new argument Cmp& uu.
 *
 * Revision 2.6  2000/09/06  10:28:42  keisuke
 * Modify the scaling for derivatives.
 *
 * Revision 2.5  2000/09/04  15:55:38  keisuke
 * Include the polar and azimuthal parts of duu_div.
 *
 * Revision 2.4  2000/09/04  13:08:39  keisuke
 * Change alpha[0] into mp_radial->dxdr.
 *
 * Revision 2.3  2000/09/01  08:55:55  keisuke
 * Change val_r into alpha[0].
 *
 * Revision 2.2  2000/08/31  15:59:51  keisuke
 * Modify the arguments.
 *
 * Revision 2.1  2000/08/28  16:11:44  keisuke
 * Add "int nzet" in the argumant.
 *
 * Revision 2.0  2000/08/25  08:48:36  keisuke
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers C
#include <math.h>

// Header Lorene
#include "tenseur.h"
#include "matrice.h"
#include "param.h"
#include "proto.h"

//******************************************************************

void Map_af::poisson_regular(const Cmp& source, int k_div, int nzet,
			     double unsgam1, Param& par, Cmp& uu,
			     Cmp& uu_regu, Cmp& uu_div, Tenseur& duu_div,
			     Cmp& source_regu, Cmp& source_div) const {


    assert(source.get_etat() != ETATNONDEF) ;
    assert(source.get_mp()->get_mg() == mg) ;
    assert(k_div > 0) ;

    double aa = unsgam1 ;  // exponent of the specific enthalpy

    int nzm1 = mg->get_nzone() - 1; 
    int nr = mg->get_nr(0) ;
    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    assert(nr-k_div > 0) ;

    // --------------------------------------------
    //      Expansion of "source" by T_i Y_l^m
    // --------------------------------------------

    // Expansion by cos(j \theta) and sin(j \theta) for theta direction
    // ----------------------------------------------------------------

    const Valeur& sourva = source.va ;
    assert(sourva.get_etat() == ETATQCQ) ;

    Valeur rho(sourva.get_mg()) ;
    sourva.coef() ;
    rho = *(sourva.c_cf) ;

    // Expansion by Legendre function for theta direction
    // --------------------------------------------------

    rho.ylm() ;

    // Obtaining the coefficients of the given source
    // ----------------------------------------------

    Tbl& ccf = *((rho.c_cf)->t[0]) ;

    Tbl nilm_cos(np/2+1, 2*nt, nr) ;
    Tbl nilm_sin(np/2+1, 2*nt, nr) ;

    nilm_cos.set_etat_qcq() ;
    nilm_sin.set_etat_qcq() ;

    for (int k=0; k<=np; k+=2) {

      int m = k / 2 ;

      for (int j=0; j<nt; j++) {

	int l ;
	if(m%2 == 0) {
	  l = 2 * j ;
	}
	else
	  l = 2 * j + 1 ;

	for (int i=0; i<nr; i++) {

	  nilm_cos.set(m, l, i) = ccf(k, j, i) ;
	  nilm_sin.set(m, l, i) = ccf(k+1, j, i) ;

	}
      }
    }

    // -----------------------------------------------------
    //      Expansion of "analytic source" by T_i Y_l^m
    // -----------------------------------------------------

    const Grille3d& grid = *(mg->get_grille3d(0)) ;

    Tbl cf_cil(2*nt, nr, k_div) ;
    cf_cil.set_etat_qcq() ;

    int deg[3] ;
    int dim[3] ;

    deg[0] = 1 ;
    deg[1] = 1 ;
    deg[2] = nr ;
    dim[0] = 1 ;
    dim[1] = 1 ;
    dim[2] = nr ;

    double* tmp1 = new double[nr] ;

    for (int k_deg=1; k_deg<=k_div; k_deg++) {

      for (int l=0; l<2*nt; l++) {

	for (int i=0; i<nr; i++) {

	  double xi = grid.x[i] ;

	  tmp1[i] = (aa + k_deg + 1.) *
		( -(4. * l + 6.) * pow(1. - xi * xi, aa + k_deg) * pow(xi, l)
		  + 4. * (aa + k_deg) * pow(1. - xi * xi, aa + k_deg - 1.) *
		  pow(xi, l + 2.) ) ;

	}

	if (l%2 == 0) {
	  cfrchebp(deg, dim, tmp1, dim, tmp1) ;
	}
	else
	  cfrchebi(deg, dim, tmp1, dim, tmp1) ;

	for (int i=0; i<nr; i++) {

	    cf_cil.set(l, i, k_deg-1) = tmp1[i] ;

	}
      }
    }

    // Calculation of the coefficients : Solve the simultaneous equations
    // -------------------------------

    Tbl alm_cos(np/2+1, 2*nt, k_div) ;
    Tbl alm_sin(np/2+1, 2*nt, k_div) ;

    alm_cos.set_etat_qcq() ;
    alm_sin.set_etat_qcq() ;

    Matrice matrix(k_div, k_div) ;
    matrix.set_etat_qcq() ;

    Tbl rhs_cos(k_div) ;
    Tbl rhs_sin(k_div) ;

    rhs_cos.set_etat_qcq() ;
    rhs_sin.set_etat_qcq() ;

    for (int k=0; k<=np; k+=2) {

      int m = k / 2 ;

      for (int j=0; j<nt; j++) {

	int l ;
	if(m%2 == 0) {
	  l = 2 * j ;
	}
	else
	  l = 2 * j + 1 ;

	for (int i=0; i<k_div; i++) {
	  for (int j2=0; j2<k_div; j2++) {
	    matrix.set(i, j2) = cf_cil(l, nr-1-i, j2) ;
	  }
	}

	matrix.set_band(k_div, k_div) ;

	matrix.set_lu() ;

	for (int i=0; i<k_div; i++) {
	  rhs_cos.set(i) = nilm_cos(m, l, nr-1-i) ;
	  rhs_sin.set(i) = nilm_sin(m, l, nr-1-i) ;
	}

	Tbl sol_cos = matrix.inverse(rhs_cos) ;
	Tbl sol_sin = matrix.inverse(rhs_sin) ;

	for (int i=0; i<k_div; i++) {
	  alm_cos.set(m, l, i) = sol_cos(i) ;
	  alm_sin.set(m, l, i) = sol_sin(i) ;
	}
      }
    }

    // -------------------------------------------------------
    //      Construction of the diverging analytic source
    // -------------------------------------------------------

    source_div.set_etat_qcq() ;
    (source_div.va).set_etat_cf_qcq() ;
    (source_div.va).c_cf->set_etat_qcq() ;
    (source_div.va).c_cf->t[0]->set_etat_qcq() ;

    Valeur& sdva = source_div.va ;
    Mtbl_cf& sdva_cf = *(sdva.c_cf) ;

    // Initialization
    for (int k=0; k<=np; k+=2) {
      for (int j=0; j<nt; j++) {
	for (int i=0; i<nr; i++) {
	  sdva_cf.set(0, k, j, i) = 0 ;
	  sdva_cf.set(0, k+1, j, i) = 0 ;
	}
      }
    }

    for (int k=0; k<=np; k+=2) {

      int m = k / 2  ;

      for (int j=0; j<nt; j++) {

	int l ;

	if (m%2 == 0) {
	  l = 2 * j ;
	}
	else
	  l = 2 * j + 1 ;

	for (int i=0; i<nr; i++) {

	  for (int k_deg=1; k_deg<=k_div; k_deg++) {

	    sdva_cf.set(0, k, j, i) = sdva_cf(0, k, j, i)
	      + alm_cos(m, l, k_deg-1) * cf_cil(l, i, k_deg-1) ;
	    sdva_cf.set(0, k+1, j, i) = sdva_cf(0, k+1, j, i)
	      + alm_sin(m, l, k_deg-1) * cf_cil(l, i, k_deg-1) ;

	  }
	}
      }
    }

    source_div.annule(nzet, nzm1) ;

    Base_val base_std = mg->std_base_scal() ; 

    Base_val& base_s_div = sdva.base ;
    for (int l=0; l<=nzm1; l++) {
      int base_s_div_r = base_std.b[l] & MSQ_R ; 
      int base_s_div_p = base_std.b[l] & MSQ_P ; 

      base_s_div.b[l] = base_s_div_r | T_LEG_P | base_s_div_p ; 	
    }

    sdva_cf.base = base_s_div ;  // copy of the base in the Mtbl_cf

    sdva.ylm_i() ;

    // ------------------------------------------------
    //      Construction of the regularized source
    // ------------------------------------------------

    source_regu.set_etat_qcq() ;
    source_regu = source - source_div ;

    // -------------------------------------------------------------
    //      Solving the Poisson equation for regularized source
    // -------------------------------------------------------------

    source_regu.set_dzpuis(4) ;

    assert(uu_regu.get_mp()->get_mg() == mg) ;

    (*this).poisson(source_regu, par, uu_regu) ;

    // ----------------------------------------------------------
    //      Construction of the diverging analytic potential
    // ----------------------------------------------------------

    Tbl cf_pil(2*nt, nr, k_div) ;
    cf_pil.set_etat_qcq() ;

    double* tmp2 = new double[nr] ;

    for (int k_deg=1; k_deg<=k_div; k_deg++) {

      for (int l=0; l<2*nt; l++) {

	for (int i=0; i<nr; i++) {

	  double xi = grid.x[i] ;
	  tmp2[i] = pow(xi, l) * pow(1. - xi * xi, aa + 1. + k_deg) ;

	}

	if (l%2 == 0) {
	  cfrchebp(deg, dim, tmp2, dim, tmp2) ;
	}
	else
	  cfrchebi(deg, dim, tmp2, dim, tmp2) ;

	for (int i=0; i<nr; i++) {

	  cf_pil.set(l, i, k_deg-1) = tmp2[i] ;

	}
      }
    }

    uu_div.set_etat_qcq() ;
    (uu_div.va).set_etat_cf_qcq() ;
    ((uu_div.va).c_cf)->set_etat_qcq() ;
    ((uu_div.va).c_cf)->t[0]->set_etat_qcq() ;

    Valeur& udva = uu_div.va ;
    Mtbl_cf& udva_cf = *(udva.c_cf) ;

    // Initialization
    for (int k=0; k<=np; k+=2) {
      for (int j=0; j<nt; j++) {
	for (int i=0; i<nr; i++) {
	  udva_cf.set(0, k, j, i) = 0 ;
	  udva_cf.set(0, k+1, j, i) = 0 ;
	}
      }
    }

    for (int k=0; k<=np; k+=2) {

      int m = k / 2  ;

      for (int j=0; j<nt; j++) {

	int l ;

	if (m%2 == 0) {
	  l = 2 * j ;
	}
	else
	  l = 2 * j + 1 ;

	for (int i=0; i<nr; i++) {

	  for (int k_deg=1; k_deg<=k_div; k_deg++) {

	    udva_cf.set(0, k, j, i) = udva_cf(0, k, j, i)
	      + alm_cos(m, l, k_deg-1) * cf_pil(l, i, k_deg-1) ;
	    udva_cf.set(0, k+1, j, i) = udva_cf(0, k+1, j, i)
	      + alm_sin(m, l, k_deg-1) * cf_pil(l, i, k_deg-1) ;

	  }
	}
      }
    }

    uu_div.annule(nzet, nzm1) ;

    Base_val& base_uu_div = (uu_div.va).base ;
    for (int l=0; l<=nzm1; l++) {
      int base_uu_r = base_std.b[l] & MSQ_R ; 
      int base_uu_p = base_std.b[l] & MSQ_P ; 

      base_uu_div.b[l] = base_uu_r | T_LEG_P | base_uu_p ; 	
    }

    udva_cf.base = base_uu_div ;  // copy of the base in the Mtbl_cf

    udva.ylm_i() ;

    // Changing the radial coordinate from "xi" to "r"
    // -----------------------------------------------

    udva = udva * alpha[0] * alpha[0] ;

    // ---------------------------------------------
    //      Construction of the total potential
    // ---------------------------------------------

    uu.set_etat_qcq() ;
    uu = uu_regu + uu_div ;

    // -------------------------------------------------------------------
    //      Construction of the derivative of the diverging potential
    // -------------------------------------------------------------------

    duu_div.set_etat_qcq() ; 

    duu_div.set(0).set_etat_qcq() ; 
    (duu_div.set(0).va).set_etat_cf_qcq() ;
    ((duu_div.set(0).va).c_cf)->set_etat_qcq() ;
    ((duu_div.set(0).va).c_cf)->t[0]->set_etat_qcq() ;

    duu_div.set(1).set_etat_qcq() ;
    (duu_div.set(1).va).set_etat_cf_qcq() ;
    ((duu_div.set(1).va).c_cf)->set_etat_qcq() ;
    ((duu_div.set(1).va).c_cf)->t[0]->set_etat_qcq() ;

    duu_div.set(2).set_etat_qcq() ;
    (duu_div.set(2).va).set_etat_cf_qcq() ;
    ((duu_div.set(2).va).c_cf)->set_etat_qcq() ;
    ((duu_div.set(2).va).c_cf)->t[0]->set_etat_qcq() ;

    Valeur& vr = duu_div.set(0).va ;
    Valeur& vt = duu_div.set(1).va ;
    Valeur& vp = duu_div.set(2).va ;

    Mtbl_cf& vr_cf = *(vr.c_cf) ;
    Mtbl_cf& vt_cf = *(vt.c_cf) ;
    Mtbl_cf& vp_cf = *(vp.c_cf) ;

    // -----------
    // Radial part
    // -----------

    Tbl cf_dril(2*nt, nr, k_div) ;
    cf_dril.set_etat_qcq() ;

    double* tmp3 = new double[nr] ;

    for (int k_deg=1; k_deg<=k_div; k_deg++) {

      for (int i=0; i<nr; i++) {

	double xi = grid.x[i] ;
	tmp3[i] = -2. * (aa + 1. + k_deg) * xi
	              * pow(1. - xi * xi, aa + k_deg) ;

      }

      cfrchebi(deg, dim, tmp3, dim, tmp3) ;

      for (int i=0; i<nr; i++) {

	cf_dril.set(0, i, k_deg-1) = tmp3[i] ;

      }

      for (int l=1; l<2*nt; l++) {

	for (int i=0; i<nr; i++) {

	  double xi = grid.x[i] ;
	  tmp3[i] = l * pow(xi, l - 1.) * pow(1. - xi * xi, aa + 1. + k_deg)
	    -2. * (aa + 1. + k_deg) * pow(xi, l + 1.)
	        * pow(1. - xi * xi, aa + k_deg) ;

	}

	if (l%2 == 0) {
	  cfrchebi(deg, dim, tmp3, dim, tmp3) ;
	}
	else
	  cfrchebp(deg, dim, tmp3, dim, tmp3) ;

	for (int i=0; i<nr; i++) {

	  cf_dril.set(l, i, k_deg-1) = tmp3[i] ;

	}
      }
    }

    // Initialization
    for (int k=0; k<=np; k+=2) {
      for (int j=0; j<nt; j++) {
	for (int i=0; i<nr; i++) {
	  vr_cf.set(0, k, j, i) = 0 ;
	  vr_cf.set(0, k+1, j, i) = 0 ;
	}
      }
    }

    for (int k=0; k<=np; k+=2) {

      int m = k / 2  ;

      for (int j=0; j<nt; j++) {

	int l ;

	if (m%2 == 0) {
	  l = 2 * j ;
	}
	else
	  l = 2 * j + 1 ;

	for (int i=0; i<nr; i++) {

	  for (int k_deg=1; k_deg<=k_div; k_deg++) {

	    vr_cf.set(0, k, j, i) = vr_cf(0, k, j, i)
	      + alm_cos(m, l, k_deg-1) * cf_dril(l, i, k_deg-1) ;
	    vr_cf.set(0, k+1, j, i) = vr_cf(0, k+1, j, i)
	      + alm_sin(m, l, k_deg-1) * cf_dril(l, i, k_deg-1) ;

	  }
	}
      }
    }

    (duu_div.set(0)).annule(nzet, nzm1) ;

    // Reconstruction of the basis of the radial part
    // ----------------------------------------------

    Base_val& base_duu_div_r = vr.base ;
    for (int l=0; l<=nzm1; l++) {
      int base_duu_r_p = base_std.b[l] & MSQ_P ; 

      base_duu_div_r.b[l] = R_CHEBPIM_I | T_LEG_P | base_duu_r_p ; 	
    }

    vr_cf.base = base_duu_div_r ;
    vr.ylm_i() ;
 
    const Coord& RR = dxdr ;

    // Changing the radial coordinate from "xi" to "r"
    // -----------------------------------------------

    Base_val sauve_base( vr.base ) ;
    vr = duu_div(0).va * alpha[0] * alpha[0] * RR ;
    vr.base = sauve_base ;

    // -------------------------
    // Polar and azimuthal parts
    // -------------------------

    Tbl cf_dpil(2*nt, nr, k_div) ;
    cf_dpil.set_etat_qcq() ;

    double* tmp4 = new double[nr] ;

    for (int k_deg=1; k_deg<=k_div; k_deg++) {

      for (int i=0; i<nr; i++) {
	tmp4[i] = 0 ;
      }

      cfrchebi(deg, dim, tmp4, dim, tmp4) ;

      for (int i=0; i<nr; i++) {

	cf_dpil.set(0, i, k_deg-1) = tmp4[i] ;

      }

      for (int l=1; l<2*nt; l++) {

	for (int i=0; i<nr; i++) {

	  double xi = grid.x[i] ;
	  tmp4[i] = pow(xi, l - 1.) * pow(1. - xi * xi, aa + 1. + k_deg) ;

	}

	if (l%2 == 0) {
	  cfrchebi(deg, dim, tmp4, dim, tmp4) ;
	}
	else
	  cfrchebp(deg, dim, tmp4, dim, tmp4) ;

	for (int i=0; i<nr; i++) {

	  cf_dpil.set(l, i, k_deg-1) = tmp4[i] ;

	}
      }
    }

    // Initialization
    for (int k=0; k<=np; k+=2) {
      for (int j=0; j<nt; j++) {
	for (int i=0; i<nr; i++) {
	  vt_cf.set(0, k, j, i) = 0 ;
	  vt_cf.set(0, k+1, j, i) = 0 ;
	  vp_cf.set(0, k, j, i) = 0 ;
	  vp_cf.set(0, k+1, j, i) = 0 ;
	}
      }
    }

    for (int k=0; k<=np; k+=2) {

      int m = k / 2  ;

      for (int j=0; j<nt; j++) {

	int l ;

	if (m%2 == 0) {
	  l = 2 * j ;
	}
	else
	  l = 2 * j + 1 ;

	for (int i=0; i<nr; i++) {

	  for (int k_deg=1; k_deg<=k_div; k_deg++) {

	    vt_cf.set(0, k, j, i) = vt_cf(0, k, j, i)
	      + alm_cos(m, l, k_deg-1) * cf_dpil(l, i, k_deg-1) ;
	    vt_cf.set(0, k+1, j, i) = vt_cf(0, k+1, j, i)
	      + alm_sin(m, l, k_deg-1) * cf_dpil(l, i, k_deg-1) ;

	    vp_cf.set(0, k, j, i) = vp_cf(0, k, j, i)
	      + alm_cos(m, l, k_deg-1) * cf_dpil(l, i, k_deg-1) ;
	    vp_cf.set(0, k+1, j, i) = vp_cf(0, k+1, j, i)
	      + alm_sin(m, l, k_deg-1) * cf_dpil(l, i, k_deg-1) ;

	  }
	}
      }
    }

    (duu_div.set(1)).annule(nzet, nzm1) ;
    (duu_div.set(2)).annule(nzet, nzm1) ;


    // Reconstruction of the basis of the polar part
    // ---------------------------------------------

    Base_val& base_duu_div_p = vt.base ;
    for (int l=0; l<=nzm1; l++) {
      int base_duu_p_p = base_std.b[l] & MSQ_P ;

      base_duu_div_p.b[l] = R_CHEBPIM_I | T_LEG_P | base_duu_p_p ;
    }

    vt_cf.base = base_duu_div_p ;
    vt.ylm_i() ;


    // Reconstruction of the basis of the azimuthal part
    // -------------------------------------------------

    Base_val& base_duu_div_t = vp.base ;
    for (int l=0; l<=nzm1; l++) {
      int base_duu_t_p = base_std.b[l] & MSQ_P ;

      base_duu_div_t.b[l] = R_CHEBPIM_I | T_LEG_P | base_duu_t_p ;
    }

    vp_cf.base = base_duu_div_t ;
    vp.ylm_i() ;


    // Calculation of the derivatives
    // ------------------------------

    vt = (duu_div(1).va).dsdt() ;

    vp = (duu_div(2).va).stdsdp() ;

    // Changing the radial coordinate from "xi" to "r"
    // -----------------------------------------------

    vt = duu_div(1).va * alpha[0] ;

    vp = duu_div(2).va * alpha[0] ;

    // Set the basis of duu_div to the spherical one
    // ---------------------------------------------

    duu_div.set_triad( (*this).get_bvect_spher() ) ;

    delete [] tmp1 ;
    delete [] tmp2 ;
    delete [] tmp3 ;
    delete [] tmp4 ;


}
