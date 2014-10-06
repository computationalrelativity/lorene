/*
 * Methods Star_bin_xcts::update_metric_der_comp
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2010 Michal Bejger
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

char star_bin_upmetr_der_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/06 15:13:17  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2010/12/09 10:44:42  m_bejger
 * Cosmetic changes
 *
 * Revision 1.3  2010/06/15 14:58:19  m_bejger
 * Companion D^j beta^i computed in companion frame and imported
 *
 * Revision 1.2  2010/06/15 08:09:43  m_bejger
 * *** empty log message ***
 *
 * Revision 1.1  2010/05/04 07:51:05  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Headers Lorene
#include "star.h"
#include "utilitaires.h"
#include "graphique.h"

void Star_bin_xcts::update_metric_der_comp(const Star_bin_xcts& comp) {

  // Derivatives of metric coefficients
  // ----------------------------------

    // Computation of dcov_Psi
    Vector temp = (comp.Psi_auto).derive_cov(comp.flat) ;

    temp.dec_dzpuis(2) ;
    temp.change_triad(mp.get_bvect_cart()) ;

    assert ( *(temp.get_triad()) == *(dcov_Psi.get_triad())) ;

	for(int i=1; i<=3; i++) {

		dcov_Psi.set(i).import(temp(i)) ;
        dcov_Psi.set(i).set_spectral_va().set_base(temp(i).get_spectral_va().get_base()) ;

	}
    dcov_Psi.inc_dzpuis(2) ;

    dcov_Psi += Psi_auto.derive_cov(flat) ;

    // Computation of dcov_chi
    temp = (comp.chi_auto).derive_cov(comp.flat) ;

    temp.dec_dzpuis(2) ;
    temp.change_triad(mp.get_bvect_cart()) ;

    assert ( *(temp.get_triad()) == *(dcov_chi.get_triad())) ;

    for(int i=1; i<=3; i++) {

		dcov_chi.set(i).import(temp(i)) ;
        dcov_chi.set(i).set_spectral_va().set_base(temp(i).get_spectral_va().get_base()) ;

    }
    dcov_chi.inc_dzpuis(2) ;

    dcov_chi += chi_auto.derive_cov(flat) ;

  // Computation of \hat{A}^{ij}_{comp}
  // ----------------------------------

  // D^j beta^i

  Tensor temp_beta = (comp.beta_auto).derive_con(comp.flat) ;
  temp_beta.dec_dzpuis(2) ;
  temp_beta.change_triad(mp.get_bvect_cart()) ;

  Tensor dbeta_comp(mp, 2, CON, mp.get_bvect_cart()) ;

  assert ( *(temp_beta.get_triad()) == *(dbeta_comp.get_triad()) ) ;

    for(int i=1; i<=3; i++)
	for(int j=1; j<=3; j++) {

        // importing
	    (dbeta_comp.set(i,j)).import( (temp_beta)(i,j) );
		// setting appropriate bases
		(dbeta_comp.set(i,j)).set_spectral_va().set_base(temp_beta(i,j).get_spectral_va().get_base()) ;

	}

   dbeta_comp.inc_dzpuis(2) ;

  // Trace of D_j beta^i  :
  Scalar divbeta_comp = beta_comp.divergence(flat) ;

  for (int i=1; i<=3; i++)
     for (int j=1; j<=i; j++) {
//##
//    for (int j=i; j<=3; j++) {

      haij_comp.set(i, j) = dbeta_comp(i, j) + dbeta_comp(j, i)
	  - double(2) /double(3) * divbeta_comp * (flat.con())(i,j) ;

    }

  // Computation of (\hat{A}_{ij}\hat{A}^{ij})_{comp}
  // ------------------------------------------------

  Sym_tensor haij_auto_cov = haij_auto.up_down(flat) ;

  haij_comp = 0.5 * pow(Psi, 7.) / chi * haij_comp ;
  //## for comparison: old formulation
  //haij_comp = 0.5 * haij_comp / nn ;
  haij_comp.std_spectral_base() ;

  hacar_comp = contract(haij_auto_cov, 0, 1, haij_comp, 0, 1, true) ;

  // The derived quantities are obsolete
  // -----------------------------------

  del_deriv() ;

}
