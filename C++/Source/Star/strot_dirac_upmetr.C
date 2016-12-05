/*
 *  Function Star_rot_Dirac::update_metric
 *
 *    (see file star_rot_dirac.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:16  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:40  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2005/03/25 13:47:26  j_novak
 * Added the update of log(Q).
 *
 * Revision 1.2  2005/02/17 17:30:55  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.1  2005/01/31 08:51:48  j_novak
 * New files for rotating stars in Dirac gauge (still under developement).
 *
 *
 * $Header$
 *
 */


// Lorene headers
#include "star_rot_dirac.h"

#include "utilitaires.h"
#include "unites.h" 

namespace Lorene {
void Star_rot_Dirac::update_metric(){

  // Lapse function 
  // ---------------

  nn = exp( logn ) ;

  nn.std_spectral_base() ; // set the bases for spectral expansions


  // Quantity log(Q)
  //----------------

  lnq = log(qqq) ;
  lnq.std_spectral_base() ;

  // Comformal factor $\Psi^4$
  // -------------------------

  psi4 = (qqq * qqq) / (nn * nn) ;

  psi4.std_spectral_base() ;

  // Factor $\Psi^2$
  // ----------------

  psi2 = sqrt( psi4 ) ;

  psi2.std_spectral_base() ;

  // Quantity $ln( \Psi )$
  // ---------------------

  ln_psi = 0.5*log( psi2 ) ;

  ln_psi.std_spectral_base() ;

  // Conformal metric $\tilde{\gamma}$ 
  // --------------------------------------

  tgamma = flat.con() + hh ;   // contravariant representation 
                               // $\tilde{\gamma}^{ij}$ 

  // Physical metric $\gamma$
  // -----------------------------

  gamma = tgamma.con() / psi4 ;    // contravariant representation
                                 //  $\gamma^{ij}$


  // Quantities $A^{ij}$, $\tilde{A}_{ij}, and $\tilde{A}_{ij} A^{ij}$
  // -----------------------------------------------------------------

  aa = ( beta.ope_killing_conf(flat) - hh.derive_lie(beta)) 
    / ( 2*nn ) ;

  taa = aa.up_down(tgamma) ;

  aa_quad = contract(taa, 0, 1, aa, 0, 1) ;

  aa_quad.std_spectral_base() ;


  // The derived quantities are no longer up to date :
  // ------------------------------------------------

  del_deriv() ;

}



}
