/*
 *  Function Star_rot_Dirac::hydro_euler
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

char strot_dirac_hydro_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2005/03/14 14:01:34  j_novak
 * u_euler is now defined on all the grid.
 *
 * Revision 1.2  2005/02/17 17:31:56  f_limousin
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


// C headers
#include<math.h>
#include<assert.h>

// Lorene headers
#include"star_rot_dirac.h"

void Star_rot_Dirac::hydro_euler(){

  // u_euler (fluid 3-velocity w.r.t. the Eulerian frame)
  // -----------------------------------------------------


  u_euler.set(1).set_etat_zero() ;
  u_euler.set(2).set_etat_zero() ;

  u_euler.set(3) = omega ;
  u_euler.set(3).std_spectral_base() ;
  u_euler.set(3).mult_rsint() ;
  u_euler.set(3) += beta(3) ;

  u_euler = u_euler / nn ;

  // v2 (square of the norm of u_euler)
  // ----------------------------------

  v2 = contract(contract(gamma.cov(), 0, u_euler, 0), 0, u_euler, 0) ;


  // gam_euler (Lorentz factor between the fluid and Eulerian observers)
  // -------------------------------------------------------------------

  gam_euler = 1. / sqrt(1. - v2) ; 

  gam_euler.std_spectral_base() ;


  // ener_euler (energy density w.r.t. the Eulerian observer)
  // ------------------------------------------------------

  ener_euler = gam_euler*gam_euler*(ener + press) - press ;

  ener_euler.std_spectral_base() ;


  // j_euler (momentum density 3-vector w.r.t. the Eulerian observer)
  // ----------------------------------------------------------------

  j_euler = (ener_euler + press)*u_euler ;

  j_euler.std_spectral_base() ;


  // s_euler (trace of the stress tensor w.r.t. the Eulerian observer)
  // ----------------------------------------------------------------

  s_euler = (ener_euler + press)*v2 + 3*press ;

  s_euler.std_spectral_base() ;


  // stress_euler (stress tensor w.r.t. the Eulerian observer)
  // ---------------------------------------------------------


  stress_euler = (ener_euler + press)*u_euler*u_euler + press*gamma.con() ;

  stress_euler.std_spectral_base() ;


  // The derived quantities are obsolete
  // ------------------------------------

  del_deriv() ;


}

  
