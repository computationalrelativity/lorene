/*
 *   Copyright (c) 2005 Philippe Grandclément
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


char et_bin_nsbh_axe_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:52:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:08  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2005/08/31 09:53:58  m_saijo
 * Delete unnecessary words around headers
 *
 * Revision 1.2  2005/08/31 09:13:45  p_grandclement
 * add math.h
 *
 * Revision 1.1  2005/08/29 15:10:16  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 *
 * $Header$
 *
 */

// C Headers
#include <cmath>

// Headers Lorene
#include "et_bin_nsbh.h"
#include "graphique.h"

namespace Lorene {
double Et_bin_nsbh::compute_axe(double omega) const {

  // On récupère les trucs qui vont servir :
  Cmp dx_mu (nnn().dsdx() / nnn());
  Cmp dx_loggam (loggam().dsdx()) ;
  double part_dx = dx_mu(0,0,0,0) + dx_loggam(0,0,0,0) ;
  
  Cmp xabs (mp) ;
  xabs = mp.xa ;
  
  // Derivee_square
  Cmp C_cmp (-pow(confpsi(), 4)/nnn()/nnn()*omega*omega) ;
  C_cmp.std_base_scal() ;
  double dC = C_cmp.dsdx()(0,0,0,0) ;
 
  // Derive single
  Cmp B_cmp (2*pow(confpsi(), 4)/nnn()/nnn()*shift(1)*omega) ;
  B_cmp.std_base_scal() ;
  double dB = B_cmp.dsdx()(0,0,0,0) ;
   
  // Derive const
  Cmp A_cmp = 1-pow(confpsi(), 4)/nnn()/nnn()*(shift(0)*shift(0) + shift(1)*shift(1) + shift(2)*shift(2)) ;
  A_cmp.std_base_scal() ;
  double dA = A_cmp.dsdx()(0,0,0,0) ;

  // coefficients de G
  double A = A_cmp (0,0,0,0) ;
  double B = B_cmp (0,0,0,0) ;
  double C = C_cmp (0,0,0,0) ;
  
  // Les coefficients du polynome :
  double a_coef = dC + 2*C*part_dx ;
  double b_coef = dB + 2*C + 2*B*part_dx ;
  double c_coef = dA  + B + 2*A*part_dx;
 
  cout << a_coef << " " << b_coef << " " << c_coef << endl ;
  
  double determinant = b_coef*b_coef - 4*a_coef*c_coef ;
  
  if (determinant <0) {
      cout << "No solution for Xabs found in Et_bin_nsbh_compute_axe !" << endl ;
      abort() ;
  }
  
  double sol_un = (-b_coef - sqrt(determinant))/2/a_coef ;
  double sol_deux = (-b_coef + sqrt(determinant))/2/a_coef ;
  
  bool signe_un = (sol_un >0) ? true : false ;
  bool signe_deux = (sol_deux >0) ? true : false ;
    
  double res ;
  
  if (signe_un == signe_deux) {
      cout << "To many solutions in Et_bin_nsbh_compute_axe !" << endl ;
      abort() ;
  }
  else {
      res = (signe_un) ? sol_un : sol_deux ;
      }

  return res ;
}
}
