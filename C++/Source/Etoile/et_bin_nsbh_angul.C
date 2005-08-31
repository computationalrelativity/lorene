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


char et_bin_nsbh_angul_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
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

#include <math.h>
 
// Headers Lorene
#include "et_bin_nsbh.h"
#include "graphique.h"

double Et_bin_nsbh::compute_angul() const {

  // On récupère les trucs qui vont servir :
  Cmp dx_mu (nnn().dsdx() / nnn());
  Cmp dx_loggam (loggam().dsdx()) ;
  double part_dx = dx_mu(0,0,0,0) + dx_loggam(0,0,0,0) ;
  
  Cmp xabs (mp) ;
  xabs = mp.xa ;
    
  Cmp yabs (mp) ;
  yabs = mp.ya ;
 
  // Derivee_square
  Cmp G_square_cmp (-pow(confpsi(), 4)/nnn()/nnn()*(xabs*xabs+yabs*yabs)) ;
  G_square_cmp.std_base_scal() ;
  double dG_square = G_square_cmp.dsdx()(0,0,0,0) ;
 
  // Derive single
  Cmp G_single_cmp (-2*pow(confpsi(), 4)/nnn()/nnn()*(shift(1)*xabs-shift(0)*yabs)) ;
  G_single_cmp.std_base_scal() ;
  double dG_single = G_single_cmp.dsdx()(0,0,0,0) ;
   
  // Derive const
  Cmp G_const_cmp = 1-pow(confpsi(), 4)/nnn()/nnn()*(shift(0)*shift(0) + shift(1)*shift(1) + shift(2)*shift(2)) ;
  G_const_cmp.std_base_scal() ;
  double dG_const = G_const_cmp.dsdx()(0,0,0,0) ;

  // coefficients de G
  double G_square = G_square_cmp (0,0,0,0) ;
  double G_single = G_single_cmp (0,0,0,0) ;
  double G_const = G_const_cmp (0,0,0,0) ;
  
  // Les coefficients du polynome :
  double a_coef = dG_square + 2*G_square*part_dx ;
  double b_coef = dG_single + 2*G_single*part_dx ;
  double c_coef = dG_const  + 2*G_const *part_dx;
 
  double determinant = b_coef*b_coef - 4*a_coef*c_coef ;
          
  bool soluce = true ;
  double res ;
  
  if (determinant <0) {
     soluce = false ;
  }
  
  else {
  double sol_un = (-b_coef - sqrt(determinant))/2/a_coef ;
  double sol_deux = (-b_coef + sqrt(determinant))/2/a_coef ;
  
  bool signe_un = (sol_un >0) ? true : false ;
  bool signe_deux = (sol_deux >0) ? true : false ;

  
  if (signe_un == signe_deux) {
      soluce = false ;
  }
  else {
      res = (signe_un) ? sol_un : sol_deux ;
      }
}
 
  if (soluce == false)
       return 0 ;
   else 
      return res ;
}
