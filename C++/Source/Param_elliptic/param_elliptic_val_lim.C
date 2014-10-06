/*
 *   Copyright (c) 2004 Philippe Grandclement
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

char param_elliptic_val_lim_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:13:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/08/24 09:14:49  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 
 * $Header$
 *
 */

#include "headcpp.h"

#include <cmath>
#include <cstdlib>

#include "base_val.h"
#include "param_elliptic.h"
#include "proto.h"

double Param_elliptic::F_plus (int zone, int k, int j) const {

  if (done_F(zone, k, j) == 0)
    compute_val_F(zone, k, j) ;

  return val_F_plus (zone, k, j) ;
}

double Param_elliptic::F_minus (int zone, int k, int j) const {

  if (done_F(zone, k, j) == 0)
    compute_val_F(zone, k, j) ;

  return val_F_minus (zone, k, j) ;
}

double Param_elliptic::dF_plus (int zone, int k, int j) const {

  if (done_F(zone, k, j) == 0)
    compute_val_F(zone, k, j) ;

  return val_dF_plus (zone, k, j) ;
}

double Param_elliptic::dF_minus (int zone, int k, int j) const {

  if (done_F(zone, k, j) == 0)
    compute_val_F(zone, k, j) ;

  return val_dF_minus (zone, k, j) ;
}


double Param_elliptic::G_plus (int zone) const {

  if (done_G(zone) == 0)
    compute_val_G(zone) ;

  return val_G_plus (zone) ;
}

double Param_elliptic::G_minus (int zone) const {

  if (done_G(zone) == 0)
    compute_val_G(zone) ;

  return val_G_minus (zone) ;
}

double Param_elliptic::dG_plus (int zone) const {

  if (done_G(zone) == 0)
    compute_val_G(zone) ;

  return val_dG_plus (zone) ;
}

double Param_elliptic::dG_minus (int zone) const {

  if (done_G(zone) == 0)
    compute_val_G(zone) ;

  return val_dG_minus (zone) ;
}


void Param_elliptic::compute_val_F (int zone, int k, int j) const {
  
  int nr = get_mp().get_mg()->get_nr(zone) ;
  Tbl coefs (nr) ;
  coefs.set_etat_qcq() ;

  bool zero ;

  if (var_F.get_spectral_va().c_cf->get_etat() == ETATZERO)
    zero = true ;
  else
    if ((*var_F.get_spectral_va().c_cf)(zone).get_etat() == ETATZERO)
      zero = true ;
    else 
      zero = false ;

  if (zero)
    coefs.annule_hard() ;
  else
    for (int i=0 ; i<nr ; i++)
      coefs.set(i) = (*var_F.get_spectral_va().c_cf)(zone, k, j, i) ;

  int lq, mq ;
  int base_r ;
  var_F.get_spectral_va().base.give_quant_numbers (zone, k,j, lq, mq, base_r) ;

  double alpha = get_alpha(zone) ;

  Tbl output (val_solp(coefs, alpha, base_r)) ;

  // On range :
  val_F_plus.set(zone, k, j) = output(0) ;
  val_F_minus.set(zone, k, j) = output(1) ;
  val_dF_plus.set(zone, k, j) = output(2) ;
  val_dF_minus.set(zone, k, j) = output(3) ;
  done_F.set(zone, k, j) = 1 ;

}

void Param_elliptic::compute_val_G (int zone) const {
  
  int nr = get_mp().get_mg()->get_nr(zone) ;
  Tbl coefs (nr) ;
  coefs.set_etat_qcq() ;
 

  bool zero ;
  if (var_G.get_spectral_va().c_cf->get_etat() == ETATZERO)
    zero = true ;
  else
    if ((*var_G.get_spectral_va().c_cf)(zone).get_etat() == ETATZERO)
      zero = true ;
    else 
      zero = false ;
  if (zero)
    coefs.annule_hard() ;
  else
    for (int i=0 ; i<nr ; i++)
      coefs.set(i) = (*var_G.get_spectral_va().c_cf)(zone, 0, 0, i) ;
  
  int lq, mq ;
  int base_r ;
  var_G.get_spectral_va().base.give_quant_numbers (zone, 0, 0, lq, mq, base_r) ;


  double alpha = get_alpha(zone) ;

  Tbl output (val_solp(coefs, alpha, base_r)) ;

  // On range :
  val_G_plus.set(zone) = output(0) ;
  val_G_minus.set(zone) = output(1) ;
  val_dG_plus.set(zone) = output(2) ;
  val_dG_minus.set(zone) = output(3) ;
  done_G.set(zone) = 1 ;
}
