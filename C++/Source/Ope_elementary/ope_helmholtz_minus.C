/*
 *   Copyright (c) 2003 Philippe Grandclement
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

char ope_helmholtz_minus_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/12/11 14:48:50  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */
#include <math.h>

#include "proto.h"
#include "ope_elementary.h"

// Standard constructor :
Ope_helmholtz_minus::Ope_helmholtz_minus (int nbr, int base, double alf, 
					  double bet, double mas): 
  Ope_elementary(nbr, base, alf, bet), masse (mas) {
}

// Constructor by copy :
Ope_helmholtz_minus::Ope_helmholtz_minus (const Ope_helmholtz_minus& so) : 
  Ope_elementary(so), masse (so.masse) {
}

// Destructor :
Ope_helmholtz_minus::~Ope_helmholtz_minus() {} 

// True functions :
void Ope_helmholtz_minus::do_ope_mat() const {
  if (ope_mat != 0x0)
    delete ope_mat ;

  ope_mat = new Matrice 
    (helmholtz_minus_mat(nr, alpha, beta, masse, base_r)) ;
}

void Ope_helmholtz_minus::do_ope_cl() const {
  if (ope_mat == 0x0)
    do_ope_mat() ;

  if (ope_cl != 0x0)
    delete ope_cl ;

  ope_cl = new Matrice 
    (cl_helmholtz_minus(*ope_mat, alpha, beta, masse, base_r)) ;
}

void Ope_helmholtz_minus::do_non_dege() const {
  if (ope_cl == 0x0)
    do_ope_cl() ;

  if (non_dege != 0x0)
    delete non_dege ;

  non_dege = new Matrice 
    (prepa_helmholtz_minus_nondege(*ope_cl, alpha, beta, masse, base_r)) ;
}
  
Tbl Ope_helmholtz_minus::get_solp (const Tbl& so) const {

  if (non_dege == 0x0)
    do_non_dege() ;

  Tbl res(solp_helmholtz_minus (*ope_mat, *non_dege, so, alpha, beta, base_r));
  
  Tbl valeurs (val_solp (res, alpha, base_r)) ;
  sp_plus = valeurs(0) ;
  sp_minus = valeurs(1) ;
  dsp_plus = valeurs(2) ;
  dsp_minus = valeurs(3) ;
  
  return res ;
}

Tbl Ope_helmholtz_minus::get_solh() const {
  
  double rlim, rminus, rplus;  

  switch (base_r) {
  case R_CHEBU:
    // SH est exp(-beta*r)/r
    rlim = -0.5 / alpha ;
    s_one_minus = exp(-masse*rlim)/rlim/sqrt(2) ;
    ds_one_minus = -s_one_minus * (masse+1./rlim) ;
    break ;

  case R_CHEB:
    // SH_one est sin(masse*r)/r :
    rminus = beta-alpha ;
    rplus = beta+alpha ;
    s_one_minus = exp(masse*rminus)/rminus/sqrt(2) ;
    ds_one_minus = exp(masse*rminus)*(masse/rminus - 1./rminus/rminus)/sqrt(2) ;
    s_one_plus = exp(masse*rplus)/rplus/sqrt(2) ;
    ds_one_plus = exp(masse*rplus)*(masse/rplus - 1./rplus/rplus)/sqrt(2) ;
    
    // Sh two est cos(masse*r)/r :
    s_two_minus = exp(-masse*rminus)/rminus/sqrt(2) ;
    ds_two_minus = exp(-masse*rminus)*(-masse/rminus - 1./rminus/rminus)/
      sqrt(2) ;
    s_two_plus =  exp(-masse*rplus)/rplus/sqrt(2) ;
    ds_two_plus = exp(-masse*rplus)*(-masse/rplus - 1./rplus/rplus)/sqrt(2) ;
    break ;
  default:
    cout << "Case unkwnown in Ope_helmholtz_minus::get_solh" << endl ;
    abort() ;
    break ;
  }

  return solh_helmholtz_minus (nr, alpha, beta, masse, base_r) ;
}


void Ope_helmholtz_minus::inc_l_quant() {

  cout << "inc_l_quant not implemented for Helmholtz operator." << endl ;
  abort() ;
}
