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

char ope_helmholtz_plus_C[] = "$Header$" ;

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
Ope_helmholtz_plus::Ope_helmholtz_plus (int nbr, int base, double alf, 
					  double bet, double mas): 
  Ope_elementary(nbr, base, alf, bet), masse (mas) {
}

// Constructor by copy :
Ope_helmholtz_plus::Ope_helmholtz_plus (const Ope_helmholtz_plus& so) : 
  Ope_elementary(so), masse (so.masse) {
}

// Destructor :
Ope_helmholtz_plus::~Ope_helmholtz_plus() {} 

// True functions :
void Ope_helmholtz_plus::do_ope_mat() const {
  if (ope_mat != 0x0)
    delete ope_mat ;

  ope_mat = new Matrice 
    (helmholtz_plus_mat(nr, alpha, beta, masse, base_r)) ;
}

void Ope_helmholtz_plus::do_ope_cl() const {
  if (ope_mat == 0x0)
    do_ope_mat() ;

  if (ope_cl != 0x0)
    delete ope_cl ;

  ope_cl = new Matrice 
    (cl_helmholtz_plus(*ope_mat, alpha, beta, masse, base_r)) ;
}

void Ope_helmholtz_plus::do_non_dege() const {
  if (ope_cl == 0x0)
    do_ope_cl() ;

  if (non_dege != 0x0)
    delete non_dege ;

  non_dege = new Matrice 
    (prepa_helmholtz_plus_nondege(*ope_cl, alpha, beta, masse, base_r)) ;
}
  
Tbl Ope_helmholtz_plus::get_solp (const Tbl& so) const {

  if (non_dege == 0x0)
    do_non_dege() ;

  Tbl res(solp_helmholtz_plus (*ope_mat, *non_dege, so, alpha, beta, base_r));
  
  Tbl valeurs (val_solp (res, alpha, base_r)) ;
  sp_plus = valeurs(0) ;
  sp_minus = valeurs(1) ;
  dsp_plus = valeurs(2) ;
  dsp_minus = valeurs(3) ;
  
  return res ;
}

Tbl Ope_helmholtz_plus::get_solh() const {

  // SH_one est sin(masse*r)/r :
  double rminus = beta - alpha ;
  double rplus = beta + alpha ;
  
  s_one_minus = sin(masse*rminus)/rminus/sqrt(2) ;
  ds_one_minus = (masse*cos(masse*rminus)-sin(masse*rminus)/rminus)/
    rminus/sqrt(2) ;
  s_one_plus = sin(masse*rplus)/rplus/sqrt(2) ;
  ds_one_plus = (masse*cos(masse*rplus)-sin(masse*rplus)/rplus)/
    rplus/sqrt(2) ;
  
  // Sh two est cos(masse*r)/r :
  s_two_minus = cos(masse*rminus)/rminus/sqrt(2) ;
  ds_two_minus = (-masse*sin(masse*rminus)-cos(masse*rminus)/rminus)/
    rminus/sqrt(2) ;
  s_two_plus = cos(masse*rplus)/rplus/sqrt(2) ;
  ds_two_plus = (-masse*sin(masse*rplus)-cos(masse*rplus)/rplus)/
    rplus/sqrt(2) ;
  
  
  return solh_helmholtz_plus (nr, alpha, beta, masse, base_r) ;
}


void Ope_helmholtz_plus::inc_l_quant() {

  cout << "inc_l_quant not implemented for Helmholtz operator." << endl ;
  abort() ;
}
