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

char ope_poisson_C[] = "$Header$" ;

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

#include "proto.h"
#include "ope_elementary.h"

// Standard constructor :
Ope_poisson::Ope_poisson (int nbr, int baser, double alf, double bet, int lq, int dz): 
  Ope_elementary(nbr, baser, alf, bet), l_quant (lq), 
  dzpuis (dz) {

  assert ((dzpuis==2) || (dzpuis==3) || (dzpuis==4)) ;
}

// Constructor by copy :
Ope_poisson::Ope_poisson (const Ope_poisson& so) : Ope_elementary(so), 
  l_quant (so.l_quant), dzpuis (so.dzpuis) {

  assert ((dzpuis==2) || (dzpuis==3) || (dzpuis==4)) ;
}

// Destructor :
Ope_poisson::~Ope_poisson() {} 

// True functions :
void Ope_poisson::do_ope_mat() const {
  if (ope_mat != 0x0)
    delete ope_mat ;

  ope_mat = new Matrice 
    (laplacien_mat(nr, l_quant, beta/alpha, dzpuis, base_r)) ;
}

void Ope_poisson::do_ope_cl() const {
  if (ope_mat == 0x0)
    do_ope_mat() ;

  if (ope_cl != 0x0)
    delete ope_cl ;

  ope_cl = new Matrice 
    (combinaison(*ope_mat, l_quant, beta/alpha, dzpuis, base_r)) ;
}

void Ope_poisson::do_non_dege() const {
  if (ope_cl == 0x0)
    do_ope_cl() ;

  if (non_dege != 0x0)
    delete non_dege ;

  non_dege = new Matrice 
    (prepa_nondege(*ope_cl, l_quant, beta/alpha, dzpuis, base_r)) ;
}
  
Tbl Ope_poisson::get_solp (const Tbl& so) const {

  if (non_dege == 0x0)
    do_non_dege() ;

  Tbl res(solp(*ope_mat, *non_dege, alpha, beta, so, dzpuis, base_r)) ;
  
  Tbl valeurs (val_solp (res, alpha, base_r)) ;
  sp_plus = valeurs(0) ;
  sp_minus = valeurs(1) ;
  dsp_plus = valeurs(2) ;
  dsp_minus = valeurs(3) ;

  return res ;
}

Tbl Ope_poisson::get_solh() const {
 
  Tbl valeurs (val_solh (l_quant, alpha, beta, base_r)) ;
  if (valeurs.get_ndim() == 2) {
    // cas 2 sh
    s_one_plus = valeurs(0,0) ;
    s_one_minus = valeurs(0,1) ;
    ds_one_plus = valeurs(0,2) ;
    ds_one_minus = valeurs(0,3) ;

    s_two_plus = valeurs(1,0) ;
    s_two_minus = valeurs(1,1) ;
    ds_two_plus = valeurs(1,2) ;
    ds_two_minus = valeurs(1,3) ;
  }
  else {
    // cas 1 sh :
    s_one_plus = valeurs(0) ;
    s_one_minus = valeurs(1) ;
    ds_one_plus = valeurs(2) ;
    ds_one_minus = valeurs(3) ;
  } 
 
  return solh(nr, l_quant, beta/alpha, base_r) ;
}

void Ope_poisson::inc_l_quant() {
  
  if (ope_mat != 0x0) {
    delete ope_mat ;
    ope_mat = 0x0 ;
  }
  
  if (ope_cl != 0x0) {
    delete ope_cl ;
    ope_cl = 0x0 ;
  }

  if (non_dege != 0x0) {
    delete non_dege ;
    non_dege = 0x0 ;
  } 
  l_quant ++ ;
}
