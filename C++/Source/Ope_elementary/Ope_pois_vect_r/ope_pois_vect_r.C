/*
 *  Methods of class Ope_pois_vect_r
 *
 *   (see file ope_elementary.h for documentation)
 *
 */

/*
 *   Copyright (c) 2004 Jerome Novak
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

char ope_pois_vect_r_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/05/10 15:28:22  j_novak
 * First version of functions for the solution of the r-component of the
 * vector Poisson equation.
 *
 * 
 * $Header$
 *
 */

#include"type_parite.h"
#include "ope_elementary.h"

Matrice ope_pvect_r_mat(int , int , double , int, int ) ;
Tbl sh_pvect_r(int, int, double, int) ;
Matrice cl_pvect_r (const Matrice&, int, double, int, int) ;
Matrice nondeg_pvect_r (const Matrice&, int, double, int, int) ;
Tbl val_solh (int, double, double, int) ;

// Standard constructor :
Ope_pois_vect_r::Ope_pois_vect_r (int nbr, int baser, double alf, double bet, int lq, int dz): 
  Ope_poisson(nbr, baser, alf, bet, lq, dz) 
{
  assert (dzpuis == 4) ;
}

// Constructor by copy :
Ope_pois_vect_r::Ope_pois_vect_r (const Ope_pois_vect_r& so) : Ope_poisson(so)
{
  assert (dzpuis == 4) ;

}

// Destructor :
Ope_pois_vect_r::~Ope_pois_vect_r() {} 

// True functions :
void Ope_pois_vect_r::do_ope_mat() const {
  if (ope_mat != 0x0)
    delete ope_mat ;

  ope_mat = new Matrice 
    (ope_pvect_r_mat(nr, l_quant, beta/alpha, dzpuis, base_r)) ;
}

void Ope_pois_vect_r::do_ope_cl() const {
  if (ope_mat == 0x0)
    do_ope_mat() ;

  if (ope_cl != 0x0)
    delete ope_cl ;

  ope_cl = new Matrice 
    (cl_pvect_r(*ope_mat, l_quant, beta/alpha, dzpuis, base_r)) ;
}

void Ope_pois_vect_r::do_non_dege() const {
  if (ope_cl == 0x0)
    do_ope_cl() ;

  if (non_dege != 0x0)
    delete non_dege ;

  non_dege = new Matrice 
    (nondeg_pvect_r(*ope_cl, l_quant, beta/alpha, dzpuis, base_r)) ;
}
  
Tbl Ope_pois_vect_r::get_solh() const {

  int l1 = ( (l_quant == 0) ? 1 : l_quant - 1 ) ;
  int l2 = l_quant + 1 ;
 
  Tbl valeurs1 (val_solh (l1, alpha, beta, base_r)) ;
  Tbl valeurs2 (val_solh (l2, alpha, beta, base_r)) ;

  assert (valeurs1.get_ndim() == valeurs2.get_ndim()) ;

  if (valeurs1.get_ndim() == 2) {
    // cas 2 sh
    s_one_plus = valeurs1(0,0) ;
    s_one_minus = valeurs1(0,1) ;
    ds_one_plus = valeurs1(0,2) ;
    ds_one_minus = valeurs1(0,3) ;

    s_two_plus = valeurs2(1,0) ;
    s_two_minus = valeurs2(1,1) ;
    ds_two_plus = valeurs2(1,2) ;
    ds_two_minus = valeurs2(1,3) ;
  }
  else {
    // cas 1 sh :
    Tbl&  valeurs = (base_r == R_CHEBU) ? valeurs2 : valeurs1 ;
    s_one_plus = valeurs(0) ;
    s_one_minus = valeurs(1) ;
    ds_one_plus = valeurs(2) ;
    ds_one_minus = valeurs(3) ;
  } 
 
  return sh_pvect_r(nr, l_quant, beta/alpha, base_r) ;
}

