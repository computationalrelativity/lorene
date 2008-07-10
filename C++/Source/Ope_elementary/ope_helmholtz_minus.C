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
 * Revision 1.6  2008/07/10 10:59:17  p_grandclement
 * forgot another ones
 *
 * Revision 1.5  2004/08/24 09:14:45  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.4  2004/01/15 09:15:38  p_grandclement
 * Modification and addition of the Helmholtz operators
 *
 * Revision 1.3  2003/12/11 16:12:10  e_gourgoulhon
 * Changed sqrt(2) to sqrt(double(2)).
 *
 * Revision 1.2  2003/12/11 15:57:27  p_grandclement
 * include stdlib.h encore ...
 *
 * Revision 1.1  2003/12/11 14:48:50  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */
#include <math.h>
#include <stdlib.h>

#include "proto.h"
#include "ope_elementary.h"

// Standard constructor :
Ope_helmholtz_minus::Ope_helmholtz_minus (int nbr, int base, int lquant, 
					  double alf, double bet, double mas): 
  Ope_elementary(nbr, base, alf, bet), lq (lquant), masse (mas) {
}

// Constructor by copy :
Ope_helmholtz_minus::Ope_helmholtz_minus (const Ope_helmholtz_minus& so) : 
  Ope_elementary(so), lq(so.lq), masse (so.masse) {
}

// Destructor :
Ope_helmholtz_minus::~Ope_helmholtz_minus() {} 

// True functions :
void Ope_helmholtz_minus::do_ope_mat() const {
  if (ope_mat != 0x0)
    delete ope_mat ;

  ope_mat = new Matrice 
    (helmholtz_minus_mat(nr, lq, alpha, beta, masse, base_r)) ;
}

void Ope_helmholtz_minus::do_ope_cl() const {
  if (ope_mat == 0x0)
    do_ope_mat() ;

  if (ope_cl != 0x0)
    delete ope_cl ;

  ope_cl = new Matrice 
    (cl_helmholtz_minus(*ope_mat, base_r)) ;
}

void Ope_helmholtz_minus::do_non_dege() const {
  if (ope_cl == 0x0)
    do_ope_cl() ;

  if (non_dege != 0x0)
    delete non_dege ;

  non_dege = new Matrice 
    (prepa_helmholtz_minus_nondege(*ope_cl, base_r)) ;
}
  
Tbl Ope_helmholtz_minus::get_solp (const Tbl& so) const {

  if (non_dege == 0x0)
    do_non_dege() ;

  Tbl res(solp_helmholtz_minus (*ope_mat, *non_dege, so, alpha, beta, lq, base_r));
  
  Tbl valeurs (val_solp (res, alpha, base_r)) ;
  sp_plus = valeurs(0) ;
  sp_minus = valeurs(1) ;
  dsp_plus = valeurs(2) ;
  dsp_minus = valeurs(3) ;
  
  return res ;
}

Tbl Ope_helmholtz_minus::get_solh() const {
  
  Tbl res (solh_helmholtz_minus (nr, lq, alpha, beta, masse, base_r)) ;
   
  // Un peu tricky...
  if (res.get_ndim() == 1) {
    Tbl val_lim (val_solp (res, alpha, base_r)) ;

    s_one_plus   = val_lim(0) ;
    s_one_minus  = val_lim(1) ; 
    ds_one_plus  = val_lim(2) ;
    ds_one_minus = val_lim(3) ;

  }
  else {
    Tbl auxi (nr) ;
    auxi.set_etat_qcq() ;
    for (int i=0 ; i<nr ; i++)
      auxi.set(i) = res(0,i) ;

    Tbl val_one  (val_solp (auxi, alpha, base_r)) ; 
   
    s_one_plus   = val_one(0) ;
    s_one_minus  = val_one(1) ; 
    ds_one_plus  = val_one(2) ;
    ds_one_minus = val_one(3) ;

    for (int i=0 ; i<nr ; i++)
      auxi.set(i) = res(1,i) ;

    Tbl val_two  (val_solp (auxi, alpha, base_r)) ;

    s_two_plus   = val_two(0) ;
    s_two_minus  = val_two(1) ; 
    ds_two_plus  = val_two(2) ;
    ds_two_minus = val_two(3) ;   
  }
  return res ;
}



void Ope_helmholtz_minus::inc_l_quant() {

  cout << "inc_l_quant not implemented for Helmholtz operator." << endl ;
  abort() ;
}
