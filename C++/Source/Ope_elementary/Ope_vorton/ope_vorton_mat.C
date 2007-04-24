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

char ope_vorton_mat_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <math.h>
#include <stdlib.h>

#include "proto.h"
#include "ope_elementary.h"
#include "diff.h"


                      //-----------------------------------
                     // Routine pour les cas non prevus -- 
                     //-----------------------------------

Matrice _vorton_mat_pas_prevu(int, double, double, int, int) {
  cout << "Vorton : base not implemented..." << endl ;
  abort() ;
  exit(-1) ;
  Matrice res(1, 1) ;
  return res;
}
                    //-------------------------
		    //--   CAS R_CHEBU   -----
		    //------------------------

Matrice _vorton_mat_r_chebu_trois (int n, int lq) {


    Diff_xdsdx2 xdd(R_CHEBU, n) ;
    Diff_dsdx d(R_CHEBU, n) ;
    Diff_sx sx(R_CHEBU, n) ;

    return (xdd+2*d-lq*(lq+1)*sx) ;
}

Matrice _vorton_mat_r_chebu (int n, double, double, int lq, int dz) {
    Matrice res(n, n) ;
    res.set_etat_qcq() ;
    switch (dz) {
	case 3 :
	    res = _vorton_mat_r_chebu_trois (n, lq) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
    }
    return res ;
}



                    //-------------------------
		    //--   CAS R_CHEB   -----
		    //------------------------

Matrice _vorton_mat_r_cheb (int n, double alpha, double beta, int lq, int) {
       Diff_dsdx2 d2(R_CHEB, n) ;
       Diff_xdsdx2 xd2(R_CHEB, n) ;
       Diff_x2dsdx2 x2d2(R_CHEB, n) ;
       Diff_dsdx d1(R_CHEB, n) ;
       Diff_xdsdx xd1(R_CHEB, n) ;
       Diff_id id(R_CHEB, n) ;

       double echelle = beta/alpha ;
       return (x2d2 + (2*echelle)*xd2 + (echelle*echelle)*d2 - (lq*(lq+1))*id) ;
}

void Ope_vorton::do_ope_mat() const {
  if (ope_mat != 0x0) 
    delete ope_mat ;

  // Routines de derivation
  static Matrice (*vorton_mat[MAX_BASE])(int, double, double, int, int);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      vorton_mat[i] = _vorton_mat_pas_prevu ;
    }
    // Les routines existantes
    vorton_mat[R_CHEB >> TRA_R] = _vorton_mat_r_cheb ;
    vorton_mat[R_CHEBU >> TRA_R] = _vorton_mat_r_chebu ;
  }
  ope_mat = new Matrice(vorton_mat[base_r](nr, alpha, beta, l_quant, dzpuis)) ;
}
