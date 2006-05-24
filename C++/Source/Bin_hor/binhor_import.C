/*
 *  Methods of class Bin_hor
 *
 *   (see file bin_hor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose Luis Jaramillo
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


char binhor_import_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2006/05/24 17:00:29  f_limousin
 * New function
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "nbr_spx.h"
#include "tenseur.h"
#include "tensor.h"
#include "isol_hor.h"
#include "proto.h"
#include "utilitaires.h"
//#include "graphique.h"

void Bin_hor::import_bh (const Bin_hor& bin_sol){

  double jjtime = hole1.jtime ;
  double ttime = hole1.the_time[jjtime] ;

  // Importation of n_auto
  // ----------------------

  Scalar temp1 (hole1.mp) ;
  temp1.import(bin_sol(1).n_auto()) ;
  temp1.std_spectral_base() ;
  hole1.n_auto_evol.update(temp1, jjtime, ttime) ;

  Scalar temp2 (hole2.mp) ;
  temp2.import(bin_sol(2).n_auto()) ;
  temp2.std_spectral_base() ;
  hole2.n_auto_evol.update(temp2, jjtime, ttime) ;
    
  // Importation of psi_auto
  // ----------------------

  temp1.import(bin_sol(1).psi_auto()) ;
  temp1.std_spectral_base() ;
  hole1.psi_auto_evol.update(temp1, jjtime, ttime) ;

  temp2.import(bin_sol(2).psi_auto()) ;
  temp2.std_spectral_base() ;
  hole2.psi_auto_evol.update(temp2, jjtime, ttime) ;

  // Importation of beta_auto
  // ----------------------

  Vector tmp_vect1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  Vector shift_sol1 (bin_sol(1).beta_auto()) ;
  shift_sol1.change_triad(bin_sol(1).mp.get_bvect_cart()) ;
  shift_sol1.change_triad(hole1.mp.get_bvect_cart()) ;
  assert (*(shift_sol1.get_triad()) == *(tmp_vect1.get_triad())) ;

  tmp_vect1.set(1).import(shift_sol1(1)) ;
  tmp_vect1.set(2).import(shift_sol1(2)) ;
  tmp_vect1.set(3).import(shift_sol1(3)) ;
  tmp_vect1.std_spectral_base() ;
  tmp_vect1.change_triad(hole1.mp.get_bvect_spher()) ;

  hole1.beta_auto_evol.update(tmp_vect1, jjtime, ttime) ;



  Vector tmp_vect2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  Vector shift_sol2 (bin_sol(2).beta_auto()) ;
  shift_sol2.change_triad(bin_sol(2).mp.get_bvect_cart()) ;
  shift_sol2.change_triad(hole2.mp.get_bvect_cart()) ;
  assert (*(shift_sol2.get_triad()) == *(tmp_vect2.get_triad())) ;

  tmp_vect2.set(1).import(shift_sol2(1)) ;
  tmp_vect2.set(2).import(shift_sol2(2)) ;
  tmp_vect2.set(3).import(shift_sol2(3)) ;
  tmp_vect2.std_spectral_base() ;
  tmp_vect2.change_triad(hole2.mp.get_bvect_spher()) ;

  hole2.beta_auto_evol.update(tmp_vect2, jjtime, ttime) ;


  // Initialisation of gamt, gamt_point, trK and trK_point
  // -------------------------------------------------------

  hole1.init_met_trK() ;
  hole2.init_met_trK() ;

  hole1.set_gamt(hole1.ff) ;
  hole2.set_gamt(hole2.ff) ;

  set_omega(bin_sol.get_omega()) ;

}

