/*
 *  Methods transverse( ) and longit_pot( ) of class Sym_tensor
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003  Eric Gourgoulhon & Jerome Novak
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

char sym_tensor_decomp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/11/26 21:57:03  e_gourgoulhon
 * First version; not ready yet.
 *
 *
 * $Header$
 *
 */


// C headers
#include <stdlib.h>

// Lorene headers
#include "tensor.h"

const Sym_tensor_trans Sym_tensor::transverse(const Metric& metre) const {

  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_transverse[j] == 0x0) {
   	cout << "Sym_tensor::transverse : not ready yet !" << endl ; 
	abort() ; 
  }

  return *p_transverse[j] ;
	

}




const Vector Sym_tensor::longit_pot(const Metric& metre) const {

  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_longit_pot[j] == 0x0) {
   	cout << "Sym_tensor::longit_pot : not ready yet !" << endl ; 
	abort() ; 
  }

  return *p_longit_pot[j] ;
	

}

