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

char ope_poisson_pseudo_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/08/24 09:14:48  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * 
 * $Header$
 *
 */

#include "proto.h"
#include "ope_elementary.h"

// Standard constructor :
Ope_poisson_pseudo_1d::Ope_poisson_pseudo_1d (int nbr, int baser,
					      double alf, double bet, int lq): 
  Ope_elementary(nbr, baser, alf, bet), l_quant (lq) {
}

// Constructor by copy :
Ope_poisson_pseudo_1d::Ope_poisson_pseudo_1d (const Ope_poisson_pseudo_1d& so) : 
  Ope_elementary(so), 
  l_quant (so.l_quant) {
}

// Destructor :
Ope_poisson_pseudo_1d::~Ope_poisson_pseudo_1d() {} 

void Ope_poisson_pseudo_1d::inc_l_quant() {

  cout << "inc_l_quant not implemented for this operator." << endl ;
  abort() ;
}

void Ope_poisson_pseudo_1d::dec_l_quant() {

  cout << "dec_l_quant not implemented for this operator." << endl ;
  abort() ;
}
