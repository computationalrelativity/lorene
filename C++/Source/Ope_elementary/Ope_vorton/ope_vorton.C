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

char ope_vorton_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "ope_elementary.h"

// Standard constructor :
Ope_vorton::Ope_vorton (int nbr, int base, double alf, 
				    double bet, int lq, int dz) : 
  Ope_elementary(nbr, base, alf, bet), l_quant(lq), dzpuis(dz) {
}

// Constructor by copy :
Ope_vorton::Ope_vorton (const Ope_vorton& so) : 
  Ope_elementary(so), l_quant(so.l_quant), dzpuis(so.dzpuis) {
}

// Destructor :
Ope_vorton::~Ope_vorton() {} 

void Ope_vorton::inc_l_quant() {

  cout << "inc_l_quant not implemented for this operator." << endl ;
  abort() ;
}

void Ope_vorton::dec_l_quant() {

  cout << "dec_l_quant not implemented for this operator." << endl ;
  abort() ;
}
