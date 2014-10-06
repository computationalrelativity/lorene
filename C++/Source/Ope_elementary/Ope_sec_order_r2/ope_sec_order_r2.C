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

char ope_sec_order_r2_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "ope_elementary.h"

// Standard constructor :
Ope_sec_order_r2::Ope_sec_order_r2 (int nbr, int base, double alf, 
				    double bet, double a, double b, double c) : 
  Ope_elementary(nbr, base, alf, bet), a_param(a), b_param(b), c_param(c) {

  assert (a!=0) ;
}

// Constructor by copy :
Ope_sec_order_r2::Ope_sec_order_r2 (const Ope_sec_order_r2& so) : 
  Ope_elementary(so), a_param (so.a_param), b_param(so.b_param), c_param(so.c_param) {
}

// Destructor :
Ope_sec_order_r2::~Ope_sec_order_r2() {} 

void Ope_sec_order_r2::inc_l_quant() {

  cout << "inc_l_quant not implemented for Helmholtz operator." << endl ;
  abort() ;
}
