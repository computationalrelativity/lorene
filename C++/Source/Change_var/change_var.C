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

char change_var_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/12/11 15:53:31  p_grandclement
 * includ stdlib
 *
 * Revision 1.1  2003/12/11 14:48:48  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */

#include <math.h>
#include <stdlib.h>
#include <iostream.h>

#include "proto.h"
#include "change_var.h"

// Les fonctions elementaires dont on a besoin
double one (double) {
  return 1 ;
}

double zero (double) {
  return 0 ;
}

double ide (double x) {
  return x ;
}

double part_ln (double x) {
  return 1+x*x*log(x)/3. ;
}

double part_ln_der (double x) {
  return 2./3.*x*log(x)+x/3. ;
}

// Construction du changement de variable ...
Change_var::Change_var (int type_change) {
  
  switch (type_change) {
  case STD:
    func_F = zero ;
    der_F = zero ;
    func_G = one ;
    der_G = zero ;
    break ;
  
  case W_BETA:
    func_F = one ;
    der_F = zero ;
    func_G = ide ;
    der_G = one ;
    break ;

  case W_BETA_INF:
    func_F = part_ln ;
    der_F = part_ln_der ;
    func_G = ide ;
    der_G = one ;
    break ;

  case H_BETA:
    func_F = one ;
    der_F = zero ;
    func_G = one ;
    der_G = zero ;
    break ;

  default:
    cout << "Unknown type in Change_var::Change_var(int)" << endl ;
    abort() ;
    break ;
  }
}

Change_var::Change_var (const Change_var& so) :  
  func_F(so.func_F), der_F(so.der_F), func_G(so.func_G), der_G(so.der_G) {}

Change_var::~Change_var() {}

double Change_var::val_F (double air) {
  return (*func_F)(air) ;
}

double Change_var::val_der_F (double air) {
  return (*der_F)(air) ;
}

double Change_var::val_G (double air) {
  return (*func_G)(air) ;
}

double Change_var::val_der_G (double air) {
  return (*der_G)(air) ;
}

