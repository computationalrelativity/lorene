/*
 *  File to compare two tensors
 *
 */

/*
 *   Copyright (c) 2003 Francois Limousin
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

char tenseur_compare_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/05/27 07:17:19  p_grandclement
 * Correction of some shadowed variables
 *
 * Revision 1.1  2003/06/20 15:01:49  f_limousin
 * *** empty log message ***
 *
 * Revision 1.2  2001/12/11 06:44:41  e_gourgoulhon
 * template files
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include "math.h"

// Lorene headers
#include "tenseur.h"


void Tenseur::compare(const Tenseur& tens, const char* name) {
  assert ( valence == tens.get_valence() ) ;

  cout << "----------------------------------------------" << endl ;
  cout << "Comparison of " << name << " : " << endl << endl ;

  if (valence == 0) {
    Cmp comp = tens() ;
    Cmp tmp = (*this)() ;
    tmp.compare(comp, name) ;
  }
    
  if (valence == 1) {
    Cmp comp0 = tens(0) ;
    Cmp comp1 = tens(1) ;
    Cmp comp2 = tens(2) ;

    Cmp tmp0 = (*this)(0) ;
    Cmp tmp1 = (*this)(1) ;
    Cmp tmp2 = (*this)(2) ;

    tmp0.compare(comp0, name, 0) ;
    tmp1.compare(comp1, name, 1) ;
    tmp2.compare(comp2, name, 2) ;
  }

   if (valence == 2) {
    Cmp comp00 = tens(0,0) ;
    Cmp comp01 = tens(0,1) ;
    Cmp comp02 = tens(0,2) ;
    Cmp comp10 = tens(1,0) ;
    Cmp comp11 = tens(1,1) ;
    Cmp comp12 = tens(1,2) ;
    Cmp comp20 = tens(2,0) ;
    Cmp comp21 = tens(2,1) ;
    Cmp comp22 = tens(2,2) ;

    Cmp tmp00 = (*this)(0,0) ;
    Cmp tmp01 = (*this)(0,1) ;
    Cmp tmp02 = (*this)(0,2) ;
    Cmp tmp10 = (*this)(1,0) ;
    Cmp tmp11 = (*this)(1,1) ;
    Cmp tmp12 = (*this)(1,2) ;
    Cmp tmp20 = (*this)(2,0) ;
    Cmp tmp21 = (*this)(2,1) ;
    Cmp tmp22 = (*this)(2,2) ;

    tmp00.compare(comp00, name, 0, 0) ;
    tmp01.compare(comp01, name, 0, 1) ;
    tmp02.compare(comp02, name, 0, 2) ;
    tmp10.compare(comp10, name, 1, 0) ;
    tmp11.compare(comp11, name, 1, 1) ;
    tmp12.compare(comp12, name, 1, 2) ;
    tmp20.compare(comp20, name, 2, 0) ;
    tmp21.compare(comp21, name, 2, 1) ;
    tmp22.compare(comp22, name, 2, 2) ;
}

   if (valence > 2 ) {
     abort() ;
   }
}

void Tenseur::compare(FILE* fich, const char* name_i) {

  Mg3d mg(fich) ;
  Map_et mpg(mg, fich) ;

  Tenseur tens(mpg, mpg.get_bvect_cart(), fich) ;

  compare(tens, name_i) ;


}
