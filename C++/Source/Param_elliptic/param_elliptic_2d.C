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

char param_elliptic_2d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/08/24 09:14:49  p_grandclement
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

#include "headcpp.h"

#include <math.h>
#include <stdlib.h>

#include "base_val.h" 
#include "map.h"
#include "ope_elementary.h"
#include "param_elliptic.h"
#include "change_var.h"
#include "scalar.h"


void Param_elliptic::set_poisson_2d(const Scalar& source, bool indic) {

  int dzpuis = source.get_dzpuis() ;

  if (type_map != MAP_AFF) {
    cout << "set_poisson_2d only defined for an affine mapping..." << endl ;
    abort() ;
  }
  else {

    int nz = get_mp().get_mg()->get_nzone() ;
    
    int nr ;
    double alpha, beta ;
    int m_quant, l_quant, base_r_1d ;

    int conte = 0 ;
    for (int l=0 ; l<nz ; l++) {
      
      nr = get_mp().get_mg()->get_nr(l) ;
      alpha = get_alpha (l) ;
      beta = get_beta (l) ;

      for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	  if (operateurs[conte] != 0x0)	    
	    delete operateurs[conte] ;
	    source.get_spectral_va().base.give_quant_numbers(l, k, j, m_quant, l_quant, base_r_1d) ;
	    if (k!=1)
	      if ((indic) || ((!indic) && (l_quant !=0)))
		operateurs[conte] = new Ope_poisson_2d (nr, base_r_1d, alpha, beta, l_quant, dzpuis) ;
	      else
		operateurs[conte] = 0x0 ;
	    else
	      operateurs[conte] = 0x0 ;
	    conte ++ ;
	  }
    }
  }
}

void Param_elliptic::set_helmholtz_minus_2d(int zone, double masse, const Scalar& source) {

  int dzpuis = source.get_dzpuis() ;
  assert (masse > 0) ;

  if (type_map != MAP_AFF) {
    cout << "set_helmholtz_minus_2d only defined for an affine mapping..." << endl ;
    abort() ;
  }
  else {

    int nz = get_mp().get_mg()->get_nzone() ;
    if (zone == nz-1)  
      source.check_dzpuis(2) ; 
    int nr ;
    double alpha, beta ;
    int m_quant, l_quant, base_r_1d ;

    int conte = 0 ;
    for (int l=0 ; l<nz ; l++) {
      
      nr = get_mp().get_mg()->get_nr(l) ;
      alpha = get_alpha (l) ;
      beta = get_beta (l) ;

      for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	  if (l==zone) {
	    if (operateurs[conte] != 0x0) {	    
	      delete operateurs[conte] ;
	      source.get_spectral_va().base.give_quant_numbers 
		(l, k, j, m_quant, l_quant, base_r_1d) ;
	      operateurs[conte] = new Ope_helmholtz_minus_2d (nr, base_r_1d, alpha, beta, l_quant, masse, dzpuis) ;
	    }
	  } 
	  conte ++ ;
	}
    }
  }
}

