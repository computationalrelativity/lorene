/*
 *   Copyright (c) 2004 Philippe Grandclement
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


char base_val_quantum_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2014/10/06 15:12:57  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2013/01/11 08:20:11  j_novak
 * New radial spectral bases with Legendre polynomials (R_LEG, R_LEGP, R_LEGI).
 *
 * Revision 1.7  2009/10/23 12:55:16  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.6  2009/10/08 16:20:13  j_novak
 * Addition of new bases T_COS and T_SIN.
 *
 * Revision 1.5  2007/12/11 15:28:09  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.4  2005/09/07 13:09:50  j_novak
 * New method for determining the highest multipole that can be described on a 3D
 * grid.
 *
 * Revision 1.3  2004/11/23 15:08:01  m_forot
 * Added the bases for the cases without any equatorial symmetry
 * (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.2  2004/08/30 16:27:59  r_prix
 * added #include <stdlib.h> (got ERROR 'abort' is undefined without this...)
 *
 * Revision 1.1  2004/08/24 09:14:41  p_grandclement
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

// Headers C
#include <cstdio>
#include <cassert>
#include <cstdlib>

// Headers Lorene
#include "grilles.h"
#include "base_val.h"
#include "utilitaires.h"

void Base_val::give_quant_numbers (int l, int k, int j, 
			 int& m_quant, int& l_quant, int& base_r_1d) const {

  int base_p = get_base_p(l) ;
  int base_t = get_base_t(l) ;
  int base_r = get_base_r(l) ;
  
  switch (base_p) {
  case P_COSSIN :
      m_quant = k/2 ;
    break;
    
  case P_COSSIN_P :
    if (k%2 == 0)
      m_quant = k ;
    else
      m_quant = (k-1) ;
    break;

  case P_COSSIN_I :
      m_quant = 2*( (k-1) / 2) + 1 ;
    break;
  default:
    cout << "Unknown basis in phi in give_quant_numbers ..." << endl ;
    abort() ;
    break ;
  }

  switch (base_t) {
  case T_COS_P :
    l_quant = 2*j ;
    break;
    
  case T_SIN_P :
    l_quant = 2*j ;
    break;

  case T_COS_I :
    l_quant = 2*j+1 ;
    break;

  case T_SIN_I :
    l_quant = 2*j+1 ;
    break;
 
  case T_COSSIN_CP :
    if (m_quant%2 == 0)
      l_quant = 2*j ;
    else
      l_quant = 2*j+1 ;
    break ;
    
  case T_COSSIN_SP :
    if (m_quant%2 == 0)
      l_quant = 2*j ;
    else
      l_quant = 2*j+1 ;
    break ;
    
  case T_COSSIN_CI :
    if (m_quant%2 == 0)
      l_quant = 2*j+1 ;
    else
      l_quant = 2*j ;
    break ;
    
  case T_COSSIN_SI :
    if (m_quant%2 == 0)
      l_quant = 2*j+1 ;
    else
      l_quant = 2*j ;
    break ;

  case T_COSSIN_C :
       l_quant = j ;
    break ;
    
  case T_COSSIN_S :
       l_quant = j ;
    break ;   

  case T_COS :
       l_quant = j ;
    break ;
    
  case T_LEG_P :
    if (m_quant%2 == 0)
      l_quant = 2*j ;
    else
      l_quant = 2*j+1 ;
    break ;

 case T_LEG_PP :
   l_quant = 2*j ;
   break ;

  case T_LEG_I :
    if (m_quant%2 == 0)
      l_quant = 2*j+1 ;
    else
      l_quant = 2*j ;
    break ;
   
  case T_LEG_IP :
    l_quant = 2*j+1 ;
    break ;

  case T_LEG_PI :
    l_quant = 2*j+1 ;
    break ; 

  case T_LEG_II :
    l_quant = 2*j ;
    break ; 

  case T_LEG :
   l_quant = j ;
   break ;

  case T_LEG_MP :
   l_quant = j ;
   break ;

  case T_LEG_MI :
   l_quant = j ;
   break ;

  case T_CL_COS_P:
    l_quant = 2*j ;
    break ;
  
  case T_CL_SIN_P:
    l_quant = 2*j ;
    break ;
   
  case T_CL_COS_I:
    l_quant = 2*j+1 ;
    break ;
    
  case T_CL_SIN_I:
    l_quant = 2*j+1 ;
    break ;
    
  default:
    cout << "Unknown basis in theta in give_quant_numbers ..." << endl ;
    abort() ;
    break ;
  }

  switch (base_r) {
  case R_CHEB :
    base_r_1d = R_CHEB ;
    break ;

  case R_CHEBP :
    base_r_1d = R_CHEBP ;
    break ;
    
   case R_CHEBI :
    base_r_1d = R_CHEBI ;
    break ;
    
  case R_LEG :
    base_r_1d = R_LEG ;
    break ;

  case R_LEGP :
    base_r_1d = R_LEGP ;
    break ;
    
   case R_LEGI :
    base_r_1d = R_LEGI ;
    break ;
    
  case R_JACO02 :
    base_r_1d = R_JACO02 ;
    break ;
    
  case R_CHEBPIM_P :    
    if (m_quant%2 == 0) 
      base_r_1d = R_CHEBP ;
    else
      base_r_1d = R_CHEBI ;
    break ;
  case R_CHEBPIM_I :
    if (m_quant%2 == 0)
      base_r_1d = R_CHEBI ;
    else
      base_r_1d = R_CHEBP ;
    break ;
  case R_CHEBPI_P :    
    if (l_quant%2 == 0) 
      base_r_1d = R_CHEBP ;
    else
      base_r_1d = R_CHEBI ;
    break ;
  case R_CHEBPI_I :
    if (l_quant%2 == 0)
      base_r_1d = R_CHEBI ;
    else
      base_r_1d = R_CHEBP ;
    break ;   
  case R_CHEBU :
    base_r_1d = R_CHEBU ;
    break ;
    
  default:
    cout << "Unknown basis in r in give_quant_numbers ..." << endl ;
    abort() ;
    break ;
  }
}

int Base_val::give_lmax(const Mg3d& mgrid, int lz) const {

#ifndef NDEBUG
    int nz = mgrid.get_nzone() ;
    assert (lz < nz) ;
#endif

    int ntm1 = mgrid.get_nt(lz) - 1;
    int base_t = get_base_t(lz) ;
    bool m_odd = (mgrid.get_np(lz) > 2) ;

    int l_max = 0 ;

    switch (base_t) {
	case T_COS_P :
	    l_max = 2*ntm1 ;
	    break;
    
	case T_SIN_P :
	    l_max = 2*ntm1 ;
	    break;

	case T_COS_I :
	    l_max = 2*ntm1+1 ;
	    break;

	case T_SIN_I :
	    l_max = 2*ntm1+1 ;
	    break;
 
	case T_COSSIN_CP :
	    if (!m_odd)
		l_max = 2*ntm1 ;
	    else
		l_max = 2*ntm1+1 ;
	    break ;
    
	case T_COSSIN_SP :
	    if (!m_odd)
		l_max = 2*ntm1 ;
	    else
		l_max = 2*ntm1+1 ;
	    break ;
	    
	case T_COSSIN_CI :
	    if (!m_odd)
		l_max = 2*ntm1+1 ;
	    else
		l_max = 2*ntm1 ;
	    break ;
	    
	case T_COSSIN_SI :
	    if (!m_odd)
		l_max = 2*ntm1+1 ;
	    else
		l_max = 2*ntm1 ;
	    break ;
	    
	case T_COSSIN_C :
	    l_max = ntm1 ;
	    break ;
	    
	case T_COSSIN_S :
	    l_max = ntm1 ;
	    break ;   
	    
	case T_COS :
	    l_max = ntm1 ;
	    break ;   
	    
	case T_LEG_P :
	    if (!m_odd)
		l_max = 2*ntm1 ;
	    else
		l_max = 2*ntm1+1 ;
	    break ;
	    
	case T_LEG_PP :
	    l_max = 2*ntm1 ;
	    break ;
	    
	case T_LEG_I :
	    if (!m_odd)
		l_max = 2*ntm1+1 ;
	    else
		l_max = 2*ntm1 ;
	    break ;
	    
	case T_LEG_IP :
	    l_max = 2*ntm1+1 ;
	    break ;
	    
	case T_LEG_PI :
	    l_max = 2*ntm1+1 ;
	    break ; 
	    
	case T_LEG_II :
	    l_max = 2*ntm1 ;
	    break ; 
	    
	case T_LEG :
	    l_max = ntm1 ;
	    break ;
	    
	case T_LEG_MP :
	    l_max = ntm1 ;
	    break ;
	    
	case T_LEG_MI :
	    l_max = ntm1 ;
	    break ;
	    
	case T_CL_COS_P:
	    l_max = 2*ntm1 ;
	    break ;
	    
	case T_CL_SIN_P:
	    l_max = 2*ntm1 ;
	    break ;
	    
	case T_CL_COS_I:
	    l_max = 2*ntm1+1 ;
	    break ;
	    
	case T_CL_SIN_I:
	    l_max = 2*ntm1+1 ;
	    break ;
	    
	default:
	    cout << "Unknown basis in theta in Base_val::get_lmax ..." 
		 << endl ;
	    abort() ;
	    break ;
    }
    return l_max ;
}
