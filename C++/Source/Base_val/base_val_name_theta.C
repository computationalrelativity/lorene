/*
 *  Method Base_val::name_theta
 *
 *	(see file base_val.h for documentation). 
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon. 
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.11  2023/06/28 10:04:32  j_novak
 * Use of C++ strings and flows instead of C types.
 *
 * Revision 1.10  2016/12/05 16:17:44  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.9  2014/10/13 08:52:39  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:12:57  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2009/10/23 12:55:16  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.6  2009/10/08 16:20:13  j_novak
 * Addition of new bases T_COS and T_SIN.
 *
 * Revision 1.5  2004/12/17 13:35:01  m_forot
 * Add the case T_LEG
 *
 * Revision 1.4  2004/11/23 15:08:01  m_forot
 * Added the bases for the cases without any equatorial symmetry
 * (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.3  2004/10/04 13:40:38  j_novak
 * Added the T_COS base case.
 *
 * Revision 1.2  2004/08/24 09:14:41  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.1  2003/10/19 19:49:40  e_gourgoulhon
 * First version
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstring>
#include <cstdlib>

// Header C++
#include <sstream>

// Lorene headers
#include "type_parite.h"
#include "base_val.h"

// Local prototypes
namespace Lorene {
  void basename_t_unknown(int, int, string&) ; 
  void basename_t_cos(int, int, string&) ; 
  void basename_t_sin(int, int, string&) ; 
  void basename_t_cos_p(int, int, string&) ; 
  void basename_t_sin_p(int, int, string&) ; 
  void basename_t_cos_i(int, int, string&) ; 
  void basename_t_sin_i(int, int, string&) ; 
  void basename_t_cossin_cp(int, int, string&) ; 
  void basename_t_cossin_sp(int, int, string&) ;
  void basename_t_cossin_c(int, int, string&) ; 
  void basename_t_cossin_s(int, int, string&) ; 
  void basename_t_cossin_ci(int, int, string&) ; 
  void basename_t_cossin_si(int, int, string&) ; 
  void basename_t_leg_p(int, int, string&) ; 
  void basename_t_leg(int, int, string&) ; 
  void basename_t_leg_mp(int, int, string&) ; 
  void basename_t_leg_mi(int, int, string&) ; 
  void basename_t_leg_pp(int, int, string&) ; 
  void basename_t_leg_i(int, int, string&) ; 
  void basename_t_leg_ip(int, int, string&) ; 
  void basename_t_leg_pi(int, int, string&) ; 
  void basename_t_leg_ii(int, int, string&) ; 
  void basename_t_cl_cos_p(int, int, string&) ; 
  void basename_t_cl_sin_p(int, int, string&) ; 
  void basename_t_cl_cos_i(int, int, string&) ; 
  void basename_t_cl_sin_i(int, int, string&) ; 


  //----------------------------//
  //      Base_val method       //
  //----------------------------//

  void Base_val::name_theta(int l, int k, int j, string& name) const {

    // Array of actual base name functions
    static void(*vbasename_t[MAX_BASE])(int, int, string&) ;  

    static bool first_call = true ;

    // Initializations at first call
    // -----------------------------
    if ( first_call ) {

      first_call = false ;

      for (int i=0 ; i<MAX_BASE ; i++) {
	vbasename_t[i] = basename_t_unknown ;
      }

      vbasename_t[T_COS >> TRA_T] = basename_t_cos ;
      vbasename_t[T_SIN >> TRA_T] = basename_t_sin ;
      vbasename_t[T_COS_P >> TRA_T] = basename_t_cos_p ;
      vbasename_t[T_SIN_P >> TRA_T] = basename_t_sin_p ;
      vbasename_t[T_COS_I >> TRA_T] = basename_t_cos_i ;
      vbasename_t[T_SIN_I >> TRA_T] = basename_t_sin_i ;
      vbasename_t[T_COSSIN_CP >> TRA_T] = basename_t_cossin_cp ;
      vbasename_t[T_COSSIN_SP >> TRA_T] = basename_t_cossin_sp ;
      vbasename_t[T_COSSIN_CI >> TRA_T] = basename_t_cossin_ci ;
      vbasename_t[T_COSSIN_SI >> TRA_T] = basename_t_cossin_si ;
      vbasename_t[T_COSSIN_C >> TRA_T] = basename_t_cossin_c ;
      vbasename_t[T_COSSIN_S >> TRA_T] = basename_t_cossin_s ;
      vbasename_t[T_LEG_P >> TRA_T] = basename_t_leg_p ;
      vbasename_t[T_LEG_MP >> TRA_T] = basename_t_leg_mp ;
      vbasename_t[T_LEG_MI >> TRA_T] = basename_t_leg_mi ;
      vbasename_t[T_LEG >> TRA_T] = basename_t_leg ;
      vbasename_t[T_LEG_PP >> TRA_T] = basename_t_leg_pp ;
      vbasename_t[T_LEG_I >> TRA_T] = basename_t_leg_i ;
      vbasename_t[T_LEG_IP >> TRA_T] = basename_t_leg_ip ;
      vbasename_t[T_LEG_PI >> TRA_T] = basename_t_leg_pi ;
      vbasename_t[T_LEG_II >> TRA_T] = basename_t_leg_ii ;
      vbasename_t[T_CL_COS_P >> TRA_T] = basename_t_cl_cos_p ;
      vbasename_t[T_CL_SIN_P >> TRA_T] = basename_t_cl_sin_p ;
      vbasename_t[T_CL_COS_I >> TRA_T] = basename_t_cl_cos_i ;
      vbasename_t[T_CL_SIN_I >> TRA_T] = basename_t_cl_sin_i ;

    }
	
    // Call to the function adapted to the basis in domain l
    //------------------------------------------------------
	
    assert( (l>=0) && (l<nzone) ) ; 
	
    int base_t = ( b[l] & MSQ_T ) >> TRA_T ;
	
    vbasename_t[base_t](k, j, name) ; 

  }
	
	
  //-------------------------------//
  //  individual basis functions   //
  //-------------------------------//
	
  void basename_t_unknown(int, int, string&) {
    cout << "Base_val::name_theta : unknwon basis !" << endl ; 
    abort() ; 
  } 


  void basename_t_cos(int , int j, string& name) {

    assert( j>=0 ) ; 

    ostringstream ostr ;
    ostr << "cos" << j << 't' << flush ;
    name = ostr.str() ;
  }	

  void basename_t_sin(int , int j, string& name) {

    assert( j>=0 ) ; 

    ostringstream ostr ;
    ostr << "sin" << j << 't' << flush ;
    name = ostr.str() ;
  }	


  void basename_t_cos_p(int , int j, string& name) {

    assert( j>=0 ) ; 

    ostringstream ostr ;
    ostr << "cos" << 2*j << 't' << flush ;
    name = ostr.str() ;
  }	


  void basename_t_sin_p(int , int j, string& name) {

    assert( j>=0 ) ; 
  
    if (j == 0) {
      name = "unused" ; 
      return ;
    }
  
    ostringstream ostr ;
    ostr << "sin" << 2*j << 't' << flush ;
    name = ostr.str() ;
  }
	
  void basename_t_cl_cos_p(int , int j, string& name) {

    assert( j>=0 ) ; 

    ostringstream ostr ;
    ostr << "cl_cos" << 2*j << 't' << flush ;
    name = ostr.str() ;
  }

  void basename_t_cl_sin_p(int , int j, string& name) {

    assert( j>=0 ) ; 

    ostringstream ostr ;
    ostr << "cl_sin" << 2*j << 't' << flush ;
    name = ostr.str() ;
  }
	
  void basename_t_cos_i(int , int j, string& name) {

    assert( j>=0 ) ; 

    ostringstream ostr ;
    ostr << "cos" << (2*j + 1) << 't' << flush ;
    name = ostr.str() ;
  }	

  void basename_t_cl_cos_i(int , int j, string& name) {

    assert( j>=0 ) ; 

    ostringstream ostr ;
    ostr << "cl_cos" << (2*j + 1) << 't' << flush ;
    name = ostr.str() ;
  }

  void basename_t_sin_i(int , int j, string& name) {

    assert( j>=0 ) ; 

    ostringstream ostr ;
    ostr << "sin" << (2*j + 1) << 't' << flush ;
    name = ostr.str() ;
  }
	
  void basename_t_cl_sin_i(int , int j, string& name) {

    assert( j>=0 ) ; 

    ostringstream ostr ;
    ostr << "cl_sin" << (2*j + 1) << 't' << flush ;
    name = ostr.str() ;
  }	

	
  void basename_t_cossin_cp(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    ostringstream ostr ;
    int m = k / 2 ; 
    int xt ; 
    if (m%2 == 0) {
      ostr << "cos" ; 
      xt = 2*j ; 
    }
    else {
      ostr << "sin" ; 
      xt = 2*j + 1 ; 
    }
    ostr << xt << 't' << flush ;
    name = ostr.str() ;
  }	


  void basename_t_cossin_sp(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    ostringstream ostr ;
    int m = k / 2 ; 
    int xt ; 
    if (m%2 == 0) {
      if (j == 0) {
	name = "unused" ;
	return ;  
      }
      else {
	ostr << "sin" ; 
	xt = 2*j ;
      } 
    }
    else {
      ostr << "cos" ; 
      xt = 2*j + 1 ; 
    }
    ostr << xt << 't' << flush ;
    name = ostr.str() ;
  }	

  void basename_t_cossin_c(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    ostringstream ostr ;
    int m = k / 2 ; 
    int xt ; 
    if (m%2 == 0) {
      ostr << "cos" ; 
      xt = j ; 
    }
    else {
      if (j == 0) {
	name = "unused" ;
	return ;  
      } else {
	ostr << "sin" ; 
	xt = j ;
      } 
    }
    ostr << xt << 't' << flush ;
    name = ostr.str() ;
  }


  void basename_t_cossin_s(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    ostringstream ostr ;
    int m = k / 2 ; 
    int xt ; 
    if (m%2 == 0) {
      if (j == 0) {
	name = "unused" ;
	return ;  
      }
      else {
	ostr << "sin" ; 
	xt = j ;
      } 
    }
    else {
      ostr << "cos" ; 
      xt = j ; 
    }
    ostr << xt << 't' << flush ;
    name = ostr.str() ;
  }	


  void basename_t_cossin_ci(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    ostringstream ostr ;
    int m = k / 2 ; 
    int xt ; 
    if (m%2 == 0) {
      ostr << "cos" ; 
      xt = 2*j + 1; 
    }
    else {
      if (j == 0) {
	name = "unused" ;
	return ;  
      }
      else {
	ostr << "sin" ; 
	xt = 2*j ; 
      }
    }
    ostr << xt << 't' << flush ;
    name = ostr.str() ;
  }	


  void basename_t_cossin_si(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    ostringstream ostr ;
    int m = k / 2 ; 
    int xt ; 
    if (m%2 == 0) {
      ostr << "sin" ; 
      xt = 2*j + 1; 
    }
    else {
      ostr << "cos" ; 
      xt = 2*j ; 
    }
    ostr << xt << 't' << flush ;
    name = ostr.str() ;
  }	

  void basename_t_leg(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    int m = k / 2 ;
	 
    if (j < m/2) {
      name = "unused" ; 
      return ; 
    }

    ostringstream ostr ;
    ostr << "P_" << j << '^' << m << flush ;
    name = ostr.str() ;
  }	

  void basename_t_leg_mp(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    int m = 2 * (k / 2) ;
	 
    if (j < m/2) {
      name = "unused" ; 
      return ; 
    }
    ostringstream ostr ;
    ostr << "P_" << j << '^' << m << flush ;
    name = ostr.str() ;
  }	

  void basename_t_leg_mi(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    int m = 2 * ((k-1) / 2) + 1 ;
	 
    if (j < m/2) {
      name = "unused" ; 
      return ; 
    }
	
    ostringstream ostr ;
    ostr << "P_" << j << '^' << m << flush ;
    name = ostr.str() ;
  }	

  void basename_t_leg_p(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    int m = k / 2 ;
	 
    if (j < m/2) {
      name = "unused" ; 
      return ; 
    }

    int xt = (m%2 == 0) ? 2*j : 2*j + 1 ; 
    ostringstream ostr ;
    ostr << "P_" << xt << '^' << m << flush ;
    name = ostr.str() ;
  }	


  void basename_t_leg_pp(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    int m = 2 * (k / 2) ; 
	 
    if (j < m/2) {
      name = "unused" ; 
      return ; 
    }
	
    ostringstream ostr ;
    ostr << "P_" << 2*j << '^' << m << flush ;
    name = ostr.str() ;
  }	


  void basename_t_leg_i(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    int m = k / 2 ;
	 
    if (j < m/2 + m%2) {
      name = "unused" ; 
      return ; 
    }
	
    int xt = (m%2 == 0) ? 2*j + 1 : 2*j ; 
    ostringstream ostr ;
    ostr << "P_" << xt << '^' << m << flush ;
    name = ostr.str() ;
  }	
	

  void basename_t_leg_ip(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    int m = 2 * (k / 2) ; 
	 
    if (j < m/2) {
      name = "unused" ; 
      return ; 
    }
	
    ostringstream ostr ;
    ostr << "P_" << (2*j + 1) << '^' << m << flush ;
    name = ostr.str() ;
  }	


  void basename_t_leg_pi(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    int m = 2 * ((k-1) / 2) + 1 ; 
	 
    if (j < m/2) {
      name = "unused" ; 
      return ; 
    }
	
    ostringstream ostr ;
    ostr << "P_" << (2*j + 1) << '^' << m << flush ;
    name = ostr.str() ;
  }	


  void basename_t_leg_ii(int k, int j, string& name) {

    assert( k>=0 ) ; 
    assert( j>=0 ) ; 

    int m = 2 * ((k-1) / 2) + 1 ; 
	 
    if (j < m/2 + 1) {
      name = "unused" ; 
      return ; 
    }
	
    ostringstream ostr ;
    ostr << "P_" << 2*j << '^' << m << flush ;
    name = ostr.str() ;

  }	

	
}
