/*
 *  Method Base_val::name_r
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
 * Revision 1.8  2023/06/28 10:04:32  j_novak
 * Use of C++ strings and flows instead of C types.
 *
 * Revision 1.7  2016/12/05 16:17:44  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:52:39  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:12:57  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2013/01/11 08:20:11  j_novak
 * New radial spectral bases with Legendre polynomials (R_LEG, R_LEGP, R_LEGI).
 *
 * Revision 1.3  2007/12/11 15:28:09  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.2  2004/11/23 15:08:01  m_forot
 * Added the bases for the cases without any equatorial symmetry
 * (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
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
#include "base_val.h"

// Local prototypes
namespace Lorene {
  void basename_r_unknown(int, int, int, string&) ; 
  void basename_r_cheb(int, int, int, string&) ; 
  void basename_r_chebp(int, int, int, string&) ; 
  void basename_r_chebi(int, int, int, string&) ; 
  void basename_r_chebpim_p(int, int, int, string&) ; 
  void basename_r_chebpim_i(int, int, int, string&) ; 
  void basename_r_chebpi_p(int, int, int, string&) ; 
  void basename_r_chebpi_i(int, int, int, string&) ; 
  void basename_r_leg(int, int, int, string&) ; 
  void basename_r_legp(int, int, int, string&) ; 
  void basename_r_legi(int, int, int, string&) ; 
  void basename_r_jaco02(int, int, int, string&) ;

  //----------------------------//
  //      Base_val method       //
  //----------------------------//

  void Base_val::name_r(int l, int k, int j, int i, string& name) const {

    // Array of actual base name functions
    static void(*vbasename_r[MAX_BASE])(int, int, int, string&) ;  

    static bool first_call = true ;

    // Initializations at first call
    // -----------------------------
    if ( first_call ) {

      first_call = false ;

      for (int ib=0 ; ib<MAX_BASE ; ib++) {
	vbasename_r[ib] = basename_r_unknown ;
      }

      vbasename_r[R_CHEB >> TRA_R] = basename_r_cheb ;
      vbasename_r[R_CHEBP >> TRA_R] = basename_r_chebp ;
      vbasename_r[R_CHEBI >> TRA_R] = basename_r_chebi ;
      vbasename_r[R_CHEBPIM_P >> TRA_R] = basename_r_chebpim_p ;
      vbasename_r[R_CHEBPIM_I >> TRA_R] = basename_r_chebpim_i ;
      vbasename_r[R_CHEBU >> TRA_R] = basename_r_cheb ;
      vbasename_r[R_CHEBPI_P >> TRA_R] = basename_r_chebpi_p ;
      vbasename_r[R_CHEBPI_I >> TRA_R] = basename_r_chebpi_i ;
      vbasename_r[R_LEG >> TRA_R] = basename_r_leg ;
      vbasename_r[R_LEGP >> TRA_R] = basename_r_legp ;
      vbasename_r[R_LEGI >> TRA_R] = basename_r_legi ;
      vbasename_r[R_JACO02 >> TRA_R] = basename_r_jaco02 ;
    }
	
    // Call to the function adapted to the basis in domain l
    //------------------------------------------------------
	
    assert( (l>=0) && (l<nzone) ) ; 
	
    int base_r = ( b[l] & MSQ_R ) >> TRA_R ;
	
    vbasename_r[base_r](k, j, i, name) ; 

  }
	
	
  //-------------------------------//
  //  individual basis functions   //
  //-------------------------------//
	
  void basename_r_unknown(int, int, int, string&) {
    cout << "Base_val::name_r : unknwon basis !" << endl ; 
    abort() ; 
  } 


  void basename_r_cheb(int, int, int i, string& name) {

    assert( i>=0 ) ;

    ostringstream ostr ;
    ostr << 'T' << i << flush ;
    name = ostr.str() ;
  }	


  void basename_r_chebp(int, int, int i, string& name) {

    assert( i>=0 ) ; 

    ostringstream ostr ;
    ostr << 'T' << 2*i << flush ;
    name = ostr.str() ;
  }	


  void basename_r_chebi(int, int, int i, string& name) {

    assert( i>=0 ) ; 

    ostringstream ostr ;
    ostr << 'T' << (2*i+1) << flush ;
    name = ostr.str() ;
  }	


  void basename_r_chebpim_p(int k, int, int i, string& name) {

    assert( k>=0 ) ; 
    assert( i>=0 ) ; 

    int m = k / 2 ; 
    int xr = (m%2 == 0) ? 2*i : 2*i + 1 ; 

    ostringstream ostr ;
    ostr << 'T' << xr << flush ;
    name = ostr.str() ;
  }	


  void basename_r_chebpim_i(int k, int, int i, string& name) {

    assert( k>=0 ) ; 
    assert( i>=0 ) ; 

    int m = k / 2 ; 
    int xr = (m%2 == 0) ? 2*i + 1 : 2*i ; 

    ostringstream ostr ;
    ostr << 'T' << xr << flush ;
    name = ostr.str() ;
  }	

  void basename_r_chebpi_p(int , int j, int i, string& name) {

    assert( j>=0 ) ; 
    assert( i>=0 ) ; 

    int xr = (j%2 == 0) ? 2*i : 2*i + 1 ; 

    ostringstream ostr ;
    ostr << 'T' << xr << flush ;
    name = ostr.str() ;
  }	


  void basename_r_chebpi_i(int , int j, int i, string& name) {

    assert( j>=0 ) ; 
    assert( i>=0 ) ; 

    int xr = (j%2 == 0) ? 2*i + 1 : 2*i ; 

    ostringstream ostr ;
    ostr << 'T' << xr << flush ;
    name = ostr.str() ;
  }	

  void basename_r_leg(int, int, int i, string& name) {

    assert( i>=0 ) ; 

    ostringstream ostr ;
    ostr << 'P' << i << flush ;
    name = ostr.str() ;
  }	


  void basename_r_legp(int, int, int i, string& name) {

    assert( i>=0 ) ; 

    ostringstream ostr ;
    ostr << 'P' << 2*i << flush ;
    name = ostr.str() ;
  }	


  void basename_r_legi(int, int, int i, string& name) {

    assert( i>=0 ) ; 

    ostringstream ostr ;
    ostr << 'T' << (2*i + 1)  << flush ;
    name = ostr.str() ;
  }	

  void basename_r_jaco02(int, int, int i, string& name) {

    assert( i>=0 ) ; 

    ostringstream ostr ;
    ostr << 'J' << i << flush ;
    name = ostr.str() ;
  }	

}
