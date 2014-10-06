/*
 *   Copyright (c) 2001 Jerome Novak
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


char base_val_manip_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2014/10/06 15:12:56  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2009/10/08 16:20:13  j_novak
 * Addition of new bases T_COS and T_SIN.
 *
 * Revision 1.8  2008/12/03 15:21:21  j_novak
 * New method mult_cost.
 *
 * Revision 1.7  2008/02/18 13:53:38  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.6  2004/11/23 15:08:00  m_forot
 * Added the bases for the cases without any equatorial symmetry
 * (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.5  2004/01/27 14:31:26  j_novak
 * New method Base_val::mult_sint()
 *
 * Revision 1.4  2004/01/27 14:13:59  j_novak
 * Added method Base_val::mult_x()
 *
 * Revision 1.3  2003/09/16 08:54:09  j_novak
 * Addition of the T_LEG_II base (odd in theta, only for odd m) and the
 * transformation functions to and from the T_SIN_P base.
 *
 * Revision 1.2  2002/10/16 14:36:30  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2001/10/12  15:06:50  novak
 * *** empty log message ***
 *
 * Revision 1.2  2001/10/12 15:05:14  novak
 * *** empty log message ***
 *
 * Revision 1.1  2001/10/12 14:57:24  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cassert>

// Headers Lorene
#include "headcpp.h"
#include "type_parite.h"
#include "base_val.h"

void Base_val::dsdx() {
    
    switch(get_base_r(0)) {
	case R_CHEBP: 
	    set_base_r(0, R_CHEBI) ;
	    break ;
	case R_CHEBI: 
	    set_base_r(0, R_CHEBP) ;
	    break ;
	case R_CHEBPIM_P: 
	    set_base_r(0, R_CHEBPIM_I) ;
	    break ;
	case R_CHEBPIM_I: 
	    set_base_r(0, R_CHEBPIM_P) ;
	    break ;
	case R_CHEBPI_P: 
	    set_base_r(0, R_CHEBPI_I) ;
	    break ;
	case R_CHEBPI_I: 
	    set_base_r(0, R_CHEBPI_P) ;
	    break ;    
	default: 
	    break ;
    }  
    return ;
}

void Base_val::sx() {
    
    switch(get_base_r(0)) {
	case R_CHEBP: 
	    set_base_r(0, R_CHEBI) ;
	    break ;
	case R_CHEBI: 
	    set_base_r(0, R_CHEBP) ;
	    break ;
	case R_CHEBPIM_P: 
	    set_base_r(0, R_CHEBPIM_I) ;
	    break ;
	case R_CHEBPIM_I: 
	    set_base_r(0, R_CHEBPIM_P) ;
	    break ;
	case R_CHEBPI_P: 
	    set_base_r(0, R_CHEBPI_I) ;
	    break ;
	case R_CHEBPI_I: 
	    set_base_r(0, R_CHEBPI_P) ;
	    break ;
	default: 
	    break ;
    }  
    return ;
}

void Base_val::mult_x() {
    
    switch(get_base_r(0)) {
	case R_CHEBP: 
	    set_base_r(0, R_CHEBI) ;
	    break ;
	case R_CHEBI: 
	    set_base_r(0, R_CHEBP) ;
	    break ;
	case R_CHEBPIM_P: 
	    set_base_r(0, R_CHEBPIM_I) ;
	    break ;
	case R_CHEBPIM_I: 
	    set_base_r(0, R_CHEBPIM_P) ;
	    break ;
	case R_CHEBPI_P: 
	    set_base_r(0, R_CHEBPI_I) ;
	    break ;
	case R_CHEBPI_I: 
	    set_base_r(0, R_CHEBPI_P) ;
	    break ;
	default: 
	    break ;
    }  
    return ;
}

void Base_val::dsdt() {
    
    switch(get_base_t(0)) {
	case T_COS_P:
	    set_base_t(T_SIN_P) ;
	    break ;
	case T_COS_I:
	    set_base_t(T_SIN_I) ;
	    break ;
	case T_SIN_P:
	    set_base_t(T_COS_P) ;
	    break ;
	case T_SIN_I:
	    set_base_t(T_COS_I) ;
	    break ;
	case T_COSSIN_CP:
	    set_base_t(T_COSSIN_SP) ;
	    break ;
	case T_COSSIN_SP:
	    set_base_t(T_COSSIN_CP) ;
	    break ;
	case T_COSSIN_CI:
	    set_base_t(T_COSSIN_SI) ;
	    break ;
	case T_COSSIN_SI:
	    set_base_t(T_COSSIN_CI) ;
	    break ;
	case T_COSSIN_C:
	    set_base_t(T_COSSIN_S) ;
	    break ;
	case T_COSSIN_S:
	    set_base_t(T_COSSIN_C) ;
	    break ;
	case T_COS:
	    set_base_t(T_SIN) ;
	    break ;
	case T_SIN:
	    set_base_t(T_COS) ;
	    break ;
	default: 
	    cout << "Wrong base in Base_val::dsdt()!" << endl ;
	    abort() ;
	    exit(-1) ;
	    break ;
    }  
    return ;
}

void Base_val::ssint() {
    
    switch(get_base_t(0)) {
	case T_COS_P:
	    set_base_t(T_SIN_I) ;
	    break ;
	case T_COS_I:
	    set_base_t(T_SIN_P) ;
	    break ;
	case T_SIN_P:
	    set_base_t(T_COS_I) ;
	    break ;
	case T_SIN_I:
	    set_base_t(T_COS_P) ;
	    break ;
	case T_COSSIN_CP:
	    set_base_t(T_COSSIN_SI) ;
	    break ;
	case T_COSSIN_SP:
	    set_base_t(T_COSSIN_CI) ;
	    break ;
	case T_COSSIN_CI:
	    set_base_t(T_COSSIN_SP) ;
	    break ;
	case T_COSSIN_SI:
	    set_base_t(T_COSSIN_CP) ;
	    break ;
	case T_COSSIN_C:
	    set_base_t(T_COSSIN_S) ;
	    break ;
	case T_COSSIN_S:
	    set_base_t(T_COSSIN_C) ;
	    break ;  
	case T_COS:
	    set_base_t(T_SIN) ;
	    break ;
	case T_SIN:
	    set_base_t(T_COS) ;
	    break ;  
	default: 
	    cout << "Wrong base in Base_val::ssint()!" << endl ;
	    abort() ;
	    exit(-1) ;
	    break ;
    }  
    return ;
}

void Base_val::mult_sint() {
    
    switch(get_base_t(0)) {
	case T_COS_P:
	    set_base_t(T_SIN_I) ;
	    break ;
	case T_COS_I:
	    set_base_t(T_SIN_P) ;
	    break ;
	case T_SIN_P:
	    set_base_t(T_COS_I) ;
	    break ;
	case T_SIN_I:
	    set_base_t(T_COS_P) ;
	    break ;
	case T_COSSIN_CP:
	    set_base_t(T_COSSIN_SI) ;
	    break ;
	case T_COSSIN_SP:
	    set_base_t(T_COSSIN_CI) ;
	    break ;
	case T_COSSIN_CI:
	    set_base_t(T_COSSIN_SP) ;
	    break ;
	case T_COSSIN_SI:
	    set_base_t(T_COSSIN_CP) ;
	    break ;
	case T_COSSIN_C:
	    set_base_t(T_COSSIN_S) ;
	    break ;
	case T_COSSIN_S:
	    set_base_t(T_COSSIN_C) ;
	    break ;   
	case T_COS:
	    set_base_t(T_SIN) ;
	    break ;
	case T_SIN:
	    set_base_t(T_COS) ;
	    break ;  
	default: 
	    cout << "Wrong base in Base_val::mult_sint()!" << endl ;
	    abort() ;
	    exit(-1) ;
	    break ;
    }  
    return ;
}

void Base_val::mult_cost() {
    
    switch(get_base_t(0)) {
	case T_COS_P:
	    set_base_t(T_COS_I) ;
	    break ;
	case T_COS_I:
	    set_base_t(T_COS_P) ;
	    break ;
	case T_SIN_P:
	    set_base_t(T_SIN_I) ;
	    break ;
	case T_SIN_I:
	    set_base_t(T_SIN_P) ;
	    break ;
	case T_COSSIN_CP:
	    set_base_t(T_COSSIN_CI) ;
	    break ;
	case T_COSSIN_SP:
	    set_base_t(T_COSSIN_SI) ;
	    break ;
	case T_COSSIN_CI:
	    set_base_t(T_COSSIN_CP) ;
	    break ;
	case T_COSSIN_SI:
	    set_base_t(T_COSSIN_SP) ;
	    break ;
	case T_COSSIN_C:
	    set_base_t(T_COSSIN_C) ;
	    break ;
	case T_COSSIN_S:
	    set_base_t(T_COSSIN_S) ;
	    break ;   
	case T_COS:
	    set_base_t(T_COS) ;
	    break ;
	case T_SIN:
	    set_base_t(T_SIN) ;
	    break ;  
	default: 
	    cout << "Wrong base in Base_val::mult_cost()!" << endl ;
	    abort() ;
	    exit(-1) ;
	    break ;
    }  
    return ;
}

void Base_val::ylm() {
    
    switch(get_base_t(0)) {
	case T_COS_P:
	    set_base_t(T_LEG_PP) ;
	    break ;
	case T_COS_I:
	    set_base_t(T_LEG_IP) ;
	    break ;
	case T_SIN_I:
	    set_base_t(T_LEG_PI) ;
	    break ;
	case T_SIN_P:
	    set_base_t(T_LEG_II) ;
	    break ;
	case T_COSSIN_CP:
	    set_base_t(T_LEG_P) ;
	    break ;
	case T_COSSIN_CI:
	    set_base_t(T_LEG_I) ;
	    break ;
	case T_COSSIN_C:
	    set_base_t(T_LEG) ;
	    break ;
	case T_COSSIN_S:
	    set_base_t(T_LEG) ;
	    break ;
	case T_COS:
	    set_base_t(T_LEG_MP) ;
	    break ;
	default: 
	    cout << "Wrong base in Base_val::ylm()!" << endl ;
	    abort() ;
	    exit(-1) ;
	    break ;
    }  
    return ;
}
