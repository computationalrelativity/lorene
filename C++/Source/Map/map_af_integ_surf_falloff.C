/*
 *  Methods of the class Map_af for obtaining a surface integral
 *   with a falloff condition at the outer boundary
 *
 *    (see file map.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Joshua A. Faber
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

char map_af_integ_surf_falloff_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/11/30 20:53:08  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

#include <stdlib.h>
#include <math.h>

#include "map.h"
#include "cmp.h"
#include "proto.h"

double Map_af::integrale_surface_falloff (const Cmp& ci) const {
    
  // Returns Surface integral/R^2 -> multiply by R^2 to get the right value!!!!!!

    assert (ci.get_etat() != ETATNONDEF) ;
    if (ci.get_etat() == ETATZERO)
	return 0 ;
    
    assert (ci.get_etat() == ETATQCQ) ;
    
    int nz = ci.get_mp()->get_mg()->get_nzone() ;
    
    ci.va.coef() ;
    int nr = get_mg()->get_nr(nz-1) ;
    int nt = get_mg()->get_nt(nz-1) ;
    
    int base_r = ci.va.base.get_base_r(nz-1) ;
    int base_t = ci.va.base.get_base_t(nz-1) ;
    int base_p = ci.va.base.get_base_p(nz-1) ;
    
    double result = 0 ;
    double* coef = new double [nr] ;
    double* auxi = new double[1] ;
     
    for (int j=0 ; j<nt ; j++) {
	for (int i=0 ; i<nr ; i++)
	    coef[i] = (*ci.va.c_cf)(nz-1, 0, j, i) ;
  
	switch (base_r) {
	
	    case R_CHEB :
		som_r_cheb (coef, nr, 1, 1, 1, auxi) ;
		break ;
	    case R_CHEBP :
		som_r_chebp (coef, nr, 1, 1, 1, auxi) ;
		break ;
	    case R_CHEBI :
		som_r_chebi (coef, nr, 1, 1, 1, auxi) ;
		break ;
	    case R_CHEBU :
		som_r_chebu (coef, nr, 1, 1, 1, auxi) ;
		break ;
	    case R_CHEBPIM_P :
		som_r_chebpim_p (coef, nr, 1, 1, 1, auxi) ;
		break ;
	    case R_CHEBPIM_I :
		som_r_chebpim_i (coef, nr, 1, 1, 1, auxi) ;
		break ;
	    default :
		som_r_pas_prevu (coef, nr, 1, 1, 1, auxi) ;
		break ;
	}
	result += 2 * (*auxi)/(1-4*j*j) ;
	}
	
    delete [] auxi ;
    delete [] coef ;
	
    switch (base_t) {
	case T_COS_P :
	    break ;
	case T_COSSIN_CP :
	    break ;
	case T_COSSIN_CI :
	    result = 0 ;
	    break ;
	default :
	    cout << "base_t cas non prevu dans Map_af::integrale_surface" << endl ;
	    abort() ;
	    break ;
    }
    
     switch (base_p) {
	case P_COSSIN :
	    result *= 2*M_PI ;
	    break ;
	case P_COSSIN_P :
	    result *= 2*M_PI ;
	    break ;
	default :
	    cout << "base_p cas non prevu dans Map_af::integrale_surface" << endl ;
	    abort() ;
	    break ;
    }

    return (result) ;
}

