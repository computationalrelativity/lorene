/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char map_af_integ_surf_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2009/10/08 16:20:47  j_novak
 * Addition of new bases T_COS and T_SIN.
 *
 * Revision 1.4  2007/10/05 15:56:19  j_novak
 * Addition of new bases for the non-symmetric case in theta.
 *
 * Revision 1.3  2004/03/10 12:43:06  jl_jaramillo
 * Treatment of case ETATUN in surface integrals for Scalar's.
 *
 * Revision 1.2  2004/01/29 08:50:03  p_grandclement
 * Modification of Map::operator==(const Map&) and addition of the surface
 * integrales using Scalar.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.7  2001/02/19  11:40:27  phil
 * correction indices
 *
 * Revision 1.6  2001/02/12  14:14:29  phil
 * *** empty log message ***
 *
 * Revision 1.5  2001/02/12  14:00:52  phil
 * on prends tous les coefficients now
 *
 * Revision 1.4  2001/02/12  12:35:34  phil
 * gestion des bases angulaires plus proprement
 *
 * Revision 1.3  2001/01/02  10:52:27  phil
 * ajout calcul a l'infini
 *
 * Revision 1.2  2000/09/19  13:54:14  phil
 * *** empty log message ***
 *
 * Revision 1.1  2000/09/19  13:09:32  phil
 * Initial revision
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
#include "scalar.h"

                      //=============
                     // Cmp version
                    //===============  

double Map_af::integrale_surface (const Cmp& ci, double rayon) const {
    
    assert (ci.get_etat() != ETATNONDEF) ;
    assert (rayon > 0) ;
    if (ci.get_etat() == ETATZERO)
	return 0 ;
    
    assert (ci.get_etat() == ETATQCQ) ;
    
    int l ;
    double xi ;
    val_lx (rayon, 0, 0, l, xi) ;
    
    if (l == get_mg()->get_nzone()-1) {
	ci.check_dzpuis(0) ;
    }
    
    ci.va.coef() ;
    int nr = get_mg()->get_nr(l) ;
    int nt = get_mg()->get_nt(l) ;
    
    int base_r = ci.va.base.get_base_r(l) ;
    int base_t = ci.va.base.get_base_t(l) ;
    int base_p = ci.va.base.get_base_p(l) ;
    
    double result = 0 ;
    double* coef = new double [nr] ;
    double* auxi = new double[1] ;
     
    bool odd_theta = false ;
    double c_cos = 2 ;
    switch (base_t) {
	case T_COS_P : case T_COSSIN_CP :
	    break ;
	case T_COS_I : case T_COSSIN_CI :
	    odd_theta = true ;
	    break ;
	case T_COS : case T_COSSIN_C :
	    c_cos = 1. ;
	    break ;
	default :
	    cout << "base_t cas non prevu dans Map_af::integrale_surface" << endl ;
	    abort() ;
	    break ;
    }
    
    if (!odd_theta) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++)
		coef[i] = (*ci.va.c_cf)(l, 0, j, i) ;
	
	    switch (base_r) {
	    
		case R_CHEB :
		    som_r_cheb (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBP :
		    som_r_chebp (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBI :
		    som_r_chebi (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBU :
		    som_r_chebu (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBPI_P :
		    som_r_chebpi_p (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBPI_I :
		    som_r_chebpi_i (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBPIM_P :
		    som_r_chebpim_p (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBPIM_I :
		    som_r_chebpim_i (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		default :
		    som_r_pas_prevu (coef, nr, 1, 1, xi, auxi) ;
		    break ;
	    }
	    result += 2 * (*auxi)/(1-c_cos*c_cos*j*j) ;
	    if (c_cos == 1.) j++ ;
	}
    }
    delete [] auxi ;
    delete [] coef ;
	
     switch (base_p) {
	case P_COSSIN :
	    result *= 2*rayon*rayon*M_PI ;
	    break ;
	case P_COSSIN_P :
	    result *= 2*rayon*rayon*M_PI ;
	    break ;
	default :
	    cout << "base_p cas non prevu dans Map_af::integrale_surface" << endl ;
	    abort() ;
	    break ;
    }

    return result ;
}


// Integrale a l'infini
double Map_af::integrale_surface_infini (const Cmp& ci) const {
    
    assert (ci.get_etat() != ETATNONDEF) ;
    assert (ci.check_dzpuis(2));
    
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
      
    bool odd_theta = false ;
    double c_cos = 2. ;
    switch (base_t) {
	case T_COS_P : case T_COSSIN_CP :
	    break ;
	case T_COS_I : case T_COSSIN_CI :
	    odd_theta = true ;
	    break ;
	case T_COS : case T_COSSIN_C :
	    c_cos = 1. ;
	    break ;
	default :
	    cout << "base_t cas non prevu dans Map_af::integrale_surface" << endl ;
	    abort() ;
	    break ;
    }
    
    if (!odd_theta) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++)
		coef[i] = (*ci.va.c_cf)(nz-1, 0, j, i) ;
	    
	    switch (base_r) {
		case R_CHEBU :
		    som_r_chebu (coef, nr, 1, 1, 1, auxi) ;
		    break ;
		default :
		    som_r_pas_prevu (coef, nr, 1, 1, 1, auxi) ;
		    break ;
	    }
	    result += 2 * (*auxi)/(1-c_cos*c_cos*j*j) ;
	    if (c_cos == 1.) j++ ;
	}
    }
    delete [] auxi ;
    delete [] coef ;
    
    switch (base_p) {
	case P_COSSIN :
	    result *= 2*M_PI ;
	    break ;
	case P_COSSIN_P :
	    result *= 2*M_PI ;
	    break ;
	default :
	    cout << "base_p cas non prevu dans Map_af::integrale_surface_infini" << endl ;
	    abort() ;
	    break ;
    }
    
    return result ;
}

                      //=============
                     // Scalar version
                    //===============  

double Map_af::integrale_surface (const Scalar& ci, double rayon) const {
    
    assert (ci.get_etat() != ETATNONDEF) ;
    assert (rayon > 0) ;
    if (ci.get_etat() == ETATZERO)
	return 0 ;
    
    assert ( (ci.get_etat() == ETATQCQ) ||  (ci.get_etat() == ETATUN) ) ;
    
    int l ;
    double xi ;
    val_lx (rayon, 0, 0, l, xi) ;
    
    if (l == get_mg()->get_nzone()-1) {
	ci.check_dzpuis(0) ;
    }
    
    ci.get_spectral_va().coef() ;
    int nr = get_mg()->get_nr(l) ;
    int nt = get_mg()->get_nt(l) ;
    
    int base_r = ci.get_spectral_va().base.get_base_r(l) ;
    int base_t = ci.get_spectral_va().base.get_base_t(l) ;
    int base_p = ci.get_spectral_va().base.get_base_p(l) ;
    
    double result = 0 ;
    double* coef = new double [nr] ;
    double* auxi = new double[1] ;

    bool odd_theta = false ;
    double c_cos = 2. ;
     
    switch (base_t) {
	case T_COS_P : case T_COSSIN_CP :
	    break ;
	case T_COS_I: case T_COSSIN_CI :
	    odd_theta = true ; 
	    break ;
	case T_COS : case T_COSSIN_C :
	    c_cos = 1. ;
	    break ;
	default :
	    cout << "base_t cas non prevu dans Map_af::integrale_surface" << endl ;
	    abort() ;
	    break ;
    }
    
    if (!odd_theta) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++)
		coef[i] = (*ci.get_spectral_va().c_cf)(l, 0, j, i) ;
	    
	    switch (base_r) {
		
		case R_CHEB :
		    som_r_cheb (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBP :
		    som_r_chebp (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBI :
		    som_r_chebi (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBU :
		    som_r_chebu (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBPI_P :
		    som_r_chebpi_p (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBPI_I :
		    som_r_chebpi_i (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBPIM_P :
		    som_r_chebpim_p (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		case R_CHEBPIM_I :
		    som_r_chebpim_i (coef, nr, 1, 1, xi, auxi) ;
		    break ;
		default :
		    som_r_pas_prevu (coef, nr, 1, 1, xi, auxi) ;
		    break ;
	    }
	    result += 2 * (*auxi)/(1-c_cos*c_cos*j*j) ;
	    if (c_cos == 1.) j++ ;
	}
    }
	
    delete [] auxi ;
    delete [] coef ;
	
     switch (base_p) {
	case P_COSSIN :
	    result *= 2*rayon*rayon*M_PI ;
	    break ;
	case P_COSSIN_P :
	    result *= 2*rayon*rayon*M_PI ;
	    break ;
	default :
	    cout << "base_p cas non prevu dans Map_af::integrale_surface" << endl ;
	    abort() ;
	    break ;
    }

    return result ;
}


// Integrale a l'infini
double Map_af::integrale_surface_infini (const Scalar& ci) const {
    
    assert (ci.get_etat() != ETATNONDEF) ;
    assert (ci.check_dzpuis(2));
    
    if (ci.get_etat() == ETATZERO)
	return 0 ;
    
    assert ( (ci.get_etat() == ETATQCQ) ||  (ci.get_etat() == ETATUN) ) ; 
   
    int nz = ci.get_mp().get_mg()->get_nzone() ;
    
    ci.get_spectral_va().coef() ;
    int nr = get_mg()->get_nr(nz-1) ;
    int nt = get_mg()->get_nt(nz-1) ;
    
    int base_r = ci.get_spectral_va().base.get_base_r(nz-1) ;
    int base_t = ci.get_spectral_va().base.get_base_t(nz-1) ;
    int base_p = ci.get_spectral_va().base.get_base_p(nz-1) ;
    
    double result = 0 ;
    double* coef = new double [nr] ;
    double* auxi = new double[1] ;
      
    bool odd_theta = false ;
    double c_cos = 2. ;
     
    switch (base_t) {
	case T_COS_P : case T_COSSIN_CP :
	    break ;
	case T_COS_I: case T_COSSIN_CI :
	    odd_theta = true ; 
	    break ;
	case T_COS : case T_COSSIN_C :
	    c_cos = 1. ;
	    break ;
	default :
	    cout << "base_t cas non prevu dans Map_af::integrale_surface" << endl ;
	    abort() ;
	    break ;
    }
    
    if (!odd_theta) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++)
		coef[i] = (*ci.get_spectral_va().c_cf)(nz-1, 0, j, i) ;
	    
	    switch (base_r) {
		case R_CHEBU :
		    som_r_chebu (coef, nr, 1, 1, 1, auxi) ;
		    break ;
		default :
		    som_r_pas_prevu (coef, nr, 1, 1, 1, auxi) ;
		    break ;
	    }
	    result += 2 * (*auxi)/(1-c_cos*c_cos*j*j) ;
	    if (c_cos == 1.) j++ ;
	}
    }
	
    delete [] auxi ;
    delete [] coef ;
    
     switch (base_p) {
	case P_COSSIN :
	    result *= 2*M_PI ;
	    break ;
	case P_COSSIN_P :
	    result *= 2*M_PI ;
	    break ;
	default :
	    cout << "base_p cas non prevu dans Map_af::integrale_surface_infini" << endl ;
	    abort() ;
	    break ;
    }
    
    return result ;
}
