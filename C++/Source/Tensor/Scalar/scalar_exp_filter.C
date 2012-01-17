/*
 *  Function applying an exponential filter to a Scalar: 
 *  sigma(n/N) = exp(alpha*(n/N)^(2p)). See scalar.h for documentation.
 */

/*
 *   Copyright (c) 2007  Jerome Novak
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

char scalar_exp_filter_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2012/01/17 10:29:27  j_penner
 * added two routines to handle generalized exponential filtering
 *
 * Revision 1.2  2007/10/31 10:50:16  j_novak
 * Testing whether the coefficients are zero in a given domain.
 *
 * Revision 1.1  2007/10/31 10:33:13  j_novak
 * Added exponential filters to smooth Gibbs-type phenomena.
 *
 *
 * $Header$
 *
 */

// C headers
#include <assert.h>
#include <math.h>

// Lorene headers
#include "tensor.h"
#include "proto.h"

void Scalar::exponential_filter_r(int lzmin, int lzmax, int p, 
			  double alpha) {
    assert(lzmin >= 0) ;
    const Mg3d& mgrid = *mp->get_mg() ;
#ifndef NDEBUG
    int nz = mgrid.get_nzone() ;
#endif
    assert(lzmax < nz) ;
    assert(etat != ETATNONDEF) ;
    if (etat == ETATZERO) return ;
    va.coef() ;
    assert(va.c_cf != 0x0) ;
    assert(alpha < 0.) ;
    double alp = log(pow(10., alpha)) ;
    
    for (int lz=lzmin; lz<=lzmax; lz++) {
	if ((*va.c_cf)(lz).get_etat() == ETATQCQ) 
	    for (int k=0; k<mgrid.get_np(lz); k++) 
		for (int j=0; j<mgrid.get_nt(lz); j++) {
		    int nr = mgrid.get_nr(lz) ;
		    for (int i=0; i<nr; i++) {
			double eta = double(i)/double(nr) ;
			va.c_cf->set(lz, k, j, i) *= exp(alp*pow(eta, 2*p)) ;
		    }
		}
    }
     if (va.c != 0x0) {
 	delete va.c ;
 	va.c = 0x0 ;
     }
     va.del_deriv() ;
     del_deriv() ;
    
    return ;
}

void Scalar::sarra_filter_r(int lzmin, int lzmax, double p, 
			  double alpha) {
    assert(lzmin >= 0) ;
    const Mg3d& mgrid = *mp->get_mg() ;
#ifndef NDEBUG
    int nz = mgrid.get_nzone() ;
#endif
    assert(lzmax < nz) ;
    assert(etat != ETATNONDEF) ;
    if (etat == ETATZERO) return ;
    va.coef() ;
    assert(va.c_cf != 0x0) ;
    assert(alpha < 0.) ;
    
    for (int lz=lzmin; lz<=lzmax; lz++) {
	if ((*va.c_cf)(lz).get_etat() == ETATQCQ) 
	    for (int k=0; k<mgrid.get_np(lz); k++) 
		for (int j=0; j<mgrid.get_nt(lz); j++) {
		    int nr = mgrid.get_nr(lz) ;
		    for (int i=0; i<nr; i++) {
			double eta = double(i)/double(nr) ;
			va.c_cf->set(lz, k, j, i) *= exp(alpha*pow(eta, p)) ;
		    }
		}
    }
     if (va.c != 0x0) {
 	delete va.c ;
 	va.c = 0x0 ;
     }
     va.del_deriv() ;
     del_deriv() ;
    
    return ;
}
void exp_filter_r_all_domains( Scalar& ss, int p, double alpha ) {
    int nz = ss.get_mp().get_mg()->get_nzone() ;
    ss.exponential_filter_r(0, nz-1, p, alpha) ;
    return ;
}

void Scalar::sarra_filter_r_all_domains( double p, double alpha ) {
    int nz = get_mp().get_mg()->get_nzone() ;
    sarra_filter_r(0, nz-1, p, alpha) ;
    return ;
}

void Scalar::exponential_filter_ylm(int lzmin, int lzmax, int p, 
			  double alpha) {
    assert(lzmin >= 0) ;
    const Mg3d& mgrid = *mp->get_mg() ;
#ifndef NDEBUG
    int nz = mgrid.get_nzone() ;
#endif
    assert(lzmax < nz) ;
    assert(etat != ETATNONDEF) ;
    if (etat == ETATZERO) return ;
    double alp = log(pow(10., alpha)) ;
    va.ylm() ;
    assert(va.c_cf != 0x0) ;
    const Base_val& base = va.base ;
    int l_q, m_q, base_r ;
    
    for (int lz=lzmin; lz<=lzmax; lz++) 
	if ((*va.c_cf)(lz).get_etat() == ETATQCQ) {
	    int np = mgrid.get_np(lz) ;
	    int nt = mgrid.get_nt(lz) ;
	    int nr = mgrid.get_nr(lz) ;
	    int lmax = base.give_lmax(mgrid, lz) ; 
	    for (int k=0; k<np; k++) 
		for (int j=0; j<nt; j++) {
		    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		    if (nullite_plm(j, nt, k, np, base) == 1 ) {
			double eta = double(l_q) / double(lmax) ;
			double sigma = exp(alp*pow(eta, 2*p)) ;
			for (int i=0; i<nr; i++) 
			    va.c_cf->set(lz, k, j, i) *= sigma ;
		    }
		}
	}
    
    va.ylm_i() ;
    if (va.c != 0x0) {
 	delete va.c ;
 	va.c = 0x0 ;
    }
    va.del_deriv() ;
    del_deriv() ;
    return ;
}

void exp_filter_ylm_all_domains(Scalar& ss, int p, double alpha ) {
    int nz = ss.get_mp().get_mg()->get_nzone() ;
    ss.exponential_filter_ylm(0, nz-1, p, alpha) ;
    return ;
}

