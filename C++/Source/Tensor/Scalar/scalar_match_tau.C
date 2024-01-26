/*
 * Function to perform a matching across domains + imposition of boundary conditions
 * at the outer domain and regularity condition at the center.
 */

/*
 *   Copyright (c) 2008  Jerome Novak
 *                 2011  Nicolas Vasset
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
 * Revision 1.8  2024/01/26 17:44:26  g_servignat
 * Updated the Pseudopolytrope_1D class to be consistent with the paper (i.e. with a GPP in the middle)
 *
 * Revision 1.7  2016/12/05 16:18:18  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:46  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:16:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2012/05/11 14:11:24  j_novak
 * Modifications to deal with the accretion of a scalar field
 *
 * Revision 1.3  2012/02/06 12:59:07  j_novak
 * Correction of some errors.
 *
 * Revision 1.2  2008/08/20 13:23:43  j_novak
 * The shift in quantum number l (e.g. for \tilde{B}) is now taken into account.
 *
 * Revision 1.1  2008/05/24 15:05:22  j_novak
 * New method Scalar::match_tau to match the output of an explicit time-marching scheme with the tau method.
 *
 *
 * $Header$
 *
 */

// C headers
#include <cassert>
#include <cmath>

// Lorene headers
#include "tensor.h"
#include "matrice.h"
#include "proto.h"
#include "param.h"

namespace Lorene {
void Scalar::match_tau(Param& par_bc, Param* par_mat) {

    const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
    assert(mp_aff != 0x0) ; //Only affine mapping for the moment

    const Mg3d& mgrid = *mp_aff->get_mg() ;
    assert(mgrid.get_type_r(0) == RARE)  ;
    assert (par_bc.get_n_tbl_mod() > 0) ;
    Tbl mu = par_bc.get_tbl_mod(0) ;
    if (etat == ETATZERO) {
	if (mu.get_etat() == ETATZERO) 
	    return ;
	else
	    annule_hard() ;
    }

    int nz = mgrid.get_nzone() ;
    int nzm1 = nz - 1 ;
    bool ced = (mgrid.get_type_r(nzm1) == UNSURR) ;
    assert(par_bc.get_n_int() >= 2);
    int domain_bc = par_bc.get_int(0) ;
    bool bc_ced = ((ced) && (domain_bc == nzm1)) ;

    int n_conditions = par_bc.get_int(1) ;
    assert ((n_conditions==1)||(n_conditions==2)) ;
    bool derivative = (n_conditions == 2) ;
    int dl = 0 ; int l_min = 0 ; int excised_i =0;
    if (par_bc.get_n_int() > 2) {
	dl = par_bc.get_int(2) ;
	l_min = par_bc.get_int(3) ;
	if(par_bc.get_n_int() >4){
	excised_i = par_bc.get_int(4);
	}
    }
    bool isexcised = (excised_i==1);
    
    int nt = mgrid.get_nt(0) ;
    int np = mgrid.get_np(0) ;
    assert (par_bc.get_n_double_mod() == 2) ;
    double c_val = par_bc.get_double_mod(0) ;
    double c_der = par_bc.get_double_mod(1) ;
    if (bc_ced) {
	c_val = 1 ;
	c_der = 0 ;
	mu.annule_hard() ;
    }
    if (mu.get_etat() == ETATZERO)
	mu.annule_hard() ;
    int system01_size =0;
    int system_size =0;
    if (isexcised ==false){
    system01_size = 1 ;
    system_size = 2 ;
    }
    else{
      system01_size = -1;
      system_size = -1;
    }
    for (int lz=1; lz<=domain_bc; lz++) {
	    system01_size += n_conditions ;
	    system_size += n_conditions ;
    }
    assert (etat != ETATNONDEF) ;

    bool need_matrices = true ;
    if (par_mat != 0x0)
	if (par_mat->get_n_matrice_mod() == 4) 
	    need_matrices = false ;
    
    Matrice system_l0(system01_size, system01_size) ;
    Matrice system_l1(system01_size, system01_size) ;
    Matrice system_even(system_size, system_size) ;
    Matrice system_odd(system_size, system_size) ;    
    
    if (need_matrices) {
	system_l0.annule_hard() ;
	system_l1.annule_hard() ;
	system_even.annule_hard() ;
	system_odd.annule_hard() ;
	int index = 0 ; int index01 = 0 ;
	int nr = mgrid.get_nr(0);
	double alpha = mp_aff->get_alpha()[0] ; 
	if (isexcised == false){
	  system_even.set(index, index) = 1. ;
	  system_even.set(index, index + 1) = -1. ; //C_{N-1} - C_N = \sum (-1)^i C_i
	  system_odd.set(index, index) = -(2.*nr - 5.)/alpha ; 
	  system_odd.set(index, index+1) = (2.*nr - 3.)/alpha ;
	  index++ ;
	  if (domain_bc == 0) {
	    system_l0.set(index01, index01) = c_val + c_der*4.*(nr-1)*(nr-1)/alpha ;
	    system_l1.set(index01, index01) = c_val + c_der*(2*nr-3)*(2*nr-3)/alpha ;
	    index01++ ;
	    system_even.set(index, index-1) = c_val + c_der*4.*(nr-2)*(nr-2)/alpha ;
	    system_even.set(index, index) = c_val + c_der*4.*(nr-1)*(nr-1)/alpha ;
	    system_odd.set(index, index-1) = c_val + c_der*(2*nr-5)*(2*nr-5)/alpha ;
	    system_odd.set(index, index) = c_val + c_der*(2*nr-3)*(2*nr-3)/alpha ;
	    index++ ;
	  }
	  else {
	    system_l0.set(index01, index01) = 1. ;
	    system_l1.set(index01, index01) = 1. ;
	    system_even.set(index, index-1) = 1. ;
	    system_even.set(index, index) = 1. ;
	    system_odd.set(index, index-1) = 1. ;
	    system_odd.set(index, index) = 1. ;
	    if (derivative) {
	      system_l0.set(index01+1, index01) = 4*(nr-1)*(nr-1)/alpha ;
	      system_l1.set(index01+1, index01) = (2*nr-3)*(2*nr-3)/alpha ;
	      system_even.set(index+1, index-1) = 4*(nr-2)*(nr-2)/alpha ;
	      system_even.set(index+1, index) = 4*(nr-1)*(nr-1)/alpha ;
	      system_odd.set(index+1, index-1) = (2*nr-5)*(2*nr-5)/alpha ;
	      system_odd.set(index+1, index) = (2*nr-3)*(2*nr-3)/alpha ;
	    }
	  }
	  if (domain_bc > 0) {
	    //	Do things for lz=1;
	    nr = mgrid.get_nr(1) ;
	    alpha = mp_aff->get_alpha()[1] ;
	    if ((1 == domain_bc)&&(bc_ced))
		alpha = -0.25/alpha ;
	    system_l0.set(index01, index01+1) = 1. ;
	    system_l0.set(index01, index01+2) = -1. ;
	    system_l1.set(index01, index01+1) = 1. ;
	    system_l1.set(index01, index01+2) = -1. ;
	    index01++ ;
	    system_even.set(index, index+1) = 1. ;
	    system_even.set(index, index+2) = -1. ;
	    system_odd.set(index, index+1) = 1. ;
	    system_odd.set(index, index+2) = -1. ;
	    index++ ;
	    if (derivative) {
		system_l0.set(index01, index01) = -(nr-2)*(nr-2)/alpha ;
		system_l0.set(index01, index01+1) = (nr-1)*(nr-1)/alpha ;
		system_l1.set(index01, index01) = -(nr-2)*(nr-2)/alpha ;
		system_l1.set(index01, index01+1) = (nr-1)*(nr-1)/alpha ;
		index01++ ;
		system_even.set(index, index) = -(nr-2)*(nr-2)/alpha ;
		system_even.set(index, index+1) = (nr-1)*(nr-1)/alpha ;
		system_odd.set(index, index) = -(nr-2)*(nr-2)/alpha ;
		system_odd.set(index, index+1) = (nr-1)*(nr-1)/alpha ;
		index++ ;
	    }
	
	    if (1 == domain_bc) {
		system_l0.set(index01, index01-1) = c_val + c_der*(nr-2)*(nr-2)/alpha ;
		system_l0.set(index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
		system_l1.set(index01, index01-1) = c_val + c_der*(nr-2)*(nr-2)/alpha ;
		system_l1.set(index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
		index01++ ; 
		system_even.set(index, index-1) = c_val + c_der*(nr-2)*(nr-2)/alpha ;
		system_even.set(index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
		system_odd.set(index, index-1) = c_val + c_der*(nr-2)*(nr-2)/alpha ;
		system_odd.set(index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
		index++ ;
	    }
	    else {
		system_l0.set(index01, index01-1) = 1. ;
		system_l0.set(index01, index01) = 1. ;
		system_l1.set(index01, index01-1) = 1. ;
		system_l1.set(index01, index01) = 1. ;
		system_even.set(index, index-1) = 1. ;
		system_even.set(index, index) = 1. ;
		system_odd.set(index, index-1) = 1. ;
		system_odd.set(index, index) = 1. ;
		if (derivative) {
		    system_l0.set(index01+1, index01-1) = (nr-2)*(nr-2)/alpha ;
		    system_l0.set(index01+1, index01) = (nr-1)*(nr-1)/alpha ;
		    system_l1.set(index01+1, index01-1) = (nr-2)*(nr-2)/alpha ;
		    system_l1.set(index01+1, index01) = (nr-1)*(nr-1)/alpha ;
		    system_even.set(index+1, index-1) = (nr-2)*(nr-2)/alpha ;
		    system_even.set(index+1, index) = (nr-1)*(nr-1)/alpha ;
		    system_odd.set(index+1, index-1) = (nr-2)*(nr-2)/alpha ;
		    system_odd.set(index+1, index) = (nr-1)*(nr-1)/alpha ;
		}		
	    }
	  }
	}
	else {
	  nr = mgrid.get_nr(1) ;
	  alpha = mp_aff->get_alpha()[1] ;
	  if ((1 == domain_bc)&&(bc_ced))
	    alpha = -0.25/alpha ;
	  if (1 == domain_bc) {
	    system_l0.set(index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
	    system_l1.set(index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
	    index01++ ; 
	    system_even.set(index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
	    system_odd.set(index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
	    index++ ;
	  }
	  else {
	    system_l0.set(index01, index01) = 1. ;
	    system_l1.set(index01, index01) = 1. ;
	    system_even.set(index, index) = 1. ;
	    system_odd.set(index, index) = 1. ;
	    if (derivative) {
	      system_l0.set(index01+1, index01) = (nr-1)*(nr-1)/alpha ;
	      system_l1.set(index01+1, index01) = (nr-1)*(nr-1)/alpha ;
	      system_even.set(index+1, index) = (nr-1)*(nr-1)/alpha ;
	      system_odd.set(index+1, index) = (nr-1)*(nr-1)/alpha ;
	    }		
	    
	  }
	}
	for (int lz=2; lz<=domain_bc; lz++) {
 	    nr = mgrid.get_nr(lz) ;
	    alpha = mp_aff->get_alpha()[lz] ;
	    if ((lz == domain_bc)&&(bc_ced))
		alpha = -0.25/alpha ;
	    system_l0.set(index01, index01+1) = 1. ;
	    system_l0.set(index01, index01+2) = -1. ;
	    system_l1.set(index01, index01+1) = 1. ;
	    system_l1.set(index01, index01+2) = -1. ;
	    index01++ ;
	    system_even.set(index, index+1) = 1. ;
	    system_even.set(index, index+2) = -1. ;
	    system_odd.set(index, index+1) = 1. ;
	    system_odd.set(index, index+2) = -1. ;
	    index++ ;
	    if (derivative) {
		system_l0.set(index01, index01) = -(nr-2)*(nr-2)/alpha ;
		system_l0.set(index01, index01+1) = (nr-1)*(nr-1)/alpha ;
		system_l1.set(index01, index01) = -(nr-2)*(nr-2)/alpha ;
		system_l1.set(index01, index01+1) = (nr-1)*(nr-1)/alpha ;
		index01++ ;
		system_even.set(index, index) = -(nr-2)*(nr-2)/alpha ;
		system_even.set(index, index+1) = (nr-1)*(nr-1)/alpha ;
		system_odd.set(index, index) = -(nr-2)*(nr-2)/alpha ;
		system_odd.set(index, index+1) = (nr-1)*(nr-1)/alpha ;
		index++ ;
	    }
	
	    if (lz == domain_bc) {
		system_l0.set(index01, index01-1) = c_val + c_der*(nr-2)*(nr-2)/alpha ;
		system_l0.set(index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
		system_l1.set(index01, index01-1) = c_val + c_der*(nr-2)*(nr-2)/alpha ;
		system_l1.set(index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
		index01++ ; 
		system_even.set(index, index-1) = c_val + c_der*(nr-2)*(nr-2)/alpha ;
		system_even.set(index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
		system_odd.set(index, index-1) = c_val + c_der*(nr-2)*(nr-2)/alpha ;
		system_odd.set(index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha ;
		index++ ;
	    }
	    else {
		system_l0.set(index01, index01-1) = 1. ;
		system_l0.set(index01, index01) = 1. ;
		system_l1.set(index01, index01-1) = 1. ;
		system_l1.set(index01, index01) = 1. ;
		system_even.set(index, index-1) = 1. ;
		system_even.set(index, index) = 1. ;
		system_odd.set(index, index-1) = 1. ;
		system_odd.set(index, index) = 1. ;
		if (derivative) {
		    system_l0.set(index01+1, index01-1) = (nr-2)*(nr-2)/alpha ;
		    system_l0.set(index01+1, index01) = (nr-1)*(nr-1)/alpha ;
		    system_l1.set(index01+1, index01-1) = (nr-2)*(nr-2)/alpha ;
		    system_l1.set(index01+1, index01) = (nr-1)*(nr-1)/alpha ;
		    system_even.set(index+1, index-1) = (nr-2)*(nr-2)/alpha ;
		    system_even.set(index+1, index) = (nr-1)*(nr-1)/alpha ;
		    system_odd.set(index+1, index-1) = (nr-2)*(nr-2)/alpha ;
		    system_odd.set(index+1, index) = (nr-1)*(nr-1)/alpha ;
		}		
	    }
	} // End of loop on lz

	assert(index01 == system01_size) ;
	assert(index == system_size) ;
	system_l0.set_lu() ;
	system_l1.set_lu() ;
	system_even.set_lu() ;
	system_odd.set_lu() ;
	if (par_mat != 0x0) {
	    Matrice* pmat0 = new Matrice(system_l0) ;
	    Matrice* pmat1 = new Matrice(system_l1) ;
	    Matrice* pmat2 = new Matrice(system_even) ;
	    Matrice* pmat3 = new Matrice(system_odd) ;
	    par_mat->add_matrice_mod(*pmat0, 0) ;
	    par_mat->add_matrice_mod(*pmat1, 1) ;
	    par_mat->add_matrice_mod(*pmat2, 2) ;
	    par_mat->add_matrice_mod(*pmat3, 3) ;
	}
    }// End of if (need_matrices) ...

    const Matrice& sys_l0 = (need_matrices ? system_l0 : par_mat->get_matrice_mod(0) ) ;
    const Matrice& sys_l1 = (need_matrices ? system_l1 : par_mat->get_matrice_mod(1) ) ;
    const Matrice& sys_even = (need_matrices ? system_even : 
			       par_mat->get_matrice_mod(2) ) ;
    const Matrice& sys_odd = (need_matrices ? system_odd : 
			      par_mat->get_matrice_mod(3) ) ;

    va.coef() ;
    va.ylm() ;
    const Base_val& base = get_spectral_base() ;
    Mtbl_cf& coef = *va.c_cf ;
    if (va.c != 0x0) {
	delete va.c ;
	va.c = 0x0 ;
    }
    int m_q, l_q, base_r ;
    for (int k=0; k<np+2; k++) {
	for (int j=0; j<nt; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;//#0 here as domain index
	    if ((nullite_plm(j, nt, k, np, base) == 1)&&(l_q >= l_min)) {
		l_q += dl ;
		int sys_size = (l_q < 2 ? system01_size : system_size) ;
		int nl = (l_q < 2 ? 1 : 2) ;
		Tbl rhs(sys_size) ; rhs.annule_hard() ;
		int index = 0 ;
		int pari = 1 ;
		double alpha = mp_aff->get_alpha()[0] ;
		int nr = mgrid.get_nr(0) ;
		double val_b = 0. ;
		double der_b =0. ;
		if (isexcised==false){
		  if (l_q>=2) {
		    if (base_r == R_CHEBP) {
		      double val_c = 0.; pari = 1 ;
		      for (int i=0; i<nr-nl; i++) {
			val_c += pari*coef(0, k, j, i) ;
			pari = -pari ;
		      }
		      rhs.set(index) = val_c ;
		    }
		    else {
		      assert( base_r == R_CHEBI) ;
		      double der_c = 0.; pari = 1 ;
		      for (int i=0; i<nr-nl-1; i++) {
			der_c += pari*(2*i+1)*coef(0, k, j, i) ;
			pari = -pari ;
		      }
		      rhs.set(index) = der_c / alpha ;
		    }
		    index++ ;
		  }
		  if (base_r == R_CHEBP) {
		    for (int i=0; i<nr-nl; i++) {
		      val_b += coef(0, k, j, i) ;
		      der_b += 4*i*i*coef(0, k, j, i) ;
		    }
		  }
		  else {
		    assert(base_r == R_CHEBI) ;
		    for (int i=0; i<nr-nl-1; i++) {
		      val_b += coef(0, k, j, i) ;
		      der_b += (2*i+1)*(2*i+1)*coef(0, k, j, i) ;		    
		    }
		  }
		  if (domain_bc==0) {
		    rhs.set(index) = mu(k,j) - c_val*val_b - c_der*der_b/alpha ;
		    index++ ;
		  }
		  else {
		    rhs.set(index) = -val_b ;
		    if (derivative) 
		      rhs.set(index+1) = -der_b/alpha ;
		    
		    // Now for l=1
		    alpha = mp_aff->get_alpha()[1] ;
		    if ((1 == domain_bc)&&(bc_ced))
		      alpha = -0.25/alpha ;
		    nr = mgrid.get_nr(1) ;
		    int nr_lim = nr - n_conditions ;
		    val_b = 0 ; pari = 1 ;
		    for (int i=0; i<nr_lim; i++) { 
		      val_b += pari*coef(1, k, j, i) ;
		      pari = -pari ;
		    }
		    rhs.set(index) += val_b ;
		    index++ ;
		    if (derivative) {
		      der_b = 0 ; pari = -1 ;
		      for (int i=0; i<nr_lim; i++) {
			der_b += pari*i*i*coef(1, k, j, i) ;
			pari = -pari ;
		      }
		      rhs.set(index) += der_b/alpha ;
		      index++ ;
		    }
		    val_b = 0 ;
		    for (int i=0; i<nr_lim; i++)
		      val_b += coef(1, k, j, i) ;
		    der_b = 0 ;
		    for (int i=0; i<nr_lim; i++) {
		      der_b += i*i*coef(1, k, j, i) ;
		    }
		    if (domain_bc==1) {
		      rhs.set(index) = mu(k,j) - c_val*val_b - c_der*der_b/alpha ;
		      index++ ;
		    }
		    else {
		      rhs.set(index) = -val_b ;
		      if (derivative) rhs.set(index+1) = -der_b / alpha ;
		    }
		  }
		}
		else{
		     alpha = mp_aff->get_alpha()[1] ;
		    if ((1 == domain_bc)&&(bc_ced))
			alpha = -0.25/alpha ;
		    nr = mgrid.get_nr(1) ;
		    int nr_lim = nr - 1 ;
		    val_b = 0 ;
		    for (int i=0; i<nr_lim; i++)
			val_b += coef(1, k, j, i) ;
		    der_b = 0 ;
		    for (int i=0; i<nr_lim; i++) {
			der_b += i*i*coef(1, k, j, i) ;
		    }
		    if (domain_bc==1) {
			rhs.set(index) = mu(k,j) - c_val*val_b - c_der*der_b/alpha ;
			index++ ;
		    }
		    else {
			rhs.set(index) = -val_b ;
			if (derivative) rhs.set(index+1) = -der_b / alpha ;
		    }
		}



		for (int lz=2; lz<=domain_bc; lz++) {
		    alpha = mp_aff->get_alpha()[lz] ;
		    if ((lz == domain_bc)&&(bc_ced))
			alpha = -0.25/alpha ;
		    nr = mgrid.get_nr(lz) ;
		    int nr_lim = nr - n_conditions ;
		    val_b = 0 ; pari = 1 ;
		    for (int i=0; i<nr_lim; i++) { 
			val_b += pari*coef(lz, k, j, i) ;
			pari = -pari ;
		    }
		    rhs.set(index) += val_b ;
		    index++ ;
		    if (derivative) {
			der_b = 0 ; pari = -1 ;
			for (int i=0; i<nr_lim; i++) {
			    der_b += pari*i*i*coef(lz, k, j, i) ;
			    pari = -pari ;
			}
			rhs.set(index) += der_b/alpha ;
			index++ ;
		    }
		    val_b = 0 ;
		    for (int i=0; i<nr_lim; i++)
			val_b += coef(lz, k, j, i) ;
		    der_b = 0 ;
		    for (int i=0; i<nr_lim; i++) {
			der_b += i*i*coef(lz, k, j, i) ;
		    }
		    if (domain_bc==lz) {
			rhs.set(index) = mu(k,j) - c_val*val_b - c_der*der_b/alpha ;
			index++ ;
		    }
		    else {
			rhs.set(index) = -val_b ;
			if (derivative) rhs.set(index+1) = -der_b / alpha ;
		    }
		}
		Tbl solut(sys_size);
		assert(l_q>=0) ;
		switch (l_q) {
		    case 0 :
			solut = sys_l0.inverse(rhs) ;
			break ;
		    case 1: 
			solut = sys_l1.inverse(rhs) ;
			break ;
		    default:
			solut = (l_q%2 == 0 ? sys_even.inverse(rhs) : 
				 sys_odd.inverse(rhs)) ; 
		}
		index = 0 ;
		int dec = (base_r == R_CHEBP ? 0 : 1) ;
		nr = mgrid.get_nr(0) ;
		if (isexcised==false){
		  if (l_q>=2) {
		    coef.set(0, k, j, nr-2-dec) = solut(index) ;
		    index++ ;
		  }
		  coef.set(0, k, j, nr-1-dec) = solut(index) ;
		  index++ ;
		  if (base_r == R_CHEBI)
		    coef.set(0, k, j, nr-1) = 0 ;
		  if (domain_bc > 0) {
		    //Pour domaine 1
		    nr = mgrid.get_nr(1) ;
		    for (nl=1; nl<=n_conditions; nl++) {
		      int ii = n_conditions - nl + 1 ;
		      coef.set(1, k, j, nr-ii) = solut(index) ;
		      index++ ;
		    }
		  }
		}
		else{
		  coef.set(1,k,j,nr-1)=solut(index);
		  index++;
		}
		for (int lz=2; lz<=domain_bc; lz++) {
		  nr = mgrid.get_nr(lz) ;
		  for (nl=1; nl<=n_conditions; nl++) {
		    int ii = n_conditions - nl + 1 ;
		    coef.set(lz, k, j, nr-ii) = solut(index) ;
		    index++ ;
		  }
		}
	    } //End of nullite_plm
	    else {
	      for (int lz=0; lz<=domain_bc; lz++) 
		for (int i=0; i<mgrid.get_nr(lz); i++) 
		  coef.set(lz, k, j, i) = 0 ;
	    }
	} //End of loop on j
    } //End of loop on k
}

void Scalar::match_tau_star(Param& par_bc, Param* par_mat) {

    const Map_star* mp_star = dynamic_cast<const Map_star*>(mp) ;
    assert(mp_star != 0x0) ;

    const Mg3d& mgrid = *mp_star->get_mg() ;
    assert(mgrid.get_type_r(0) == RARE)  ;
    assert (par_bc.get_n_tbl_mod() > 0) ;
    Tbl mu = par_bc.get_tbl_mod(0) ;
    if (etat == ETATZERO) {
	if (mu.get_etat() == ETATZERO) 
	    return ;
	else
	    annule_hard() ;
    }

    int nz = mgrid.get_nzone() ;
    int nzm1 = nz - 1 ;
    bool ced = (mgrid.get_type_r(nzm1) == UNSURR) ;
    assert(par_bc.get_n_int() >= 2);
    int domain_bc = par_bc.get_int(0) ;
    bool bc_ced = ((ced) && (domain_bc == nzm1)) ;

    int n_conditions = par_bc.get_int(1) ;
    assert ((n_conditions==1)||(n_conditions==2)) ;
    bool derivative = (n_conditions == 2) ;
    int dl = 0 ; int l_min = 0 ; int excised_i =0;
    if (par_bc.get_n_int() > 2) {
	dl = par_bc.get_int(2) ;
	l_min = par_bc.get_int(3) ;
	if(par_bc.get_n_int() >4){
	excised_i = par_bc.get_int(4);
	}
    }
    bool isexcised = (excised_i==1);
    
    int nt = mgrid.get_nt(0) ;
    int np = mgrid.get_np(0) ;
    assert (par_bc.get_n_double_mod() == 2) ;
    double c_val = par_bc.get_double_mod(0) ;
    double c_der = par_bc.get_double_mod(1) ;
    if (bc_ced) {
	c_val = 1 ;
	c_der = 0 ;
	mu.annule_hard() ;
    }
    if (mu.get_etat() == ETATZERO)
	mu.annule_hard() ;
    int system01_size =0;
    int system_size =0;
    if (isexcised ==false){
    system01_size = 1 ;
    system_size = 2 ;
    }
    else{
      system01_size = -1;
      system_size = -1;
    }
    for (int lz=1; lz<=domain_bc; lz++) {
	    system01_size += n_conditions ;
	    system_size += n_conditions ;
    }
    assert (etat != ETATNONDEF) ;

    bool need_matrices = true ;
    if (par_mat != 0x0)
	if (par_mat->get_n_matrice_mod() == 4) 
	    need_matrices = false ;
    
	// int *sizes_01, *sizes ; sizes01 = new double[4] ; sizes = new double[4] ;
	// sizes[0] = np ; sizes[1] = nt ; sizes[2] = system_size ; sizes[3] = system_size ;
	// sizes_01[0] = np ; sizes_01[1] = nt ; sizes_01[2] = system01_size ; sizes_01[3] = system01_size ;

	// Dim_tbl dim_01(4,sizes_01) ; Dim_tbl dimdim(4,sizes) ;

	// ITbl itbl_01(dim_01) ; Itbl itbl_dim(dimdim) ;

	int npp2 = np ;
	Mg3d mg_01(npp2, system01_size, system01_size, nt, SYM, SYM) ; Mg3d mgmg(npp2, system_size, system_size, nt, SYM, SYM) ;

    Valeur system_l0(mg_01) ;
    Valeur system_l1(mg_01) ;
    Valeur system_even(mgmg) ;
    Valeur system_odd(mgmg) ;  


    
    if (need_matrices) {
	system_l0.annule_hard() ;
	system_l1.annule_hard() ;
	system_even.annule_hard() ;
	system_odd.annule_hard() ;
	
	int nr = mgrid.get_nr(0);
	Valeur alpha = mp_star->get_alpha() ; 
	for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++){
		int index = 0 ; int index01 = 0 ;
	if (isexcised == false){
	  system_even.set(k,j,index, index) = 1. ;
	  system_even.set(k,j,index, index + 1) = -1. ; //C_{N-1} - C_N = \sum (-1)^i C_i
	  system_odd.set(k,j,index, index) = -(2.*nr - 5.)/alpha(0,k,j,0) ; 
	  system_odd.set(k,j,index, index+1) = (2.*nr - 3.)/alpha(0,k,j,0) ;
	  index++ ;
	  if (domain_bc == 0) {
	    system_l0.set(k,j,index01, index01) = c_val + c_der*4.*(nr-1)*(nr-1)/alpha(0,k,j,0) ;
	    system_l1.set(k,j,index01, index01) = c_val + c_der*(2*nr-3)*(2*nr-3)/alpha(0,k,j,0) ;
	    index01++ ;
	    system_even.set(k,j,index, index-1) = c_val + c_der*4.*(nr-2)*(nr-2)/alpha(0,k,j,0) ;
	    system_even.set(k,j,index, index) = c_val + c_der*4.*(nr-1)*(nr-1)/alpha(0,k,j,0) ;
	    system_odd.set(k,j,index, index-1) = c_val + c_der*(2*nr-5)*(2*nr-5)/alpha(0,k,j,0) ;
	    system_odd.set(k,j,index, index) = c_val + c_der*(2*nr-3)*(2*nr-3)/alpha(0,k,j,0) ;
	    index++ ;
	  }
	  else {
	    system_l0.set(k,j,index01, index01) = 1. ;
	    system_l1.set(k,j,index01, index01) = 1. ;
	    system_even.set(k,j,index, index-1) = 1. ;
	    system_even.set(k,j,index, index) = 1. ;
	    system_odd.set(k,j,index, index-1) = 1. ;
	    system_odd.set(k,j,index, index) = 1. ;
	    if (derivative) {
	      system_l0.set(k,j,index01+1, index01) = 4*(nr-1)*(nr-1)/alpha(0,k,j,0) ;
	      system_l1.set(k,j,index01+1, index01) = (2*nr-3)*(2*nr-3)/alpha(0,k,j,0) ;
	      system_even.set(k,j,index+1, index-1) = 4*(nr-2)*(nr-2)/alpha(0,k,j,0) ;
	      system_even.set(k,j,index+1, index) = 4*(nr-1)*(nr-1)/alpha(0,k,j,0) ;
	      system_odd.set(k,j,index+1, index-1) = (2*nr-5)*(2*nr-5)/alpha(0,k,j,0) ;
	      system_odd.set(k,j,index+1, index) = (2*nr-3)*(2*nr-3)/alpha(0,k,j,0) ;
	    }
	  }
	  if (domain_bc > 0) {
	    //	Do things for lz=1;
	    nr = mgrid.get_nr(1) ;
	    // alpha = mp_star->get_alpha() ;
	    // if ((1 == domain_bc)&&(bc_ced))
		// alpha = -0.25/alpha ; // No ced in Map_star
	    system_l0.set(k,j,index01, index01+1) = 1. ;
	    system_l0.set(k,j,index01, index01+2) = -1. ;
	    system_l1.set(k,j,index01, index01+1) = 1. ;
	    system_l1.set(k,j,index01, index01+2) = -1. ;
	    index01++ ;
	    system_even.set(k,j,index, index+1) = 1. ;
	    system_even.set(k,j,index, index+2) = -1. ;
	    system_odd.set(k,j,index, index+1) = 1. ;
	    system_odd.set(k,j,index, index+2) = -1. ;
	    index++ ;
	    if (derivative) {
		system_l0.set(k,j,index01, index01) = -(nr-2)*(nr-2)/alpha(1,k,j,0) ;
		system_l0.set(k,j,index01, index01+1) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
		system_l1.set(k,j,index01, index01) = -(nr-2)*(nr-2)/alpha(1,k,j,0) ;
		system_l1.set(k,j,index01, index01+1) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
		index01++ ;
		system_even.set(k,j,index, index) = -(nr-2)*(nr-2)/alpha(1,k,j,0) ;
		system_even.set(k,j,index, index+1) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
		system_odd.set(k,j,index, index) = -(nr-2)*(nr-2)/alpha(1,k,j,0) ;
		system_odd.set(k,j,index, index+1) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
		index++ ;
	    }
	
	    if (1 == domain_bc) {
		system_l0.set(k,j,index01, index01-1) = c_val + c_der*(nr-2)*(nr-2)/alpha(1,k,j,0) ;
		system_l0.set(k,j,index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha(1,k,j,0) ;
		system_l1.set(k,j,index01, index01-1) = c_val + c_der*(nr-2)*(nr-2)/alpha(1,k,j,0) ;
		system_l1.set(k,j,index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha(1,k,j,0) ;
		index01++ ; 
		system_even.set(k,j,index, index-1) = c_val + c_der*(nr-2)*(nr-2)/alpha(1,k,j,0) ;
		system_even.set(k,j,index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha(1,k,j,0) ;
		system_odd.set(k,j,index, index-1) = c_val + c_der*(nr-2)*(nr-2)/alpha(1,k,j,0) ;
		system_odd.set(k,j,index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha(1,k,j,0) ;
		index++ ;
	    }
	    else {
		system_l0.set(k,j,index01, index01-1) = 1. ;
		system_l0.set(k,j,index01, index01) = 1. ;
		system_l1.set(k,j,index01, index01-1) = 1. ;
		system_l1.set(k,j,index01, index01) = 1. ;
		system_even.set(k,j,index, index-1) = 1. ;
		system_even.set(k,j,index, index) = 1. ;
		system_odd.set(k,j,index, index-1) = 1. ;
		system_odd.set(k,j,index, index) = 1. ;
		if (derivative) {
		    system_l0.set(k,j,index01+1, index01-1) = (nr-2)*(nr-2)/alpha(1,k,j,0) ;
		    system_l0.set(k,j,index01+1, index01) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
		    system_l1.set(k,j,index01+1, index01-1) = (nr-2)*(nr-2)/alpha(1,k,j,0) ;
		    system_l1.set(k,j,index01+1, index01) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
		    system_even.set(k,j,index+1, index-1) = (nr-2)*(nr-2)/alpha(1,k,j,0) ;
		    system_even.set(k,j,index+1, index) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
		    system_odd.set(k,j,index+1, index-1) = (nr-2)*(nr-2)/alpha(1,k,j,0) ;
		    system_odd.set(k,j,index+1, index) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
		}		
	    }
	  }
	}
	else {
	  nr = mgrid.get_nr(1) ;
	  alpha = mp_star->get_alpha() ;
	//   if ((1 == domain_bc)&&(bc_ced))
	//     alpha = -0.25/alpha ; // No ced in Map_star
	  if (1 == domain_bc) {
	    system_l0.set(k,j,index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha(1,k,j,0) ;
	    system_l1.set(k,j,index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha(1,k,j,0) ;
	    index01++ ; 
	    system_even.set(k,j,index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha(1,k,j,0) ;
	    system_odd.set(k,j,index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha(1,k,j,0) ;
	    index++ ;
	  }
	  else {
	    system_l0.set(k,j,index01, index01) = 1. ;
	    system_l1.set(k,j,index01, index01) = 1. ;
	    system_even.set(k,j,index, index) = 1. ;
	    system_odd.set(k,j,index, index) = 1. ;
	    if (derivative) {
	      system_l0.set(k,j,index01+1, index01) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
	      system_l1.set(k,j,index01+1, index01) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
	      system_even.set(k,j,index+1, index) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
	      system_odd.set(k,j,index+1, index) = (nr-1)*(nr-1)/alpha(1,k,j,0) ;
	    }		
	    
	  }
	}
	for (int lz=2; lz<=domain_bc; lz++) {
 	    nr = mgrid.get_nr(lz) ;
	    alpha = mp_star->get_alpha() ;
	    // if ((lz == domain_bc)&&(bc_ced))
		// alpha = -0.25/alpha ; // No ced in Map_star
	    system_l0.set(k,j,index01, index01+1) = 1. ;
	    system_l0.set(k,j,index01, index01+2) = -1. ;
	    system_l1.set(k,j,index01, index01+1) = 1. ;
	    system_l1.set(k,j,index01, index01+2) = -1. ;
	    index01++ ;
	    system_even.set(k,j,index, index+1) = 1. ;
	    system_even.set(k,j,index, index+2) = -1. ;
	    system_odd.set(k,j,index, index+1) = 1. ;
	    system_odd.set(k,j,index, index+2) = -1. ;
	    index++ ;
	    if (derivative) {
		system_l0.set(k,j,index01, index01) = -(nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		system_l0.set(k,j,index01, index01+1) = (nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		system_l1.set(k,j,index01, index01) = -(nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		system_l1.set(k,j,index01, index01+1) = (nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		index01++ ;
		system_even.set(k,j,index, index) = -(nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		system_even.set(k,j,index, index+1) = (nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		system_odd.set(k,j,index, index) = -(nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		system_odd.set(k,j,index, index+1) = (nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		index++ ;
	    }
	
	    if (lz == domain_bc) {
		system_l0.set(k,j,index01, index01-1) = c_val + c_der*(nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		system_l0.set(k,j,index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		system_l1.set(k,j,index01, index01-1) = c_val + c_der*(nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		system_l1.set(k,j,index01, index01) = c_val + c_der*(nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		index01++ ; 
		system_even.set(k,j,index, index-1) = c_val + c_der*(nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		system_even.set(k,j,index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		system_odd.set(k,j,index, index-1) = c_val + c_der*(nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		system_odd.set(k,j,index, index) = c_val + c_der*(nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		index++ ;
	    }
	    else {
		system_l0.set(k,j,index01, index01-1) = 1. ;
		system_l0.set(k,j,index01, index01) = 1. ;
		system_l1.set(k,j,index01, index01-1) = 1. ;
		system_l1.set(k,j,index01, index01) = 1. ;
		system_even.set(k,j,index, index-1) = 1. ;
		system_even.set(k,j,index, index) = 1. ;
		system_odd.set(k,j,index, index-1) = 1. ;
		system_odd.set(k,j,index, index) = 1. ;
		if (derivative) {
		    system_l0.set(k,j,index01+1, index01-1) = (nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		    system_l0.set(k,j,index01+1, index01) = (nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		    system_l1.set(k,j,index01+1, index01-1) = (nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		    system_l1.set(k,j,index01+1, index01) = (nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		    system_even.set(k,j,index+1, index-1) = (nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		    system_even.set(k,j,index+1, index) = (nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		    system_odd.set(k,j,index+1, index-1) = (nr-2)*(nr-2)/alpha(lz,k,j,0) ;
		    system_odd.set(k,j,index+1, index) = (nr-1)*(nr-1)/alpha(lz,k,j,0) ;
		}		
	    }
	} // End of loop on lz

	assert(index01 == system01_size) ;
	assert(index == system_size) ;
	// system_l0.set_lu() ;
	// system_l1.set_lu() ;
	// system_even.set_lu() ;
	// system_odd.set_lu() ;
	// if (par_mat != 0x0) {
	//     Matrice* pmat0 = new Matrice(system_l0) ;
	//     Matrice* pmat1 = new Matrice(system_l1) ;
	//     Matrice* pmat2 = new Matrice(system_even) ;
	//     Matrice* pmat3 = new Matrice(system_odd) ;
	//     par_mat->add_matrice_mod(*pmat0, 0) ;
	//     par_mat->add_matrice_mod(*pmat1, 1) ;
	//     par_mat->add_matrice_mod(*pmat2, 2) ;
	//     par_mat->add_matrice_mod(*pmat3, 3) ;
	// }
	}
    }// End of if (need_matrices) ...

    const Valeur& sys_l0 = system_l0 ; //(need_matrices ? system_l0 : par_mat->get_matrice_mod(0) ) ;
    const Valeur& sys_l1 = system_l1 ; //(need_matrices ? system_l1 : par_mat->get_matrice_mod(1) ) ;
    const Valeur& sys_even = system_even ; //(need_matrices ? system_even : 
			       //par_mat->get_matrice_mod(2) ) ;
    const Valeur& sys_odd = system_odd ; //(need_matrices ? system_odd : 
			      //par_mat->get_matrice_mod(3) ) ;

    va.coef() ;
    va.ylm() ;
    const Base_val& base = get_spectral_base() ;
    Mtbl_cf& coef = *va.c_cf ; cout << "COEF : " << coef << endl;
    if (va.c != 0x0) {
	delete va.c ;
	va.c = 0x0 ;
    }
    int m_q, l_q, base_r ;
    for (int k=0; k<np+2; k++) {
	for (int j=0; j<nt; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;//#0 here as domain index
	    if ((nullite_plm(j, nt, k, np, base) == 1)&&(l_q >= l_min)) {
		l_q += dl ;
		int sys_size = (l_q < 2 ? system01_size : system_size) ;
		int nl = (l_q < 2 ? 1 : 2) ;
		Tbl rhs(np+2,nt,sys_size) ; rhs.annule_hard() ;
		int index = 0 ;
		int pari = 1 ;
		Valeur alpha = mp_star->get_alpha() ;
		int nr = mgrid.get_nr(0) ;
		double val_b = 0. ;
		double der_b =0. ;
		if (isexcised==false){
		  if (l_q>=2) {
		    if (base_r == R_CHEBP) {
		      double val_c = 0.; pari = 1 ;
		      for (int i=0; i<nr-nl; i++) {
			val_c += pari*coef(0, k, j, i) ;
			pari = -pari ;
		      }
		      rhs.set(k,j,index) = val_c ;
		    }
		    else {
		      assert( base_r == R_CHEBI) ;
		      double der_c = 0.; pari = 1 ;
		      for (int i=0; i<nr-nl-1; i++) {
			der_c += pari*(2*i+1)*coef(0, k, j, i) ;
			pari = -pari ;
		      }
		      rhs.set(k,j,index) = der_c / alpha(0, k, j, 0) ;
		    }
		    index++ ;
		  }
		  if (base_r == R_CHEBP) {
		    for (int i=0; i<nr-nl; i++) {
		      val_b += coef(0, k, j, i) ;
		      der_b += 4*i*i*coef(0, k, j, i) ;
		    }
		  }
		  else {
		    assert(base_r == R_CHEBI) ;
		    for (int i=0; i<nr-nl-1; i++) {
		      val_b += coef(0, k, j, i) ;
		      der_b += (2*i+1)*(2*i+1)*coef(0, k, j, i) ;		    
		    }
		  }
		  if (domain_bc==0) {
		    rhs.set(k,j,index) = mu(k,j) - c_val*val_b - c_der*der_b/alpha(0, k, j, 0) ;
		    index++ ;
		  }
		  else {
		    rhs.set(k,j,index) = -val_b ;
		    if (derivative) 
		      rhs.set(k,j,index+1) = -der_b/alpha(0, k, j, 0) ;
		    
		    // Now for l=1
		    alpha = mp_star->get_alpha() ;
		    // if ((1 == domain_bc)&&(bc_ced))
		    //   alpha = -0.25/alpha ; // No ced in Map_star
		    nr = mgrid.get_nr(1) ;
		    int nr_lim = nr - n_conditions ;
		    val_b = 0 ; pari = 1 ;
		    for (int i=0; i<nr_lim; i++) { 
		      val_b += pari*coef(1, k, j, i) ;
		      pari = -pari ;
		    }
		    rhs.set(k,j,index) += val_b ;
		    index++ ;
		    if (derivative) {
		      der_b = 0 ; pari = -1 ;
		      for (int i=0; i<nr_lim; i++) {
			der_b += pari*i*i*coef(1, k, j, i) ;
			pari = -pari ;
		      }
		      rhs.set(k,j,index) += der_b/alpha(1,k,j,0) ;
		      index++ ;
		    }
		    val_b = 0 ;
		    for (int i=0; i<nr_lim; i++)
		      val_b += coef(1, k, j, i) ;
		    der_b = 0 ;
		    for (int i=0; i<nr_lim; i++) {
		      der_b += i*i*coef(1, k, j, i) ;
		    }
		    if (domain_bc==1) {
		      rhs.set(k,j,index) = mu(k,j) - c_val*val_b - c_der*der_b/alpha(1,k,j,0) ;
		      index++ ;
		    }
		    else {
		      rhs.set(k,j,index) = -val_b ;
		      if (derivative) rhs.set(k,j,index+1) = -der_b / alpha(1,k,j,0) ;
		    }
		  }
		}
		else{
		     alpha = mp_star->get_alpha() ;
		    // if ((1 == domain_bc)&&(bc_ced))
			// alpha = -0.25/alpha ; // No ced in Map_star
		    nr = mgrid.get_nr(1) ;
		    int nr_lim = nr - 1 ;
		    val_b = 0 ;
		    for (int i=0; i<nr_lim; i++)
			val_b += coef(1, k, j, i) ;
		    der_b = 0 ;
		    for (int i=0; i<nr_lim; i++) {
			der_b += i*i*coef(1, k, j, i) ;
		    }
		    if (domain_bc==1) {
			rhs.set(k,j,index) = mu(k,j) - c_val*val_b - c_der*der_b/alpha(1,k,j,0) ;
			index++ ;
		    }
		    else {
			rhs.set(k,j,index) = -val_b ;
			if (derivative) rhs.set(k,j,index+1) = -der_b / alpha(1,k,j,0) ;
		    }
		}



		for (int lz=2; lz<=domain_bc; lz++) {
		    alpha = mp_star->get_alpha() ;
		    // if ((lz == domain_bc)&&(bc_ced))
			// alpha = -0.25/alpha ; No ced in Map_star
		    nr = mgrid.get_nr(lz) ;
		    int nr_lim = nr - n_conditions ;
		    val_b = 0 ; pari = 1 ;
		    for (int i=0; i<nr_lim; i++) { 
			val_b += pari*coef(lz, k, j, i) ;
			pari = -pari ;
		    }
		    rhs.set(k,j,index) += val_b ;
		    index++ ;
		    if (derivative) {
			der_b = 0 ; pari = -1 ;
			for (int i=0; i<nr_lim; i++) {
			    der_b += pari*i*i*coef(lz, k, j, i) ;
			    pari = -pari ;
			}
			rhs.set(k,j,index) += der_b/alpha(lz, k, j, 0) ;
			index++ ;
		    }
		    val_b = 0 ;
		    for (int i=0; i<nr_lim; i++)
			val_b += coef(lz, k, j, i) ;
		    der_b = 0 ;
		    for (int i=0; i<nr_lim; i++) {
			der_b += i*i*coef(lz, k, j, i) ;
		    }
		    if (domain_bc==lz) {
			rhs.set(k,j,index) = mu(k,j) - c_val*val_b - c_der*der_b/alpha(lz, k, j, 0) ;
			index++ ;
		    }
		    else {
			rhs.set(k,j,index) = -val_b ;
			if (derivative) rhs.set(k,j,index+1) = -der_b / alpha(lz, k, j, 0) ;
		    }
		}
		Tbl solut(sys_size);
		assert(l_q>=0) ;
		Tbl rhs_tmp(sys_size) ; rhs_tmp.set_etat_qcq() ;
		switch (l_q) {
		    case 0 :{
			Matrice sys_l0_tmp(sys_size, sys_size) ; sys_l0_tmp.set_etat_qcq() ;
			for (int ind1=0; ind1<sys_size; ind1++){
				rhs_tmp.set(ind1) = rhs(k,j,ind1) ;
				for (int ind2=0; ind2<sys_size; ind2++)
					sys_l0_tmp.set(ind1,ind2) = sys_l0(k,j,ind1,ind2) ;
			}
			sys_l0_tmp.set_lu() ;
			solut = sys_l0_tmp.inverse(rhs_tmp) ;
			}
			break ;
		    case 1:{
			Matrice sys_l1_tmp(sys_size, sys_size) ; sys_l1_tmp.set_etat_qcq() ;
			for (int ind1=0; ind1<sys_size; ind1++){
				rhs_tmp.set(ind1) = rhs(k,j,ind1) ;
				for (int ind2=0; ind2<sys_size; ind2++)
					sys_l1_tmp.set(ind1,ind2) = sys_l1(k,j,ind1,ind2) ;
			}
			sys_l1_tmp.set_lu() ;
			solut = sys_l1_tmp.inverse(rhs_tmp) ;
			}
			break ;
		    default:
			{
			if (l_q%2 == 0){
				Matrice sys_even_tmp(sys_size, sys_size) ; sys_even_tmp.set_etat_qcq() ;
				for (int ind1=0; ind1<sys_size; ind1++){
					rhs_tmp.set(ind1) = rhs(k,j,ind1) ;
					for (int ind2=0; ind2<sys_size; ind2++)
						sys_even_tmp.set(ind1,ind2) = sys_even(k,j,ind1,ind2) ;
				}
				sys_even_tmp.set_lu() ;
				solut = sys_even_tmp.inverse(rhs_tmp) ;
			}
			else{
				Matrice sys_odd_tmp(sys_size, sys_size) ; sys_odd_tmp.set_etat_qcq() ;
				for (int ind1=0; ind1<sys_size; ind1++){
					rhs_tmp.set(ind1) = rhs(k,j,ind1) ;
					for (int ind2=0; ind2<sys_size; ind2++)
						sys_odd_tmp.set(ind1,ind2) = sys_odd(k,j,ind1,ind2) ;
				}
				sys_odd_tmp.set_lu() ;
				solut = sys_odd_tmp.inverse(rhs_tmp) ;
			} 
			}
		}
		index = 0 ;
		int dec = (base_r == R_CHEBP ? 0 : 1) ;
		nr = mgrid.get_nr(0) ;
		if (isexcised==false){
		  if (l_q>=2) {
		    coef.set(0, k, j, nr-2-dec) = solut(index) ;
		    index++ ;
		  }
		  coef.set(0, k, j, nr-1-dec) = solut(index) ;
		  index++ ;
		  if (base_r == R_CHEBI)
		    coef.set(0, k, j, nr-1) = 0 ;
		  if (domain_bc > 0) {
		    //Pour domaine 1
		    nr = mgrid.get_nr(1) ;
		    for (nl=1; nl<=n_conditions; nl++) {
		      int ii = n_conditions - nl + 1 ;
		      coef.set(1, k, j, nr-ii) = solut(index) ;
		      index++ ;
		    }
		  }
		}
		else{
		  coef.set(1,k,j,nr-1)=solut(index);
		  index++;
		}
		for (int lz=2; lz<=domain_bc; lz++) {
		  nr = mgrid.get_nr(lz) ;
		  for (nl=1; nl<=n_conditions; nl++) {
		    int ii = n_conditions - nl + 1 ;
		    coef.set(lz, k, j, nr-ii) = solut(index) ;
		    index++ ;
		  }
		}
	    } //End of nullite_plm
	    else {
	      for (int lz=0; lz<=domain_bc; lz++) 
		for (int i=0; i<mgrid.get_nr(lz); i++) 
		  coef.set(lz, k, j, i) = 0 ;
	    }
	} //End of loop on j
    } //End of loop on k
}

}
