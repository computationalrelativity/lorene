/*
 *  Method Map_af::primr
 *
 *    (see file map.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon
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

char map_af_primr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/13 08:53:03  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:13:12  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2013/04/25 15:46:05  j_novak
 * Added special treatment in the case np = 1, for type_p = NONSYM.
 *
 * Revision 1.4  2007/12/20 09:11:05  jl_cornou
 * Correction of an error in op_sxpun about Jacobi(0,2) polynomials
 *
 * Revision 1.3  2004/07/26 16:02:23  j_novak
 * Added a flag to specify whether the primitive should be zero either at r=0
 * or at r going to infinity.
 *
 * Revision 1.1  2004/06/14 15:25:34  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */


// C headers
#include <cstdlib>

// Lorene headers
#include "map.h"
#include "tensor.h"

namespace Lorene {
void _primr_pas_prevu(const Tbl&, int, const Tbl&, Tbl&, int&, Tbl&) ; 
void _primr_r_cheb(const Tbl&, int, const Tbl&, Tbl&, int&, Tbl&) ; 
void _primr_r_chebp(const Tbl&, int, const Tbl&, Tbl&, int&, Tbl&) ; 
void _primr_r_chebi(const Tbl&, int, const Tbl&, Tbl&, int&, Tbl&) ; 
void _primr_r_chebpim_p(const Tbl&, int, const Tbl&, Tbl&, int&, Tbl&) ; 
void _primr_r_chebpim_i(const Tbl&, int, const Tbl&, Tbl&, int&, Tbl&) ; 
void _primr_r_jaco02(const Tbl&, int, const Tbl&, Tbl&, int&, Tbl&) ;

void Map_af::primr(const Scalar& uu, Scalar& resu, bool null_infty) const {

    static void (*prim_domain[MAX_BASE])(const Tbl&, int bin, const Tbl&, 
        Tbl&, int&, Tbl& ) ; 
    static bool first_call = true ; 

    // Initialisation at first call of the array of primitive functions 
    // depending upon the basis in r
    // ----------------------------------------------------------------
    if (first_call) {
	    for (int i=0 ; i<MAX_BASE ; i++) prim_domain[i] = _primr_pas_prevu ;
        
	    prim_domain[R_CHEB >> TRA_R] = _primr_r_cheb ;
	    prim_domain[R_CHEBU >> TRA_R] = _primr_r_cheb ;
	    prim_domain[R_CHEBP >> TRA_R] = _primr_r_chebp ;
	    prim_domain[R_CHEBI >> TRA_R] = _primr_r_chebi ;
	    prim_domain[R_CHEBPIM_P >> TRA_R] = _primr_r_chebpim_p ;
	    prim_domain[R_CHEBPIM_I >> TRA_R] = _primr_r_chebpim_i ;
	    prim_domain[R_JACO02 >> TRA_R] = _primr_r_jaco02 ;

        first_call = false ; 
    }

    // End of first call operations
    // ----------------------------
    
    assert(uu.get_etat() != ETATNONDEF) ;
    assert(uu.get_mp().get_mg() == mg) ;  
    assert(resu.get_mp().get_mg() == mg) ;  
    
    // Special case of vanishing input:
    if (uu.get_etat() == ETATZERO) {
	    resu.set_etat_zero() ; 
        return ; 
    }

    // General case
    assert( (uu.get_etat() == ETATQCQ) || (uu.get_etat() == ETATUN) ) ; 
    assert(uu.check_dzpuis(2)) ; 

    int nz = mg->get_nzone() ; 
    int nzm1 = nz - 1 ; 
    int np = mg->get_np(0) ;    
    int nt = mg->get_nt(0) ;
#ifndef NDEBUG
	for (int l=1; l<nz; l++) {     
	  assert (mg->get_np(l) == np) ;
	  assert (mg->get_nt(l) == nt) ;
	}
#endif

    const Valeur& vuu = uu.get_spectral_va() ; 
    vuu.coef() ; 

    const Mtbl_cf& cuu = *(vuu.c_cf) ; 
    assert(cuu.t != 0x0) ; 

    const Base_val& buu = vuu.get_base() ; // spectral bases of the input
    
    resu.set_etat_qcq() ;   // result in ordinary state
    Valeur& vprim = resu.set_spectral_va() ; 
    vprim.set_etat_cf_qcq() ;  // allocates the Mtbl_cf for the coefficients
                               // of the result
    Mtbl_cf& cprim = *(vprim.c_cf) ;   
    cprim.set_etat_qcq() ;  // allocates the Tbl's to store the coefficients
                            // of the result in each domain   
    
    Base_val& bprim = cprim.base ;    // spectral bases of the result
    
    Tbl val_rmin(np+2,nt) ;  // Values of primitive at the left boundary 
                             // in the current domain 
    Tbl val_rmax(np+2,nt) ;  // same but for the right boundary
    
    val_rmin.set_etat_zero() ;  // initialisation: primitive = 0 at r=0
    
    int lmax = (mg->get_type_r(nzm1) == UNSURR) ? nz-2 : nzm1 ;  
    
    for (int l=0; l<=lmax; l++) {
        assert(cuu.t[l] != 0x0) ; 
        assert(cprim.t[l] != 0x0) ; 
        const Tbl& cfuu = *(cuu.t[l]) ; 
        Tbl& cfprim = *(cprim.t[l]) ; 
	
        int buu_dom = buu.get_b(l) ; 
        int base_r = (buu_dom & MSQ_R) >> TRA_R ;
        
        prim_domain[base_r](cfuu, buu_dom, val_rmin, cfprim, bprim.b[l],
            val_rmax) ; 
            
        cfprim *= alpha[l] ; 
        val_rmin = alpha[l] * val_rmax / alpha[l+1] ;  // for next domain
    }     
    
    // Special case of compactified external domain (CED)
    // --------------------------------------------------
    if (mg->get_type_r(nzm1) == UNSURR) {
        val_rmin = - val_rmin ; 
        const Tbl& cfuu = *(cuu.t[nzm1]) ; 
        Tbl& cfprim = *(cprim.t[nzm1]) ; 
	
        int buu_dom = buu.get_b(nzm1) ; 
        int base_r = (buu_dom & MSQ_R) >> TRA_R ;
        assert(base_r == R_CHEBU) ; 
        
        prim_domain[base_r](cfuu, buu_dom, val_rmin, cfprim, bprim.b[nzm1],
            val_rmax) ;
        
        cfprim *= - alpha[nzm1] ;  
    }
    
    if (null_infty) 
      for (int k=0; k<np; k++)  //## not very elegant!
	for(int j=0; j<nt; j++) 
	  val_rmax.set(k,j) = cprim.val_out_bound_jk(nzm1, j, k) ;
    
    // The output spectral bases (set on the Mtbl_cf) are copied to the Valeur:
    vprim.set_base(bprim) ; 

    if (null_infty)
      for (int l=0; l<nz; l++) //## not very elegant!
	for (int k=0; k<np; k++) 
	  for(int j=0; j<nt; j++) 
	    for (int i=0; i<mg->get_nr(l); i++) 
	      vprim.set(l, k, j, i) -= val_rmax(k,j) ;
	 
    
}
}
