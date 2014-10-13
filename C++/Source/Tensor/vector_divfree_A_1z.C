/*
 *  Methods to impose the Dirac gauge: divergence-free condition.
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2006  Jerome Novak
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

char sol_Dirac_A_1z_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:53:45  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:20  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2009/10/23 13:18:46  j_novak
 * Minor modifications
 *
 * Revision 1.1  2008/08/27 09:01:27  jl_cornou
 * Methods for solving Dirac systems for divergence free vectors
 *
 *
 * $Header$
 *
 */


// C headers
#include <cstdlib>
#include <cassert>
#include <cmath>

// Lorene headers
#include "metric.h"
#include "diff.h"
#include "proto.h"
#include "param.h"

//----------------------------------------------------------------------------------
//
//                               sol_Dirac_A
//				1 seule zone !
//----------------------------------------------------------------------------------

namespace Lorene {
void Vector_divfree::sol_Dirac_A_1z(const Scalar& aaa, Scalar& tilde_vr, Scalar& tilde_eta,
				   const Param* par_bc) const {

    const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
    assert(mp_aff != 0x0) ; //Only affine mapping for the moment

    const Mg3d& mgrid = *mp_aff->get_mg() ;
    assert(mgrid.get_type_r(0) == RARE)  ;
    if (aaa.get_etat() == ETATZERO) {
	tilde_vr = 0 ;
	tilde_eta = 0 ;
	return ;
    }
    assert(aaa.get_etat() != ETATNONDEF) ;
    int nz = mgrid.get_nzone() ;
    int nzm1 = nz - 1 ;
    bool ced = (mgrid.get_type_r(nzm1) == UNSURR) ;
    int n_shell = ced ? nz-2 : nzm1 ;
    int nz_bc = nzm1 ;
    if (par_bc != 0x0)
	if (par_bc->get_n_int() > 0) nz_bc = par_bc->get_int() ;
    n_shell = (nz_bc < n_shell ? nz_bc : n_shell) ;
//#ifndef NDEBUG
//    if (!cedbc) {
//	assert(par_bc != 0x0) ;
//	assert(par_bc->get_n_tbl_mod() >= 3) ;
//    }
//#endif
    int nt = mgrid.get_nt(0) ;
    int np = mgrid.get_np(0) ;

    Scalar source = aaa ;
    Scalar source_coq = aaa ;
    source_coq.annule_domain(0) ;
    if (ced) source_coq.annule_domain(nzm1) ;
    source_coq.mult_r() ;
    source.set_spectral_va().ylm() ;
    source_coq.set_spectral_va().ylm() ;
    Base_val base = source.get_spectral_base() ;
    base.mult_x() ;

    tilde_vr.annule_hard() ;
    tilde_vr.set_spectral_base(base) ;
    tilde_vr.set_spectral_va().set_etat_cf_qcq() ;
    tilde_vr.set_spectral_va().c_cf->annule_hard() ;   
    tilde_eta.annule_hard() ;
    tilde_eta.set_spectral_base(base) ;
    tilde_eta.set_spectral_va().set_etat_cf_qcq() ;
    tilde_eta.set_spectral_va().c_cf->annule_hard() ;   
 
    Mtbl_cf sol_vr(mgrid, base) ; sol_vr.annule_hard() ;
    Mtbl_cf sol_eta(mgrid, base) ; sol_eta.annule_hard() ;
    
    int l_q, m_q, base_r ;

    //---------------
    //--  NUCLEUS ---
    //---------------
    {int lz = 0 ;  
    int nr = mgrid.get_nr(lz) ;
    double alpha = mp_aff->get_alpha()[lz] ;
    Matrice ope(2*nr, 2*nr) ;
    ope.set_etat_qcq() ;
	
    for (int k=0 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    // quantic numbers and spectral bases
	    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 0)) {
		Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		Diff_sx os(base_r, nr) ; const Matrice& ms = os.get_matrice() ;

		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col) = md(lin,col) + 2*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col+nr) = -l_q*(l_q+1)*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col) = -ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col+nr) = md(lin,col) + ms(lin,col) ;

		ope *= 1./alpha ;
		int ind1 = nr ;
		
		if (l_q==1) {
		ind1 += 1 ;
		int pari = 1 ;
		for (int col=0 ; col<nr; col++) {
		ope.set(nr-1,col) = pari ;
		ope.set(nr-1,col+nr) = -pari ;
		pari = - pari ;
		}
		for (int col=0 ; col<nr ; col++) {
		ope.set(2*nr-1,col+nr)=1 ;
			}
		}
		
		else{
		for (int col=0; col<2*nr; col++) {
		    ope.set(ind1+nr-2, col) = 0 ;
		}
		for (int col=nr; col<2*nr; col++) 
		    ope.set(ind1+nr-2, col) = 1 ;
		for (int col=0; col<2*nr; col++) {
		    ope.set(nr-1, col) = 0 ;
		    ope.set(2*nr-1, col) = 0 ;
		}
		int pari = 1 ;
		if (base_r == R_CHEBP) {
		    for (int col=0; col<nr; col++) {
			ope.set(nr-1, col) = pari ;
			ope.set(2*nr-1, col+nr) = pari ;
			pari = - pari ;
		    }
		}
		else { //In the odd case, the last coefficient must be zero!
		    ope.set(nr-1, nr-1) = 1 ;
		    ope.set(2*nr-1, 2*nr-1) = 1 ;
		}
				
		}

		ope.set_lu() ;
		
		Tbl sec(2*nr) ;
		sec.set_etat_qcq() ;
		for (int lin=0; lin<nr; lin++)
		    sec.set(lin) = 0 ;
		for (int lin=0; lin<nr; lin++) 
		    sec.set(nr+lin) = (*source.get_spectral_va().c_cf)
			(lz, k, j, lin) ;
		sec.set(2*nr-1) = 0 ;

		

/*		// BC is here
		if ((l_q==2)&&(k==3))  {
		sec.set(ind1+nr-2) = -5./2. ; }
		else { sec.set(ind1+nr-2) = 0 ; }*/


 
		Tbl sol = ope.inverse(sec) ;
		 		
		for (int i=0; i<nr; i++) {
		    sol_vr.set(lz, k, j, i) = sol(i) ;
		    sol_eta.set(lz, k, j, i) = sol(i+nr) ;
		}
		if ((l_q==2)&&(k==3))  {
		cout << " ========================== " << endl ;
		cout << " Operateur                  " << endl ;
		cout << " ========================== " << endl ;
		cout << ope << endl ;
		cout << " ========================== " << endl ;
		cout << " Second membre              " << endl ;
		cout << " ========================== " << endl ;
		cout << sec << endl ;
		cout << " ========================== " << endl ;
		cout << " Solution                   " << endl ;
		cout << " ========================== " << endl ;
		cout << sol << endl ;

		}
	    }
	}
    }
  }
  
    Mtbl_cf& mvr = *tilde_vr.set_spectral_va().c_cf ;
    Mtbl_cf& meta = *tilde_eta.set_spectral_va().c_cf ;
	
    Mtbl_cf d_vr = sol_vr ; 
    Mtbl_cf d_eta = sol_eta ; 
    
	
    // Loop on l and m
    //----------------
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	    if ((nullite_plm(j, nt, k, np, base) == 1) && (l_q > 0)) {
		// everything is put to the right place...
		//----------------------------------------
 		int nr = mgrid.get_nr(0) ; //nucleus
 		for (int i=0 ; i<nr ; i++) {
		    mvr.set(0, k, j, i) = sol_vr(0, k, j, i) ;
		    meta.set(0, k, j, i) = sol_eta(0, k, j, i) ;
		}
	    } // End of nullite_plm  
	} //End of loop on theta
	
	    
    if (tilde_vr.set_spectral_va().c != 0x0) 
	delete tilde_vr.set_spectral_va().c ;
    tilde_vr.set_spectral_va().c = 0x0 ;
    tilde_vr.set_spectral_va().ylm_i() ;

    if (tilde_eta.set_spectral_va().c != 0x0) 
	delete tilde_eta.set_spectral_va().c ;
    tilde_eta.set_spectral_va().c = 0x0 ;
    tilde_eta.set_spectral_va().ylm_i() ;

} 
}
