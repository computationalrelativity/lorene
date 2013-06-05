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

char sol_Dirac_A_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2013/06/05 15:10:43  j_novak
 * Suppression of FINJAC sampling in r. This Jacobi(0,2) base is now
 * available by setting colloc_r to BASE_JAC02 in the Mg3d constructor.
 *
 * Revision 1.3  2009/10/23 13:18:46  j_novak
 * Minor modifications
 *
 * Revision 1.2  2008/08/27 10:55:15  jl_cornou
 * Added the case of one zone, which is a nucleus for BC
 *
 * Revision 1.1  2008/08/27 09:01:27  jl_cornou
 * Methods for solving Dirac systems for divergence free vectors
 *
 *
 *
 * $Header$
 *
 */


// C headers
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Lorene headers
#include "metric.h"
#include "diff.h"
#include "proto.h"
#include "param.h"

//----------------------------------------------------------------------------------
//
//                               sol_Dirac_A
//
//----------------------------------------------------------------------------------

void Vector_divfree::sol_Dirac_A(const Scalar& aaa, Scalar& tilde_vr, Scalar& tilde_eta,
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
    bool cedbc = (ced && (nz_bc == nzm1)) ; 
#ifndef NDEBUG
    if (!cedbc) {
	assert(par_bc != 0x0) ;
	assert(par_bc->get_n_tbl_mod() >= 3) ;
    }
#endif
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
 
    Mtbl_cf sol_part_vr(mgrid, base) ; sol_part_vr.annule_hard() ;
    Mtbl_cf sol_part_eta(mgrid, base) ; sol_part_eta.annule_hard() ;
    Mtbl_cf sol_hom1_vr(mgrid, base) ; sol_hom1_vr.annule_hard() ;
    Mtbl_cf sol_hom1_eta(mgrid, base) ; sol_hom1_eta.annule_hard() ;
    Mtbl_cf sol_hom2_vr(mgrid, base) ; sol_hom2_vr.annule_hard() ;
    Mtbl_cf sol_hom2_eta(mgrid, base) ; sol_hom2_eta.annule_hard() ;

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
		ope.set(2*nr-1, nr) = 1;
		}
		else{
		for (int col=0; col<2*nr; col++) 
		    ope.set(ind1+nr-2, col) = 0 ;
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
		
		ope.set(ind1+nr-2, ind1) = 1 ; 
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
		sec.set(ind1+nr-2) = 0 ;
		Tbl sol = ope.inverse(sec) ;
		
		

		for (int i=0; i<nr; i++) {
		    sol_part_vr.set(lz, k, j, i) = sol(i) ;
		    sol_part_eta.set(lz, k, j, i) = sol(i+nr) ;
		}
		sec.annule_hard() ;
		sec.set(ind1+nr-2) = 1 ;
		
		sol = ope.inverse(sec) ;
		
		for (int i=0; i<nr; i++) {
		    sol_hom2_vr.set(lz, k, j, i) = sol(i) ;
		    sol_hom2_eta.set(lz, k, j, i) = sol(i+nr) ;
		}
	    }
	}
    }
    }

    //-------------
    // -- Shells --
    //-------------

    for (int lz=1; lz <= n_shell; lz++) {
	int nr = mgrid.get_nr(lz) ;
	assert(mgrid.get_nt(lz) == nt) ;
	assert(mgrid.get_np(lz) == np) ;
	double alpha = mp_aff->get_alpha()[lz] ;
	double ech = mp_aff->get_beta()[lz] / alpha ;
	Matrice ope(2*nr, 2*nr) ;
	ope.set_etat_qcq() ;
	
	for (int k=0 ; k<np+1 ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		// quantic numbers and spectral bases
		base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 0)) {
		    Diff_xdsdx oxd(base_r, nr) ; const Matrice& mxd = oxd.get_matrice() ;
		    Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		    Diff_id oid(base_r, nr) ; const Matrice& mid = oid.get_matrice() ;

		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col) = mxd(lin,col) + ech*md(lin,col) 
				+ 2*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col+nr) = -l_q*(l_q+1)*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col) = -mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col+nr) = mxd(lin,col) + ech*md(lin,col) + mid(lin,col) ;

		    int ind0 = 0 ;
		    int ind1 = nr ;
		    for (int col=0; col<2*nr; col++) {
			ope.set(ind0+nr-1, col) = 0 ;
			ope.set(ind1+nr-1, col) = 0 ;
		    }
		    ope.set(ind0+nr-1, ind0) = 1 ;
		    ope.set(ind1+nr-1, ind1) = 1 ;

		    ope.set_lu() ;

		    Tbl sec(2*nr) ;
		    sec.set_etat_qcq() ;
		    for (int lin=0; lin<nr; lin++)
			sec.set(lin) = 0 ;
		    for (int lin=0; lin<nr; lin++)
			sec.set(nr+lin) = (*source_coq.get_spectral_va().c_cf)
			    (lz, k, j, lin) ;
		    sec.set(ind0+nr-1) = 0 ;
		    sec.set(ind1+nr-1) = 0 ;
		    Tbl sol = ope.inverse(sec) ;
		    for (int i=0; i<nr; i++) {
 			sol_part_vr.set(lz, k, j, i) = sol(i) ;
 			sol_part_eta.set(lz, k, j, i) = sol(i+nr) ;
		    }
		
		
		    sec.annule_hard() ;
		    sec.set(ind0+nr-1) = 1 ;
		    sol = ope.inverse(sec) ;

		
		    for (int i=0; i<nr; i++) {
			sol_hom1_vr.set(lz, k, j, i) = sol(i) ;
			sol_hom1_eta.set(lz, k, j, i) = sol(i+nr) ;
		    }			
		    sec.set(ind0+nr-1) = 0 ;
		    sec.set(ind1+nr-1) = 1 ;
		    sol = ope.inverse(sec) ;

		
		
		    for (int i=0; i<nr; i++) {
			sol_hom2_vr.set(lz, k, j, i) = sol(i) ;
			sol_hom2_eta.set(lz, k, j, i) = sol(i+nr) ;
		    }			
		}
	    }
	}
    }

    //------------------------------
    // Compactified external domain
    //------------------------------
    if (cedbc) {int lz = nzm1 ;  
    int nr = mgrid.get_nr(lz) ;
    assert(mgrid.get_nt(lz) == nt) ;
    assert(mgrid.get_np(lz) == np) ;
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
			ope.set(lin,col) = - md(lin,col) + 2*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col+nr) = -l_q*(l_q+1)*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col) = -ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col+nr) = -md(lin,col) + ms(lin,col) ;

		ope *= 1./alpha ;
		int ind0 = 0 ;
		int ind1 = nr ;
		for (int col=0; col<2*nr; col++) {
		    ope.set(ind0+nr-1, col) = 0 ;
		    ope.set(ind1+nr-2, col) = 0 ;
		    ope.set(ind1+nr-1, col) = 0 ;
		}
		for (int col=0; col<nr; col++) {
		    ope.set(ind0+nr-1, col+ind0) = 1 ;
		    ope.set(ind1+nr-1, col+ind1) = 1 ;
		}
		ope.set(ind1+nr-2, ind1+1) = 1 ;

		ope.set_lu() ;

		Tbl sec(2*nr) ;
		sec.set_etat_qcq() ;
		for (int lin=0; lin<nr; lin++)
		    sec.set(lin) = 0 ;
		for (int lin=0; lin<nr; lin++)
		    sec.set(nr+lin) = (*source.get_spectral_va().c_cf)
			(lz, k, j, lin) ;
		sec.set(ind0+nr-1) = 0 ;
		sec.set(ind1+nr-2) = 0 ;
		sec.set(ind1+nr-1) = 0 ;
 		Tbl sol = ope.inverse(sec) ;

		for (int i=0; i<nr; i++) {
		    sol_part_vr.set(lz, k, j, i) = sol(i) ;
		    sol_part_eta.set(lz, k, j, i) = sol(i+nr) ;
		}
		sec.annule_hard() ;
		sec.set(ind1+nr-2) = 1 ;
		sol = ope.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_hom1_vr.set(lz, k, j, i) = sol(i) ;
		    sol_hom1_eta.set(lz, k, j, i) = sol(i+nr) ;
		}
	    }
	}
    }
    }

    int taille = 2*nz_bc + 1;
    if (cedbc) taille-- ;
    Mtbl_cf& mvr = *tilde_vr.set_spectral_va().c_cf ;
    Mtbl_cf& meta = *tilde_eta.set_spectral_va().c_cf ;
	
    Tbl sec_membre(taille) ; 
    Matrice systeme(taille, taille) ; 
    int ligne ;  int colonne ;
    Tbl pipo(1) ;
    const Tbl& mub = (cedbc ? pipo : par_bc->get_tbl_mod(2) );
    double c_vr = (cedbc ? 0 : par_bc->get_tbl_mod(0)(0) ) ;
    double d_vr = (cedbc ? 0 : par_bc->get_tbl_mod(0)(1) ) ;
    double c_eta = (cedbc ? 0 : par_bc->get_tbl_mod(0)(2) ) ;
    double d_eta = (cedbc ? 0 : par_bc->get_tbl_mod(0)(3) ) ;

	

    Mtbl_cf dhom1_vr = sol_hom1_vr ; 
    Mtbl_cf dhom2_vr = sol_hom2_vr ; 
    Mtbl_cf dpart_vr = sol_part_vr ; 
    Mtbl_cf dhom1_eta = sol_hom1_eta ; 
    Mtbl_cf dhom2_eta = sol_hom2_eta ; 
    Mtbl_cf dpart_eta = sol_part_eta ; 
    if (!cedbc) {
	dhom1_vr.dsdx() ;
	dhom2_vr.dsdx() ;
	dpart_vr.dsdx() ;
	dhom1_eta.dsdx() ;
	dhom2_eta.dsdx() ;
	dpart_eta.dsdx() ;
    }
	
    // Loop on l and m
    //----------------
    int nr ;
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	    if ((nullite_plm(j, nt, k, np, base) == 1) && (l_q > 0)) {
		ligne = 0 ;
		colonne = 0 ;
		systeme.annule_hard() ;
		sec_membre.annule_hard() ;


		if ((nz==1)&&(mgrid.get_type_r(0) == RARE)) {
		// Only one zone, which is a nucleus
		double alpha = mp_aff->get_alpha()[nz_bc] ;
		systeme.set(ligne, colonne) = 
		    c_vr*sol_hom2_vr.val_out_bound_jk(nz_bc, j, k) 
		    + d_vr*dhom2_vr.val_out_bound_jk(nz_bc, j, k) / alpha 
		    + c_eta*sol_hom2_eta.val_out_bound_jk(nz_bc, j, k) 
		    + d_eta*dhom2_eta.val_out_bound_jk(nz_bc, j, k) / alpha ;
		
		sec_membre.set(ligne) -= c_vr*sol_part_vr.val_out_bound_jk(nz_bc, j, k) 
		    + d_vr*dpart_vr.val_out_bound_jk(nz_bc, j, k)/alpha
		    + c_eta*sol_part_eta.val_out_bound_jk(nz_bc, j, k) 
		    + d_eta*dpart_eta.val_out_bound_jk(nz_bc, j, k)/alpha
		    - mub(k, j) ;
		}
		else {
		// General case, two zones at least
		//Nucleus 
		systeme.set(ligne, colonne) = sol_hom2_vr.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_vr.val_out_bound_jk(0, j, k) ;
		ligne++ ;

		systeme.set(ligne, colonne) = sol_hom2_eta.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_eta.val_out_bound_jk(0, j, k) ;
		colonne++ ;

		//shells
		for (int zone=1 ; zone<nz_bc ; zone++) {
		    nr = mgrid.get_nr(zone) ;
		    ligne-- ;

		    //Condition at x = -1
		    systeme.set(ligne, colonne) = 
			- sol_hom1_vr.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			- sol_hom2_vr.val_in_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) += sol_part_vr.val_in_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			- sol_hom1_eta.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			- sol_hom2_eta.val_in_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) += sol_part_eta.val_in_bound_jk(zone, j, k) ;
		    ligne++ ;

		    // Condition at x=1
		    systeme.set(ligne, colonne) = 
			sol_hom1_vr.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			sol_hom2_vr.val_out_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) -= sol_part_vr.val_out_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			sol_hom1_eta.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			sol_hom2_eta.val_out_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) -= sol_part_eta.val_out_bound_jk(zone, j, k) ;
		    
		    colonne += 2 ;
		}
    
		//Last  domain	 
		nr = mgrid.get_nr(nz_bc) ;
		double alpha = mp_aff->get_alpha()[nz_bc] ;
		ligne-- ;
		    
		//Condition at x = -1
		systeme.set(ligne, colonne) = 
		    - sol_hom1_vr.val_in_bound_jk(nz_bc, j, k) ;
		if (!cedbc) systeme.set(ligne, colonne+1) = 
				- sol_hom2_vr.val_in_bound_jk(nz_bc, j, k) ;
		
		sec_membre.set(ligne) += sol_part_vr.val_in_bound_jk(nz_bc, j, k) ;
		ligne++ ;
		    
		systeme.set(ligne, colonne) = 
		    - sol_hom1_eta.val_in_bound_jk(nz_bc, j, k) ;
		if (!cedbc) systeme.set(ligne, colonne+1) = 
				- sol_hom2_eta.val_in_bound_jk(nz_bc, j, k) ;
		    
		sec_membre.set(ligne) += sol_part_eta.val_in_bound_jk(nz_bc, j, k) ;
		ligne++ ;

		if (!cedbc) {// Special condition at x=1
		systeme.set(ligne, colonne) = 
		    c_vr*sol_hom1_vr.val_out_bound_jk(nz_bc, j, k) 
		    + d_vr*dhom1_vr.val_out_bound_jk(nz_bc, j, k) / alpha 
		    + c_eta*sol_hom1_eta.val_out_bound_jk(nz_bc, j, k) 
		    + d_eta*dhom1_eta.val_out_bound_jk(nz_bc, j, k) / alpha ;
		systeme.set(ligne, colonne+1) = 
		    c_vr*sol_hom2_vr.val_out_bound_jk(nz_bc, j, k) 
		    + d_vr*dhom2_vr.val_out_bound_jk(nz_bc, j, k) / alpha
		    + c_eta*sol_hom2_eta.val_out_bound_jk(nz_bc, j, k) 
		    + d_eta*dhom2_eta.val_out_bound_jk(nz_bc, j, k) / alpha ;

		sec_membre.set(ligne) -= c_vr*sol_part_vr.val_out_bound_jk(nz_bc, j, k) 
		    + d_vr*dpart_vr.val_out_bound_jk(nz_bc, j, k)/alpha
		    + c_eta*sol_part_eta.val_out_bound_jk(nz_bc, j, k) 
		    + d_eta*dpart_eta.val_out_bound_jk(nz_bc, j, k)/alpha
		    - mub(k, j) ;
		}
		}

		// Solution of the system giving the coefficients for the homogeneous 
		// solutions
		//-------------------------------------------------------------------
		systeme.set_lu() ;
		Tbl facteur = systeme.inverse(sec_membre) ;
		int conte = 0 ;

		
		// everything is put to the right place...
		//----------------------------------------
 		nr = mgrid.get_nr(0) ; //nucleus
 		for (int i=0 ; i<nr ; i++) {
		    mvr.set(0, k, j, i) = sol_part_vr(0, k, j, i)
			+ facteur(conte)*sol_hom2_vr(0, k, j, i) ;
		    meta.set(0, k, j, i) = sol_part_eta(0, k, j, i)
			+ facteur(conte)*sol_hom2_eta(0, k, j, i) ;
 		}
 		conte++ ;
 		for (int zone=1 ; zone<=n_shell ; zone++) { //shells
 		    nr = mgrid.get_nr(zone) ;
 		    for (int i=0 ; i<nr ; i++) {
		    mvr.set(zone, k, j, i) = sol_part_vr(zone, k, j, i)
			+ facteur(conte)*sol_hom1_vr(zone, k, j, i) 
			+ facteur(conte+1)*sol_hom2_vr(zone, k, j, i) ;
			
		    meta.set(zone, k, j, i) = sol_part_eta(zone, k, j, i)
			+ facteur(conte)*sol_hom1_eta(zone, k, j, i) 
			+ facteur(conte+1)*sol_hom2_eta(zone, k, j, i) ;
 		    }
 		    conte+=2 ;
 		}
		if (cedbc) {
		    nr = mgrid.get_nr(nzm1) ; //compactified external domain
		    for (int i=0 ; i<nr ; i++) {
			mvr.set(nzm1, k, j, i) = sol_part_vr(nzm1, k, j, i)
			    + facteur(conte)*sol_hom1_vr(nzm1, k, j, i) ;
			
			meta.set(nzm1, k, j, i) = sol_part_eta(nzm1, k, j, i)
			    + facteur(conte)*sol_hom1_eta(nzm1, k, j, i) ;
		    }
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
