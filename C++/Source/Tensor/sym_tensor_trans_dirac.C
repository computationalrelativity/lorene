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

char sym_tensor_trans_dirac_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2015/08/10 15:32:26  j_novak
 * Better calls to Param::add_int(), to avoid weird problems (e.g. with g++ 4.8).
 *
 * Revision 1.7  2014/10/13 08:53:43  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2008/08/29 13:15:22  j_novak
 * Corrected a mistake in the case of no CED.
 *
 * Revision 1.4  2008/08/29 05:33:18  j_novak
 * Minor modif.
 *
 * Revision 1.3  2008/08/27 08:13:20  j_novak
 * Correction of a mistake in the imposition of BCs in sol_Dirac_A. + Treatment of
 * the case of BCs imposed on a nucleus (nz_bc = 0).
 *
 * Revision 1.2  2006/10/24 13:03:19  j_novak
 * New methods for the solution of the tensor wave equation. Perhaps, first
 * operational version...
 *
 * Revision 1.1  2006/09/05 15:38:45  j_novak
 * The fuctions sol_Dirac... are in a seperate file, with new parameters to
 * control the boundary conditions.
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
#include "diff.h"
#include "proto.h"
#include "param.h"

//----------------------------------------------------------------------------------
//
//                               sol_Dirac_A
//
//----------------------------------------------------------------------------------

namespace Lorene {
void Sym_tensor_trans::sol_Dirac_A(const Scalar& aaa, Scalar& tilde_mu, Scalar& x_new,
				   const Param* par_bc) const {

    const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
    assert(mp_aff != 0x0) ; //Only affine mapping for the moment

    const Mg3d& mgrid = *mp_aff->get_mg() ;
    assert(mgrid.get_type_r(0) == RARE)  ;
    if (aaa.get_etat() == ETATZERO) {
	tilde_mu = 0 ;
	x_new = 0 ;
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

    tilde_mu.annule_hard() ;
    tilde_mu.set_spectral_base(base) ;
    tilde_mu.set_spectral_va().set_etat_cf_qcq() ;
    tilde_mu.set_spectral_va().c_cf->annule_hard() ;   
    x_new.annule_hard() ;
    x_new.set_spectral_base(base) ;
    x_new.set_spectral_va().set_etat_cf_qcq() ;
    x_new.set_spectral_va().c_cf->annule_hard() ;   
 
    Mtbl_cf sol_part_mu(mgrid, base) ; sol_part_mu.annule_hard() ;
    Mtbl_cf sol_part_x(mgrid, base) ; sol_part_x.annule_hard() ;
    Mtbl_cf sol_hom1_mu(mgrid, base) ; sol_hom1_mu.annule_hard() ;
    Mtbl_cf sol_hom1_x(mgrid, base) ; sol_hom1_x.annule_hard() ;
    Mtbl_cf sol_hom2_mu(mgrid, base) ; sol_hom2_mu.annule_hard() ;
    Mtbl_cf sol_hom2_x(mgrid, base) ; sol_hom2_x.annule_hard() ;

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
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)) {
		Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		Diff_sx os(base_r, nr) ; const Matrice& ms = os.get_matrice() ;

		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col) = md(lin,col) + 3*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col+nr) = (2-l_q*(l_q+1))*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col) = -ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col+nr) = md(lin,col) ;

		ope *= 1./alpha ;
		int ind1 = nr ;
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
		    sol_part_mu.set(lz, k, j, i) = sol(i) ;
		    sol_part_x.set(lz, k, j, i) = sol(i+nr) ;
		}
		sec.annule_hard() ;
		sec.set(ind1+nr-2) = 1 ;
		sol = ope.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_hom2_mu.set(lz, k, j, i) = sol(i) ;
		    sol_hom2_x.set(lz, k, j, i) = sol(i+nr) ;
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
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)) {
		    Diff_xdsdx oxd(base_r, nr) ; const Matrice& mxd = oxd.get_matrice() ;
		    Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		    Diff_id oid(base_r, nr) ; const Matrice& mid = oid.get_matrice() ;

		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col) = mxd(lin,col) + ech*md(lin,col) 
				+ 3*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col+nr) = (2-l_q*(l_q+1))*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col) = -mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col+nr) = mxd(lin,col) + ech*md(lin,col) ;

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
 			sol_part_mu.set(lz, k, j, i) = sol(i) ;
 			sol_part_x.set(lz, k, j, i) = sol(i+nr) ;
		    }
		    sec.annule_hard() ;
		    sec.set(ind0+nr-1) = 1 ;
		    sol = ope.inverse(sec) ;
		    for (int i=0; i<nr; i++) {
			sol_hom1_mu.set(lz, k, j, i) = sol(i) ;
			sol_hom1_x.set(lz, k, j, i) = sol(i+nr) ;
		    }			
		    sec.set(ind0+nr-1) = 0 ;
		    sec.set(ind1+nr-1) = 1 ;
		    sol = ope.inverse(sec) ;
		    for (int i=0; i<nr; i++) {
			sol_hom2_mu.set(lz, k, j, i) = sol(i) ;
			sol_hom2_x.set(lz, k, j, i) = sol(i+nr) ;
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
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)) {
		Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		Diff_sx os(base_r, nr) ; const Matrice& ms = os.get_matrice() ;

		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col) = - md(lin,col) + 3*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col+nr) = (2-l_q*(l_q+1))*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col) = -ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col+nr) = -md(lin,col) ;

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
		    sol_part_mu.set(lz, k, j, i) = sol(i) ;
		    sol_part_x.set(lz, k, j, i) = sol(i+nr) ;
		}
		sec.annule_hard() ;
		sec.set(ind1+nr-2) = 1 ;
		sol = ope.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_hom1_mu.set(lz, k, j, i) = sol(i) ;
		    sol_hom1_x.set(lz, k, j, i) = sol(i+nr) ;
		}			
	    }
	}
    }
    }

    //---------------------------------------------------------------
    // Matching of the solutions across the domains and imposition of 
    // boundary conditions (if no compactified domain)
    //---------------------------------------------------------------
    int taille = 2*nz_bc + 1;
    if (cedbc) taille-- ;
    Mtbl_cf& mmu = *tilde_mu.set_spectral_va().c_cf ;
    Mtbl_cf& mw = *x_new.set_spectral_va().c_cf ;
	
    Tbl sec_membre(taille) ; 
    Matrice systeme(taille, taille) ; 
    int ligne ;  int colonne ;
    Tbl pipo(1) ;
    const Tbl& mub = (cedbc ? pipo : par_bc->get_tbl_mod(2) );
    double c_mu = (cedbc ? 0 : par_bc->get_tbl_mod(0)(0) ) ;
    double d_mu = (cedbc ? 0 : par_bc->get_tbl_mod(0)(1) ) ;
    double c_x = (cedbc ? 0 : par_bc->get_tbl_mod(0)(2) ) ;
    double d_x = (cedbc ? 0 : par_bc->get_tbl_mod(0)(3) ) ;
    Mtbl_cf dhom1_mu = sol_hom1_mu ; 
    Mtbl_cf dhom2_mu = sol_hom2_mu ; 
    Mtbl_cf dpart_mu = sol_part_mu ; 
    Mtbl_cf dhom1_x = sol_hom1_x ; 
    Mtbl_cf dhom2_x = sol_hom2_x ; 
    Mtbl_cf dpart_x = sol_part_x ; 
    if (!cedbc) {
	dhom1_mu.dsdx() ;
	dhom2_mu.dsdx() ;
	dpart_mu.dsdx() ;
	dhom1_x.dsdx() ;
	dhom2_x.dsdx() ;
	dpart_x.dsdx() ;
    }
	
    // Loop on l and m
    //----------------
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	    if ((nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)) {
		ligne = 0 ;
		colonne = 0 ;
		systeme.annule_hard() ;
		sec_membre.annule_hard() ;

		//Nucleus 
		int nr = mgrid.get_nr(0) ;
		
		if (nz_bc >0) {
		    systeme.set(ligne, colonne) = sol_hom2_mu.val_out_bound_jk(0, j, k) ;
		    sec_membre.set(ligne) = -sol_part_mu.val_out_bound_jk(0, j, k) ;
		    ligne++ ;
		    
		    systeme.set(ligne, colonne) = sol_hom2_x.val_out_bound_jk(0, j, k) ;
		    sec_membre.set(ligne) = -sol_part_x.val_out_bound_jk(0, j, k) ;
		    colonne++ ;
		}
		//shells
		for (int zone=1 ; zone<nz_bc ; zone++) {
		    nr = mgrid.get_nr(zone) ;
		    ligne-- ;

		    //Condition at x = -1
		    systeme.set(ligne, colonne) = 
			- sol_hom1_mu.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			- sol_hom2_mu.val_in_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) += sol_part_mu.val_in_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			- sol_hom1_x.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			- sol_hom2_x.val_in_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) += sol_part_x.val_in_bound_jk(zone, j, k) ;
		    ligne++ ;

		    // Condition at x=1
		    systeme.set(ligne, colonne) = 
			sol_hom1_mu.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			sol_hom2_mu.val_out_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) -= sol_part_mu.val_out_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			sol_hom1_x.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			sol_hom2_x.val_out_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) -= sol_part_x.val_out_bound_jk(zone, j, k) ;
		    
		    colonne += 2 ;
		}
    
		//Last  domain	 
		nr = mgrid.get_nr(nz_bc) ;
		double alpha = mp_aff->get_alpha()[nz_bc] ;
		if (nz_bc>0) {
		    ligne-- ;
		    
		    //Condition at x = -1
		    systeme.set(ligne, colonne) = 
			- sol_hom1_mu.val_in_bound_jk(nz_bc, j, k) ;
		    if (!cedbc) systeme.set(ligne, colonne+1) = 
				    - sol_hom2_mu.val_in_bound_jk(nz_bc, j, k) ;
		    
		    sec_membre.set(ligne) += sol_part_mu.val_in_bound_jk(nz_bc, j, k) ;
		    ligne++ ;
		    
		    systeme.set(ligne, colonne) = 
			- sol_hom1_x.val_in_bound_jk(nz_bc, j, k) ;
		    if (!cedbc) systeme.set(ligne, colonne+1) = 
				    - sol_hom2_x.val_in_bound_jk(nz_bc, j, k) ;
		    
		    sec_membre.set(ligne) += sol_part_x.val_in_bound_jk(nz_bc, j, k) ;
		    ligne++ ;
		}
		if (!cedbc) {// Special condition at x=1
		    if (nz_bc>0) {
		systeme.set(ligne, colonne) = 
		    c_mu*sol_hom1_mu.val_out_bound_jk(nz_bc, j, k) 
		    + d_mu*dhom1_mu.val_out_bound_jk(nz_bc, j, k) / alpha 
		    + c_x*sol_hom1_x.val_out_bound_jk(nz_bc, j, k) 
		    + d_x*dhom1_x.val_out_bound_jk(nz_bc, j, k) / alpha ;
		    }
		    else {
			assert(ligne == 0) ;
			colonne = -1 ;
		    }
		systeme.set(ligne, colonne+1) = 
		    c_mu*sol_hom2_mu.val_out_bound_jk(nz_bc, j, k) 
		    + d_mu*dhom2_mu.val_out_bound_jk(nz_bc, j, k) / alpha
		    + c_x*sol_hom2_x.val_out_bound_jk(nz_bc, j, k) 
		    + d_x*dhom2_x.val_out_bound_jk(nz_bc, j, k) / alpha ;

		sec_membre.set(ligne) -= c_mu*sol_part_mu.val_out_bound_jk(nz_bc, j, k) 
		    + d_mu*dpart_mu.val_out_bound_jk(nz_bc, j, k)/alpha
		    + c_x*sol_part_x.val_out_bound_jk(nz_bc, j, k) 
		    + d_x*dpart_x.val_out_bound_jk(nz_bc, j, k)/alpha
		    - mub(k, j) ;
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
		    mmu.set(0, k, j, i) = sol_part_mu(0, k, j, i)
			+ facteur(conte)*sol_hom2_mu(0, k, j, i) ;
		    mw.set(0, k, j, i) = sol_part_x(0, k, j, i)
			+ facteur(conte)*sol_hom2_x(0, k, j, i) ;
 		}
 		conte++ ;
 		for (int zone=1 ; zone<=n_shell ; zone++) { //shells
 		    nr = mgrid.get_nr(zone) ;
 		    for (int i=0 ; i<nr ; i++) {
		    mmu.set(zone, k, j, i) = sol_part_mu(zone, k, j, i)
			+ facteur(conte)*sol_hom1_mu(zone, k, j, i) 
			+ facteur(conte+1)*sol_hom2_mu(zone, k, j, i) ;
			
		    mw.set(zone, k, j, i) = sol_part_x(zone, k, j, i)
			+ facteur(conte)*sol_hom1_x(zone, k, j, i) 
			+ facteur(conte+1)*sol_hom2_x(zone, k, j, i) ;
 		    }
 		    conte+=2 ;
 		}
		if (cedbc) {
		    nr = mgrid.get_nr(nzm1) ; //compactified external domain
		    for (int i=0 ; i<nr ; i++) {
			mmu.set(nzm1, k, j, i) = sol_part_mu(nzm1, k, j, i)
			    + facteur(conte)*sol_hom1_mu(nzm1, k, j, i) ;
			
			mw.set(nzm1, k, j, i) = sol_part_x(nzm1, k, j, i)
			    + facteur(conte)*sol_hom1_x(nzm1, k, j, i) ;
		    }
		}

	    } // End of nullite_plm  
	} //End of loop on theta
		    
    if (tilde_mu.set_spectral_va().c != 0x0) 
	delete tilde_mu.set_spectral_va().c ;
    tilde_mu.set_spectral_va().c = 0x0 ;
    tilde_mu.set_spectral_va().ylm_i() ;

    if (x_new.set_spectral_va().c != 0x0) 
	delete x_new.set_spectral_va().c ;
    x_new.set_spectral_va().c = 0x0 ;
    x_new.set_spectral_va().ylm_i() ;

} 

//----------------------------------------------------------------------------------
//
//                               sol_Dirac_tilde_B
//
//----------------------------------------------------------------------------------

void Sym_tensor_trans::sol_Dirac_tilde_B(const Scalar& tilde_b, const Scalar& hh, 
					 Scalar& hrr, Scalar& tilde_eta, Scalar& ww,
					 Param* par_bc, Param* par_mat) const {

    const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
    assert(mp_aff != 0x0) ; //Only affine mapping for the moment

    const Mg3d& mgrid = *mp_aff->get_mg() ;
    assert(mgrid.get_type_r(0) == RARE)  ;
    if ( (tilde_b.get_etat() == ETATZERO) && (hh.get_etat() == ETATZERO) ) {
	hrr = 0 ;
	tilde_eta = 0 ;
	ww = 0 ;
	return ;
    }
    int nz = mgrid.get_nzone() ;
    int nzm1 = nz - 1 ;
    bool ced = (mgrid.get_type_r(nzm1) == UNSURR) ;
    int n_shell = ced ? nz-2 : nzm1 ;
    int nz_bc = nzm1 ;
    if (par_bc != 0x0)
	if (par_bc->get_n_int() > 0)
	    nz_bc = par_bc->get_int() ;
    n_shell = (nz_bc < n_shell ? nz_bc : n_shell) ;
    bool cedbc = (ced && (nz_bc == nzm1)) ; 
#ifndef NDEBUG
    if (!cedbc) {
	assert(par_bc != 0x0) ;
	assert(par_bc->get_n_tbl_mod() >= 2) ;
    }
#endif
    int nt = mgrid.get_nt(0) ;
    int np = mgrid.get_np(0) ;

    assert (tilde_b.get_etat() != ETATNONDEF) ;
    assert (hh.get_etat() != ETATNONDEF) ;

    Scalar source = tilde_b ;
    Scalar source_coq = tilde_b ;
    source_coq.annule_domain(0) ;
    if (ced)
	source_coq.annule_domain(nzm1) ;
    source_coq.mult_r() ;
    source.set_spectral_va().ylm() ;
    source_coq.set_spectral_va().ylm() ;
    bool bnull = (tilde_b.get_etat() == ETATZERO) ;

    assert(hh.check_dzpuis(0)) ;
    Scalar hoverr = hh ;
    hoverr.div_r_dzpuis(2) ;
    hoverr.set_spectral_va().ylm() ;
    Scalar dhdr = hh.dsdr() ;
    dhdr.set_spectral_va().ylm() ;
    Scalar h_coq = hh ;
    h_coq.set_spectral_va().ylm() ;
    Scalar dh_coq = hh.dsdr() ;
    dh_coq.mult_r_dzpuis(0) ;
    dh_coq.set_spectral_va().ylm() ;    
    bool hnull = (hh.get_etat() == ETATZERO) ;

    Base_val base = (bnull ? hoverr.get_spectral_base() : source.get_spectral_base()) ;
    base.mult_x() ;
    int lmax = base.give_lmax(mgrid, 0) + 1;

    bool need_calculation = true ;
    if (par_mat != 0x0) {
	bool param_new = false ;
	if ((par_mat->get_n_int_mod() >= 4)
	    &&(par_mat->get_n_tbl_mod()>=1)
	    &&(par_mat->get_n_matrice_mod()>=1)
	    &&(par_mat->get_n_itbl_mod()>=1)) {
	    if (par_mat->get_int_mod(0) < nz_bc) param_new = true ;
	    if (par_mat->get_int_mod(1) != lmax) param_new = true ;
	    if (par_mat->get_int_mod(2) != mgrid.get_type_t() ) param_new = true ;
	    if (par_mat->get_int_mod(3) != mgrid.get_type_p() ) param_new = true ;
	    if (par_mat->get_itbl_mod(0)(0) != mgrid.get_nr(0)) param_new = true ;
	    if (fabs(par_mat->get_tbl_mod(0)(0) - mp_aff->get_alpha()[0]) > 2.e-15)
		param_new = true ; 
	    for (int l=1; l<= n_shell; l++) {
		if (par_mat->get_itbl_mod(0)(l) != mgrid.get_nr(l)) param_new = true ;
		if (fabs(par_mat->get_tbl_mod(0)(l) - mp_aff->get_beta()[l] / 
		    mp_aff->get_alpha()[l]) > 2.e-15) param_new = true ;
	    }
	    if (ced) {
	    if (par_mat->get_itbl_mod(0)(nzm1) != mgrid.get_nr(nzm1)) param_new = true ;
	    if (fabs(par_mat->get_tbl_mod(0)(nzm1) - mp_aff->get_alpha()[nzm1]) > 2.e-15)
		param_new = true ; 
	    }
	}
	else{
	    param_new = true ;
	}
	if (param_new) {
	    par_mat->clean_all() ;
	    int* nz_bc_new = new int(nz_bc) ;
	    par_mat->add_int_mod(*nz_bc_new, 0) ;
	    int* lmax_new = new int(lmax) ;
	    par_mat->add_int_mod(*lmax_new, 1) ;
	    int* type_t_new = new int(mgrid.get_type_t()) ;
	    par_mat->add_int_mod(*type_t_new, 2) ;
	    int* type_p_new = new int(mgrid.get_type_p()) ;
	    par_mat->add_int_mod(*type_p_new, 3) ;
	    Itbl* pnr = new Itbl(nz) ;
	    pnr->set_etat_qcq() ;
	    par_mat->add_itbl_mod(*pnr) ;
	    for (int l=0; l<nz; l++)
		pnr->set(l) = mgrid.get_nr(l) ;
	    Tbl* palpha = new Tbl(nz) ;
	    palpha->set_etat_qcq() ;
	    par_mat->add_tbl_mod(*palpha) ;
	    palpha->set(0) = mp_aff->get_alpha()[0] ;
	    for (int l=1; l<nzm1; l++)
		palpha->set(l) = mp_aff->get_beta()[l] / mp_aff->get_alpha()[l] ;
	    palpha->set(nzm1) = mp_aff->get_alpha()[nzm1] ;
	 }
	else need_calculation = false ;
    }
	    
    hrr.set_etat_qcq() ;
    hrr.set_spectral_base(base) ;
    hrr.set_spectral_va().set_etat_cf_qcq() ;
    hrr.set_spectral_va().c_cf->annule_hard() ;   
    tilde_eta.annule_hard() ;
    tilde_eta.set_spectral_base(base) ;
    tilde_eta.set_spectral_va().set_etat_cf_qcq() ;
    tilde_eta.set_spectral_va().c_cf->annule_hard() ;   
    ww.annule_hard() ;
    ww.set_spectral_base(base) ;
    ww.set_spectral_va().set_etat_cf_qcq() ;
    ww.set_spectral_va().c_cf->annule_hard() ;   

    sol_Dirac_l01(hh, hrr, tilde_eta, par_mat) ;
    tilde_eta.annule_l(0,0, true) ;
 
    Mtbl_cf sol_part_hrr(mgrid, base) ; sol_part_hrr.annule_hard() ;
    Mtbl_cf sol_part_eta(mgrid, base) ; sol_part_eta.annule_hard() ;
    Mtbl_cf sol_part_w(mgrid, base) ; sol_part_w.annule_hard() ;
    Mtbl_cf sol_hom1_hrr(mgrid, base) ; sol_hom1_hrr.annule_hard() ;
    Mtbl_cf sol_hom1_eta(mgrid, base) ; sol_hom1_eta.annule_hard() ;
    Mtbl_cf sol_hom1_w(mgrid, base) ; sol_hom1_w.annule_hard() ;
    Mtbl_cf sol_hom2_hrr(mgrid, base) ; sol_hom2_hrr.annule_hard() ;
    Mtbl_cf sol_hom2_eta(mgrid, base) ; sol_hom2_eta.annule_hard() ;
    Mtbl_cf sol_hom2_w(mgrid, base) ; sol_hom2_w.annule_hard() ;
    Mtbl_cf sol_hom3_hrr(mgrid, base) ; sol_hom3_hrr.annule_hard() ;
    Mtbl_cf sol_hom3_eta(mgrid, base) ; sol_hom3_eta.annule_hard() ;
    Mtbl_cf sol_hom3_w(mgrid, base) ; sol_hom3_w.annule_hard() ;

    int l_q, m_q, base_r ;
    Itbl mat_done(lmax) ;

    //---------------
    //--  NUCLEUS ---
    //---------------
    {int lz = 0 ;  
    int nr = mgrid.get_nr(lz) ;
    double alpha = mp_aff->get_alpha()[lz] ;
    Matrice ope(3*nr, 3*nr) ;
    int ind2 = 2*nr ;
    if (need_calculation && (par_mat != 0x0)) mat_done.annule_hard() ;
		
    for (int k=0 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    // quantic numbers and spectral bases
	    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)) {
		if (need_calculation) {
		    ope.set_etat_qcq() ;
		    Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		    Diff_sx os(base_r, nr) ; const Matrice& ms = os.get_matrice() ;
		    
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col) = md(lin,col) + 3*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col+nr) = -l_q*(l_q+1)*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col+2*nr) = 0 ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col) = -0.5*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col+nr) = md(lin,col) + 3*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col+2*nr) = (2. - l_q*(l_q+1))*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+2*nr,col) = -0.5*md(lin,col)/double(l_q+1) 
				- 0.5*double(l_q+4)/double(l_q+1)*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+2*nr,col+nr) = -2*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+2*nr,col+2*nr) =  (l_q+2)*md(lin,col) 
				+ l_q*(l_q+2)*ms(lin,col) ;
		    
		    ope *= 1./alpha ;
		    for (int col=0; col<3*nr; col++) 
			if (l_q>2) ope.set(ind2+nr-2, col) = 0 ;
		    for (int col=0; col<3*nr; col++) {
			ope.set(nr-1, col) = 0 ;
			ope.set(2*nr-1, col) = 0 ;
			ope.set(3*nr-1, col) = 0 ;
		    }
		    int pari = 1 ;
		    if (base_r == R_CHEBP) {
			for (int col=0; col<nr; col++) {
			    ope.set(nr-1, col) = pari ;
			    ope.set(2*nr-1, col+nr) = pari ;
			    ope.set(3*nr-1, col+2*nr) = pari ;
			    pari = - pari ;
			}
		    }
		    else { //In the odd case, the last coefficient must be zero!
			ope.set(nr-1, nr-1) = 1 ;
			ope.set(2*nr-1, 2*nr-1) = 1 ;
			ope.set(3*nr-1, 3*nr-1) = 1 ;
		    }			
		    if (l_q>2) 
			ope.set(ind2+nr-2, ind2) = 1 ;
		    
		    ope.set_lu() ;
		    if ((par_mat != 0x0) && (mat_done(l_q) == 0)) {
			Matrice* pope = new Matrice(ope) ;
			par_mat->add_matrice_mod(*pope, lz*lmax + l_q) ;
			mat_done.set(l_q) = 1 ;
		    }
		} //End of case when a calculation is needed

		const Matrice& oper = (par_mat == 0x0 ? ope : 
				       par_mat->get_matrice_mod(lz*lmax + l_q) ) ;
		Tbl sec(3*nr) ;
		sec.set_etat_qcq() ;
		if (hnull) {
		    for (int lin=0; lin<2*nr; lin++)
			sec.set(lin) = 0 ;
		    for (int lin=0; lin<nr; lin++)
			sec.set(2*nr+lin) = (*source.get_spectral_va().c_cf)
			    (lz, k, j, lin) ;
		}
		else {
		    for (int lin=0; lin<nr; lin++)
			sec.set(lin) = (*hoverr.get_spectral_va().c_cf)(lz, k, j, lin) ;
		    for (int lin=0; lin<nr; lin++)
			sec.set(lin+nr) = -0.5*(*hoverr.get_spectral_va().c_cf)
			    (lz, k, j, lin) ;
		    if (bnull) {
			for (int lin=0; lin<nr; lin++)
			    sec.set(2*nr+lin) = -0.5/double(l_q+1)*(
				(*dhdr.get_spectral_va().c_cf)(lz, k, j, lin)
			    + (l_q+2)*(*hoverr.get_spectral_va().c_cf)(lz, k, j, lin) ) ;
		    }
		    else {
			for (int lin=0; lin<nr; lin++)
			    sec.set(2*nr+lin) = -0.5/double(l_q+1)*(
				(*dhdr.get_spectral_va().c_cf)(lz, k, j, lin)
			    + (l_q+2)*(*hoverr.get_spectral_va().c_cf)(lz, k, j, lin) )
				+ (*source.get_spectral_va().c_cf)(lz, k, j, lin) ;
		    }			
		}
		if (l_q>2) sec.set(ind2+nr-2) = 0 ;
		sec.set(3*nr-1) = 0 ;
		Tbl sol = oper.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_part_hrr.set(lz, k, j, i) = sol(i) ;
		    sol_part_eta.set(lz, k, j, i) = sol(i+nr) ;
		    sol_part_w.set(lz, k, j, i) = sol(i+2*nr) ;
		}
		sec.annule_hard() ;
		if (l_q>2) {
		    sec.set(ind2+nr-2) = 1 ;
		    sol = oper.inverse(sec) ;
		}
		else { //Homogeneous solution put in by hand in the case l=2
		    sol.annule_hard() ;
		    sol.set(0) = 4 ;
		    sol.set(nr) = 2 ;
		    sol.set(2*nr) = 1 ;
		}
		for (int i=0; i<nr; i++) {
		    sol_hom3_hrr.set(lz, k, j, i) = sol(i) ;
		    sol_hom3_eta.set(lz, k, j, i) = sol(i+nr) ;
		    sol_hom3_w.set(lz, k, j, i) = sol(i+2*nr) ;
		}
	    }
	}
    }
    }

    //-------------
    // -- Shells --
    //-------------

    for (int lz=1; lz<= n_shell; lz++) {
	if (need_calculation && (par_mat != 0x0)) mat_done.annule_hard() ;
	int nr = mgrid.get_nr(lz) ;
	int ind0 = 0 ;
	int ind1 = nr ;
	int ind2 = 2*nr ;
	double alpha = mp_aff->get_alpha()[lz] ;
	double ech = mp_aff->get_beta()[lz] / alpha ;
	Matrice ope(3*nr, 3*nr) ;
	
	for (int k=0 ; k<np+1 ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		// quantic numbers and spectral bases
		base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)) {
		    if (need_calculation) {
		    ope.set_etat_qcq() ;
		    Diff_xdsdx oxd(base_r, nr) ; const Matrice& mxd = oxd.get_matrice() ;
		    Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		    Diff_id oid(base_r, nr) ; const Matrice& mid = oid.get_matrice() ;

		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col) = mxd(lin,col) + ech*md(lin,col) 
				+ 3*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col+nr) = -l_q*(l_q+1)*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col+2*nr) = 0 ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col) = -0.5*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col+nr) = mxd(lin,col) + ech*md(lin,col) 
				+ 3*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col+2*nr) = (2. - l_q*(l_q+1))*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+2*nr,col) = 
				-0.5/double(l_q+1)*(mxd(lin,col) + ech*md(lin,col)
						    + double(l_q+4)*mid(lin,col)) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+2*nr,col+nr) = -2*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+2*nr,col+2*nr) =  
				double(l_q+2)*(mxd(lin,col) + ech*md(lin,col) 
					       + l_q*mid(lin,col)) ;
		    for (int col=0; col<3*nr; col++) {
			ope.set(ind0+nr-1, col) = 0 ;
			ope.set(ind1+nr-1, col) = 0 ;
			ope.set(ind2+nr-1, col) = 0 ;
		    }
		    ope.set(ind0+nr-1, ind0) = 1 ;
		    ope.set(ind1+nr-1, ind1) = 1 ;
		    ope.set(ind2+nr-1, ind2) = 1 ;

		    ope.set_lu() ;
		    if ((par_mat != 0x0) && (mat_done(l_q) == 0)) {
			Matrice* pope = new Matrice(ope) ;
			par_mat->add_matrice_mod(*pope, lz*lmax + l_q) ;
			mat_done.set(l_q) = 1 ;
		    }
		    } //End of case when a calculation is needed
		    const Matrice& oper = (par_mat == 0x0 ? ope : 
				       par_mat->get_matrice_mod(lz*lmax + l_q) ) ;
		    Tbl sec(3*nr) ;
		    sec.set_etat_qcq() ;
		    if (hnull) {
			for (int lin=0; lin<2*nr; lin++)
			    sec.set(lin) = 0 ;
			for (int lin=0; lin<nr; lin++)
			    sec.set(2*nr+lin) = (*source_coq.get_spectral_va().c_cf)
				(lz, k, j, lin) ;
		    }
		    else {
		    for (int lin=0; lin<nr; lin++)
			sec.set(lin) = (*h_coq.get_spectral_va().c_cf)(lz, k, j, lin) ;
		    for (int lin=0; lin<nr; lin++)
			sec.set(lin+nr) = -0.5*(*h_coq.get_spectral_va().c_cf)
			    (lz, k, j, lin) ;
		    if (bnull) {
		    for (int lin=0; lin<nr; lin++)
			sec.set(2*nr+lin) = -0.5/double(l_q+1)*(
			    (*dh_coq.get_spectral_va().c_cf)(lz, k, j, lin)
			    + (l_q+2)*(*h_coq.get_spectral_va().c_cf)(lz, k, j, lin) ) ;
		    }
		    else {
		    for (int lin=0; lin<nr; lin++)
			sec.set(2*nr+lin) = -0.5/double(l_q+1)*(
			    (*dh_coq.get_spectral_va().c_cf)(lz, k, j, lin)
			    + (l_q+2)*(*h_coq.get_spectral_va().c_cf)(lz, k, j, lin) )
			    + (*source_coq.get_spectral_va().c_cf)(lz, k, j, lin) ;
		    }
		    }
		    sec.set(ind0+nr-1) = 0 ;
		    sec.set(ind1+nr-1) = 0 ;
		    sec.set(ind2+nr-1) = 0 ;
		    Tbl sol = oper.inverse(sec) ;
		    for (int i=0; i<nr; i++) {
 			sol_part_hrr.set(lz, k, j, i) = sol(i) ;
 			sol_part_eta.set(lz, k, j, i) = sol(i+nr) ;
 			sol_part_w.set(lz, k, j, i) = sol(i+2*nr) ;
		    }
		    sec.annule_hard() ;
		    sec.set(ind0+nr-1) = 1 ;
		    sol = oper.inverse(sec) ;
		    for (int i=0; i<nr; i++) {
			sol_hom1_hrr.set(lz, k, j, i) = sol(i) ;
			sol_hom1_eta.set(lz, k, j, i) = sol(i+nr) ;
			sol_hom1_w.set(lz, k, j, i) = sol(i+2*nr) ;
		    }			
		    sec.set(ind0+nr-1) = 0 ;
		    sec.set(ind1+nr-1) = 1 ;
		    sol = oper.inverse(sec) ;
		    for (int i=0; i<nr; i++) {
			sol_hom2_hrr.set(lz, k, j, i) = sol(i) ;
			sol_hom2_eta.set(lz, k, j, i) = sol(i+nr) ;
			sol_hom2_w.set(lz, k, j, i) = sol(i+2*nr) ;
		    }			
		    sec.set(ind1+nr-1) = 0 ;
		    sec.set(ind2+nr-1) = 1 ;
		    sol = oper.inverse(sec) ;
		    for (int i=0; i<nr; i++) {
			sol_hom3_hrr.set(lz, k, j, i) = sol(i) ;
			sol_hom3_eta.set(lz, k, j, i) = sol(i+nr) ;
			sol_hom3_w.set(lz, k, j, i) = sol(i+2*nr) ;
		    }	
		}
	    }
	}
    }

    //------------------------------
    // Compactified external domain
    //------------------------------
    if (cedbc) {int lz = nzm1 ;  
    if (need_calculation && (par_mat != 0x0)) mat_done.annule_hard() ;
    int nr = mgrid.get_nr(lz) ;
    int ind0 = 0 ;
    int ind1 = nr ;
    int ind2 = 2*nr ;
    double alpha = mp_aff->get_alpha()[lz] ;
    Matrice ope(3*nr, 3*nr) ;
	
    for (int k=0 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    // quantic numbers and spectral bases
	    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)) {
		if (need_calculation) {
		ope.set_etat_qcq() ;
		Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		Diff_sx os(base_r, nr) ; const Matrice& ms = os.get_matrice() ;

		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col) = - md(lin,col) + 3*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col+nr) = -l_q*(l_q+1)*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col+2*nr) = 0 ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col) = -0.5*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col+nr) = -md(lin,col) + 3*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col+2*nr) = (2. - l_q*(l_q+1))*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+2*nr,col) =  0.5*md(lin,col)/double(l_q+1) 
			    - 0.5*double(l_q+4)/double(l_q+1)*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+2*nr,col+nr) = -2*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+2*nr,col+2*nr) =  -(l_q+2)*md(lin,col) 
			    + l_q*(l_q+2)*ms(lin,col) ;
		ope *= 1./alpha ;
		for (int col=0; col<3*nr; col++) {
		    ope.set(ind0+nr-2, col) = 0 ;
		    ope.set(ind0+nr-1, col) = 0 ;
		    ope.set(ind1+nr-2, col) = 0 ;
		    ope.set(ind1+nr-1, col) = 0 ;
		    ope.set(ind2+nr-1, col) = 0 ;
		}
		for (int col=0; col<nr; col++) {
		    ope.set(ind0+nr-1, col+ind0) = 1 ;
		    ope.set(ind1+nr-1, col+ind1) = 1 ;
		    ope.set(ind2+nr-1, col+ind2) = 1 ;
		}
		ope.set(ind0+nr-2, ind0+1) = 1 ;
		ope.set(ind1+nr-2, ind1+2) = 1 ;

		ope.set_lu() ;
		if ((par_mat != 0x0) && (mat_done(l_q) == 0)) {
		    Matrice* pope = new Matrice(ope) ;
		    par_mat->add_matrice_mod(*pope, lz*lmax + l_q) ;
		    mat_done.set(l_q) = 1 ;
		}
		} //End of case when a calculation is needed
		const Matrice& oper = (par_mat == 0x0 ? ope : 
				       par_mat->get_matrice_mod(lz*lmax + l_q) ) ;

		Tbl sec(3*nr) ;
		sec.set_etat_qcq() ;
		if (hnull) {
		    for (int lin=0; lin<2*nr; lin++)
			sec.set(lin) = 0 ;
		    for (int lin=0; lin<nr; lin++)
			sec.set(2*nr+lin) = (*source.get_spectral_va().c_cf)
			    (lz, k, j, lin) ;
		}
		else {
		    for (int lin=0; lin<nr; lin++)
			sec.set(lin) = (*hoverr.get_spectral_va().c_cf)(lz, k, j, lin) ;
		    for (int lin=0; lin<nr; lin++)
			sec.set(lin+nr) = -0.5*(*hoverr.get_spectral_va().c_cf)
			    (lz, k, j, lin) ;
		    if (bnull) {
		    for (int lin=0; lin<nr; lin++)
			sec.set(2*nr+lin) = -0.5/double(l_q+1)*(
			    (*dhdr.get_spectral_va().c_cf)(lz, k, j, lin)
			    + (l_q+2)*(*hoverr.get_spectral_va().c_cf)(lz, k, j, lin) ) ;
		    }
		    else {
		    for (int lin=0; lin<nr; lin++)
			sec.set(2*nr+lin) = -0.5/double(l_q+1)*(
			    (*dhdr.get_spectral_va().c_cf)(lz, k, j, lin)
			    + (l_q+2)*(*hoverr.get_spectral_va().c_cf)(lz, k, j, lin) )
			    + (*source.get_spectral_va().c_cf)(lz, k, j, lin) ;
		    }
		}
		sec.set(ind0+nr-2) = 0 ;
		sec.set(ind0+nr-1) = 0 ;
		sec.set(ind1+nr-1) = 0 ;
		sec.set(ind1+nr-2) = 0 ;
		sec.set(ind2+nr-1) = 0 ;
		Tbl sol = oper.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_part_hrr.set(lz, k, j, i) = sol(i) ;
		    sol_part_eta.set(lz, k, j, i) = sol(i+nr) ;
		    sol_part_w.set(lz, k, j, i) = sol(i+2*nr) ;
		}
		sec.annule_hard() ;
		sec.set(ind0+nr-2) = 1 ;
		sol = oper.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_hom1_hrr.set(lz, k, j, i) = sol(i) ;
		    sol_hom1_eta.set(lz, k, j, i) = sol(i+nr) ;
		    sol_hom1_w.set(lz, k, j, i) = sol(i+2*nr) ;
		}			
		sec.set(ind0+nr-2) = 0 ;
		sec.set(ind1+nr-2) = 1 ;
		sol = oper.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_hom2_hrr.set(lz, k, j, i) = sol(i) ;
		    sol_hom2_eta.set(lz, k, j, i) = sol(i+nr) ;
		    sol_hom2_w.set(lz, k, j, i) = sol(i+2*nr) ;
		}
	    }
	}
    }
    }

    int taille = 3*nz_bc + 1 ;
    if (cedbc) taille-- ;
    Mtbl_cf& mhrr = *hrr.set_spectral_va().c_cf ;
    Mtbl_cf& meta = *tilde_eta.set_spectral_va().c_cf ;
    Mtbl_cf& mw = *ww.set_spectral_va().c_cf ;
	
    Tbl sec_membre(taille) ; 
    Matrice systeme(taille, taille) ; 
    int ligne ;  int colonne ;
    Tbl pipo(1) ;
    const Tbl& hrrb = (cedbc ? pipo : par_bc->get_tbl_mod(1) );
    double chrr = (cedbc ? 0 : par_bc->get_tbl_mod()(4) ) ;
    double dhrr = (cedbc ? 0 : par_bc->get_tbl_mod()(5) ) ;
    double ceta = (cedbc ? 0 : par_bc->get_tbl_mod()(6) ) ;
    double deta = (cedbc ? 0 : par_bc->get_tbl_mod()(7) ) ;
    double cw = (cedbc ? 0 : par_bc->get_tbl_mod()(8) ) ;
    double dw = (cedbc ? 0 : par_bc->get_tbl_mod()(9) ) ;
    Mtbl_cf dhom1_hrr = sol_hom1_hrr ; 
    Mtbl_cf dhom2_hrr = sol_hom2_hrr ; 
    Mtbl_cf dhom3_hrr = sol_hom3_hrr ; 
    Mtbl_cf dpart_hrr = sol_part_hrr ; 
    Mtbl_cf dhom1_eta = sol_hom1_eta ; 
    Mtbl_cf dhom2_eta = sol_hom2_eta ; 
    Mtbl_cf dhom3_eta = sol_hom3_eta ; 
    Mtbl_cf dpart_eta = sol_part_eta ; 
    Mtbl_cf dhom1_w = sol_hom1_w ; 
    Mtbl_cf dhom2_w = sol_hom2_w ; 
    Mtbl_cf dhom3_w = sol_hom3_w ; 
    Mtbl_cf dpart_w = sol_part_w ; 
    if (!cedbc) {
	dhom1_hrr.dsdx() ; dhom1_eta.dsdx() ; dhom1_w.dsdx() ;
	dhom2_hrr.dsdx() ; dhom2_eta.dsdx() ; dhom2_w.dsdx() ;
	dhom3_hrr.dsdx() ; dhom3_eta.dsdx() ; dhom3_w.dsdx() ;
	dpart_hrr.dsdx() ; dpart_eta.dsdx() ; dpart_w.dsdx() ;
    }
    // Loop on l and m
    //----------------
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	    if ((nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)) {
		ligne = 0 ;
		colonne = 0 ;
		systeme.annule_hard() ;
		sec_membre.annule_hard() ;

		//Nucleus 
		int nr = mgrid.get_nr(0) ;
		if (nz_bc >0 ) {
		systeme.set(ligne, colonne) = sol_hom3_hrr.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_hrr.val_out_bound_jk(0, j, k) ;
		ligne++ ;

		systeme.set(ligne, colonne) = sol_hom3_eta.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_eta.val_out_bound_jk(0, j, k) ;
		ligne++ ;

		systeme.set(ligne, colonne) = sol_hom3_w.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_w.val_out_bound_jk(0, j, k) ;
		colonne++ ;
		}
		//shells
		for (int zone=1 ; zone<nz_bc ; zone++) {
		    nr = mgrid.get_nr(zone) ;
		    ligne -= 2 ;

		    //Condition at x = -1
		    systeme.set(ligne, colonne) = 
			- sol_hom1_hrr.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			- sol_hom2_hrr.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+2) = 
			- sol_hom3_hrr.val_in_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) += sol_part_hrr.val_in_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			- sol_hom1_eta.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			- sol_hom2_eta.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+2) = 
			- sol_hom3_eta.val_in_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) += sol_part_eta.val_in_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			- sol_hom1_w.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			- sol_hom2_w.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+2) = 
			- sol_hom3_w.val_in_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) += sol_part_w.val_in_bound_jk(zone, j, k) ;
		    ligne++ ;

		    // Condition at x=1
		    systeme.set(ligne, colonne) = 
			sol_hom1_hrr.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			sol_hom2_hrr.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+2) = 
			sol_hom3_hrr.val_out_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) -= sol_part_hrr.val_out_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			sol_hom1_eta.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			sol_hom2_eta.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+2) = 
			sol_hom3_eta.val_out_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) -= sol_part_eta.val_out_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			sol_hom1_w.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			sol_hom2_w.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+2) = 
			sol_hom3_w.val_out_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) -= sol_part_w.val_out_bound_jk(zone, j, k) ;
		    
		    colonne += 3 ;
		}
    
		//Last domain
		nr = mgrid.get_nr(nz_bc) ;
		double alpha = mp_aff->get_alpha()[nz_bc] ;
		if (nz_bc>0) {
		ligne -= 2 ;

		//Condition at x = -1
		systeme.set(ligne, colonne) = 
		    - sol_hom1_hrr.val_in_bound_jk(nz_bc, j, k) ;
		systeme.set(ligne, colonne+1) = 
		    - sol_hom2_hrr.val_in_bound_jk(nz_bc, j, k) ;
		if (!cedbc) systeme.set(ligne, colonne+2) = 
				- sol_hom3_hrr.val_in_bound_jk(nz_bc, j, k) ;

		sec_membre.set(ligne) += sol_part_hrr.val_in_bound_jk(nz_bc, j, k) ;
		ligne++ ;

		systeme.set(ligne, colonne) = 
		    - sol_hom1_eta.val_in_bound_jk(nz_bc, j, k) ;
		systeme.set(ligne, colonne+1) = 
		    - sol_hom2_eta.val_in_bound_jk(nz_bc, j, k) ;
		if (!cedbc) systeme.set(ligne, colonne+2) = 
		    - sol_hom3_eta.val_in_bound_jk(nz_bc, j, k) ;
		
		sec_membre.set(ligne) += sol_part_eta.val_in_bound_jk(nz_bc, j, k) ;
		ligne++ ;
		
		systeme.set(ligne, colonne) = 
		    - sol_hom1_w.val_in_bound_jk(nz_bc, j, k) ;
		systeme.set(ligne, colonne+1) = 
		    - sol_hom2_w.val_in_bound_jk(nz_bc, j, k) ;
		if (!cedbc) systeme.set(ligne, colonne+2) = 
				- sol_hom3_w.val_in_bound_jk(nz_bc, j, k) ;
		
		sec_membre.set(ligne) += sol_part_w.val_in_bound_jk(nz_bc, j, k) ;
		ligne++ ;
		}
		if (!cedbc) {//Special condition at x=1
		    if (nz_bc > 0) {
		systeme.set(ligne, colonne) = 
		    chrr*sol_hom1_hrr.val_out_bound_jk(nz_bc, j, k) 
		    + dhrr*dhom1_hrr.val_out_bound_jk(nz_bc, j, k) / alpha 
		    + ceta*sol_hom1_eta.val_out_bound_jk(nz_bc, j, k) 
		    + deta*dhom1_eta.val_out_bound_jk(nz_bc, j, k) / alpha 
		    + cw*sol_hom1_w.val_out_bound_jk(nz_bc, j, k) 
		    + dw*dhom1_w.val_out_bound_jk(nz_bc, j, k) / alpha ;
		systeme.set(ligne, colonne+1) = 
		    chrr*sol_hom2_hrr.val_out_bound_jk(nz_bc, j, k) 
		    + dhrr*dhom2_hrr.val_out_bound_jk(nz_bc, j, k) / alpha 
		    + ceta*sol_hom2_eta.val_out_bound_jk(nz_bc, j, k) 
		    + deta*dhom2_eta.val_out_bound_jk(nz_bc, j, k) / alpha 
		    + cw*sol_hom2_w.val_out_bound_jk(nz_bc, j, k) 
		    + dw*dhom2_w.val_out_bound_jk(nz_bc, j, k) / alpha ;
		    }
		else { // Only nucleus
		    assert(ligne == 0) ;
		    colonne = -2 ;
		}
		systeme.set(ligne, colonne+2) = 
		    chrr*sol_hom3_hrr.val_out_bound_jk(nz_bc, j, k) 
		    + dhrr*dhom3_hrr.val_out_bound_jk(nz_bc, j, k) / alpha 
		    + ceta*sol_hom3_eta.val_out_bound_jk(nz_bc, j, k) 
		    + deta*dhom3_eta.val_out_bound_jk(nz_bc, j, k) / alpha 
		    + cw*sol_hom3_w.val_out_bound_jk(nz_bc, j, k) 
		    + dw*dhom3_w.val_out_bound_jk(nz_bc, j, k) / alpha ;

		sec_membre.set(ligne) -= chrr*sol_part_hrr.val_out_bound_jk(nz_bc, j, k) 
		    + dhrr*dpart_hrr.val_out_bound_jk(nz_bc, j, k)/alpha 
		    + ceta*sol_part_eta.val_out_bound_jk(nz_bc, j, k) 
		    + deta*dpart_eta.val_out_bound_jk(nz_bc, j, k)/alpha 
		    + cw*sol_part_w.val_out_bound_jk(nz_bc, j, k) 
		    + dw*dpart_w.val_out_bound_jk(nz_bc, j, k)/alpha 
		    - hrrb(k, j)  ;
		}
			
		// Solution of the system giving the coefficients for the homogeneous 
		// solutions
		//-------------------------------------------------------------------
		systeme.set_lu() ;
		Tbl facteur = systeme.inverse(sec_membre) ;
		int conte = 0 ;

		// everything is put to the right place, 
		//---------------------------------------
 		nr = mgrid.get_nr(0) ; //nucleus
 		for (int i=0 ; i<nr ; i++) {
		    mhrr.set(0, k, j, i) = sol_part_hrr(0, k, j, i)
			+ facteur(conte)*sol_hom3_hrr(0, k, j, i) ;
		    meta.set(0, k, j, i) = sol_part_eta(0, k, j, i)
			+ facteur(conte)*sol_hom3_eta(0, k, j, i) ;
		    mw.set(0, k, j, i) = sol_part_w(0, k, j, i)
			+ facteur(conte)*sol_hom3_w(0, k, j, i) ;
 		}
 		conte++ ;
 		for (int zone=1 ; zone<=n_shell ; zone++) { //shells
 		    nr = mgrid.get_nr(zone) ;
 		    for (int i=0 ; i<nr ; i++) {
		    mhrr.set(zone, k, j, i) = sol_part_hrr(zone, k, j, i)
			+ facteur(conte)*sol_hom1_hrr(zone, k, j, i) 
			+ facteur(conte+1)*sol_hom2_hrr(zone, k, j, i) 
			+ facteur(conte+2)*sol_hom3_hrr(zone, k, j, i) ;
			
		    meta.set(zone, k, j, i) = sol_part_eta(zone, k, j, i)
			+ facteur(conte)*sol_hom1_eta(zone, k, j, i) 
			+ facteur(conte+1)*sol_hom2_eta(zone, k, j, i) 
			+ facteur(conte+2)*sol_hom3_eta(zone, k, j, i) ;
			
		    mw.set(zone, k, j, i) = sol_part_w(zone, k, j, i)
			+ facteur(conte)*sol_hom1_w(zone, k, j, i) 
			+ facteur(conte+1)*sol_hom2_w(zone, k, j, i) 
			+ facteur(conte+2)*sol_hom3_w(zone, k, j, i) ;			
 		    }
 		    conte+=3 ;
 		}
		if (cedbc) {
		    nr = mgrid.get_nr(nzm1) ; //compactified external domain
		    for (int i=0 ; i<nr ; i++) {
			mhrr.set(nzm1, k, j, i) = sol_part_hrr(nzm1, k, j, i)
			+ facteur(conte)*sol_hom1_hrr(nzm1, k, j, i) 
			+ facteur(conte+1)*sol_hom2_hrr(nzm1, k, j, i) ;

		    meta.set(nzm1, k, j, i) = sol_part_eta(nzm1, k, j, i)
			+ facteur(conte)*sol_hom1_eta(nzm1, k, j, i) 
			+ facteur(conte+1)*sol_hom2_eta(nzm1, k, j, i) ; 
			
		    mw.set(nzm1, k, j, i) = sol_part_w(nzm1, k, j, i)
			+ facteur(conte)*sol_hom1_w(nzm1, k, j, i) 
			+ facteur(conte+1)*sol_hom2_w(nzm1, k, j, i) ;
		    }
		}
	    } // End of nullite_plm  
	} //End of loop on theta
		    

    if (hrr.set_spectral_va().c != 0x0) 
	delete hrr.set_spectral_va().c ;
    hrr.set_spectral_va().c = 0x0 ;
    hrr.set_spectral_va().ylm_i() ;

    if (tilde_eta.set_spectral_va().c != 0x0) 
	delete tilde_eta.set_spectral_va().c ;
    tilde_eta.set_spectral_va().c = 0x0 ;
    tilde_eta.set_spectral_va().ylm_i() ;

    if (ww.set_spectral_va().c != 0x0) 
	delete ww.set_spectral_va().c ;
    ww.set_spectral_va().c = 0x0 ;
    ww.set_spectral_va().ylm_i() ;

}

void Sym_tensor_trans::sol_Dirac_l01(const Scalar& hh, Scalar& hrr, Scalar& tilde_eta,
				     Param* par_mat) const {

    const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
    assert(mp_aff != 0x0) ; //Only affine mapping for the moment

    if (hh.get_etat() == ETATZERO) {
	hrr.annule_hard() ;
	tilde_eta.annule_hard() ;
	return ;
    }

    const Mg3d& mgrid = *mp_aff->get_mg() ;
    int nz = mgrid.get_nzone() ;
    assert(mgrid.get_type_r(0) == RARE)  ;
    assert(mgrid.get_type_r(nz-1) == UNSURR) ;

    int nt = mgrid.get_nt(0) ;
    int np = mgrid.get_np(0) ;

    Scalar source = hh ;
    source.div_r_dzpuis(2) ;
    Scalar source_coq = hh ;
    source.set_spectral_va().ylm() ;
    source_coq.set_spectral_va().ylm() ;
    Base_val base = source.get_spectral_base() ;
    base.mult_x() ;
    int lmax = base.give_lmax(mgrid, 0) + 1;

    assert (hrr.get_spectral_base() == base) ;
    assert (tilde_eta.get_spectral_base() == base) ;
    assert (hrr.get_spectral_va().c_cf != 0x0) ;
    assert (tilde_eta.get_spectral_va().c_cf != 0x0) ;
 
    Mtbl_cf sol_part_hrr(mgrid, base) ; sol_part_hrr.annule_hard() ;
    Mtbl_cf sol_part_eta(mgrid, base) ; sol_part_eta.annule_hard() ;
    Mtbl_cf sol_hom1_hrr(mgrid, base) ; sol_hom1_hrr.annule_hard() ;
    Mtbl_cf sol_hom1_eta(mgrid, base) ; sol_hom1_eta.annule_hard() ;
    Mtbl_cf sol_hom2_hrr(mgrid, base) ; sol_hom2_hrr.annule_hard() ;
    Mtbl_cf sol_hom2_eta(mgrid, base) ; sol_hom2_eta.annule_hard() ;

    bool need_calculation = true ;
    if (par_mat != 0x0)
	if (par_mat->get_n_matrice_mod() > 0) 
	    if (&par_mat->get_matrice_mod(0) != 0x0) need_calculation = false ;

    int l_q, m_q, base_r ;
    Itbl mat_done(lmax) ;

    //---------------
    //--  NUCLEUS ---
    //---------------
    {int lz = 0 ;  
    int nr = mgrid.get_nr(lz) ;
    double alpha = mp_aff->get_alpha()[lz] ;
    Matrice ope(2*nr, 2*nr) ;
    if (need_calculation && (par_mat != 0x0)) mat_done.annule_hard() ;

    for (int k=0 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    // quantic numbers and spectral bases
	    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q < 2)) {
		if (need_calculation) {
		    ope.set_etat_qcq() ;
		    Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		    Diff_sx os(base_r, nr) ; const Matrice& ms = os.get_matrice() ;

		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col) = md(lin,col) + 3*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col+nr) = -l_q*(l_q+1)*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col) = -0.5*ms(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col+nr) = md(lin,col) + 3*ms(lin, col);

		    ope *= 1./alpha ;
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
		    
		    ope.set_lu() ;
		    if ((par_mat != 0x0) && (mat_done(l_q) == 0)) {
			Matrice* pope = new Matrice(ope) ;
			par_mat->add_matrice_mod(*pope, lz*lmax + l_q) ;
			mat_done.set(l_q) = 1 ;
		    }
		} //End of case when a calculation is needed

		const Matrice& oper = (par_mat == 0x0 ? ope : 
				       par_mat->get_matrice_mod(lz*lmax + l_q) ) ;
		Tbl sec(2*nr) ;
		sec.set_etat_qcq() ;
		for (int lin=0; lin<nr; lin++)
		    sec.set(lin) = (*source.get_spectral_va().c_cf)(lz, k, j, lin) ;
		for (int lin=0; lin<nr; lin++)
		    sec.set(nr+lin) = -0.5*(*source.get_spectral_va().c_cf)
			(lz, k, j, lin) ;
		sec.set(nr-1) = 0 ;
		if (base_r == R_CHEBP) {
		    double h0 = 0 ; //In the l=0 case:  3*hrr(r=0) = h(r=0) 
		    int pari = 1 ;
		    for (int col=0; col<nr; col++) {
			h0 += pari*
			    (*source_coq.get_spectral_va().c_cf)(lz, k, j, col) ;
			pari = - pari ;
		    }
		    sec.set(nr-1) = h0 / 3. ;
		}
		sec.set(2*nr-1) = 0 ;
		Tbl sol = oper.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_part_hrr.set(lz, k, j, i) = sol(i) ;
		    sol_part_eta.set(lz, k, j, i) = sol(i+nr) ;
		}
		sec.annule_hard() ;
	    }
	}
    }
    }

    //-------------
    // -- Shells --
    //-------------

    for (int lz=1; lz<nz-1; lz++) {
	if (need_calculation && (par_mat != 0x0)) mat_done.annule_hard() ;
	int nr = mgrid.get_nr(lz) ;
	int ind0 = 0 ;
	int ind1 = nr ;
	assert(mgrid.get_nt(lz) == nt) ;
	assert(mgrid.get_np(lz) == np) ;
	double alpha = mp_aff->get_alpha()[lz] ;
	double ech = mp_aff->get_beta()[lz] / alpha ;
	Matrice ope(2*nr, 2*nr) ;
	
	for (int k=0 ; k<np+1 ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		// quantic numbers and spectral bases
		base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q < 2)) {
		    if (need_calculation) {
		    ope.set_etat_qcq() ;
		    Diff_xdsdx oxd(base_r, nr) ; const Matrice& mxd = oxd.get_matrice() ;
		    Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		    Diff_id oid(base_r, nr) ; const Matrice& mid = oid.get_matrice() ;

		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col) = mxd(lin,col) + ech*md(lin,col) 
				+ 3*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin,col+nr) = -l_q*(l_q+1)*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col) = -0.5*mid(lin,col) ;
		    for (int lin=0; lin<nr; lin++) 
			for (int col=0; col<nr; col++) 
			    ope.set(lin+nr,col+nr) = mxd(lin,col) + ech*md(lin,col) 
				+ 3*mid(lin, col) ;

		    for (int col=0; col<2*nr; col++) {
			ope.set(ind0+nr-1, col) = 0 ;
			ope.set(ind1+nr-1, col) = 0 ;
		    }
		    ope.set(ind0+nr-1, ind0) = 1 ;
		    ope.set(ind1+nr-1, ind1) = 1 ;

		    ope.set_lu() ;
		    if ((par_mat != 0x0) && (mat_done(l_q) == 0)) {
			Matrice* pope = new Matrice(ope) ;
			par_mat->add_matrice_mod(*pope, lz*lmax + l_q) ;
			mat_done.set(l_q) = 1 ;
		    }
		    } //End of case when a calculation is needed
		    const Matrice& oper = (par_mat == 0x0 ? ope : 
				       par_mat->get_matrice_mod(lz*lmax + l_q) ) ;
		    Tbl sec(2*nr) ;
		    sec.set_etat_qcq() ;
		    for (int lin=0; lin<nr; lin++)
			sec.set(lin) = (*source_coq.get_spectral_va().c_cf)
			    (lz, k, j, lin) ; 
		    for (int lin=0; lin<nr; lin++)
			sec.set(nr+lin) = -0.5*(*source_coq.get_spectral_va().c_cf)
			    (lz, k, j, lin) ;
		    sec.set(ind0+nr-1) = 0 ;
		    sec.set(ind1+nr-1) = 0 ;
		    Tbl sol = oper.inverse(sec) ;

		    for (int i=0; i<nr; i++) {
 			sol_part_hrr.set(lz, k, j, i) = sol(i) ;
 			sol_part_eta.set(lz, k, j, i) = sol(i+nr) ;
		    }
		    sec.annule_hard() ;
		    sec.set(ind0+nr-1) = 1 ;
		    sol = oper.inverse(sec) ;
		    for (int i=0; i<nr; i++) {
			sol_hom1_hrr.set(lz, k, j, i) = sol(i) ;
			sol_hom1_eta.set(lz, k, j, i) = sol(i+nr) ;
		    }			
		    sec.set(ind0+nr-1) = 0 ;
		    sec.set(ind1+nr-1) = 1 ;
		    sol = oper.inverse(sec) ;
		    for (int i=0; i<nr; i++) {
			sol_hom2_hrr.set(lz, k, j, i) = sol(i) ;
			sol_hom2_eta.set(lz, k, j, i) = sol(i+nr) ;
		    }			
		}
	    }
	}
    }

    //------------------------------
    // Compactified external domain
    //------------------------------
    {int lz = nz-1 ;  
    if (need_calculation && (par_mat != 0x0)) mat_done.annule_hard() ;
    int nr = mgrid.get_nr(lz) ;
    int ind0 = 0 ;
    int ind1 = nr ;
    assert(mgrid.get_nt(lz) == nt) ;
    assert(mgrid.get_np(lz) == np) ;
    double alpha = mp_aff->get_alpha()[lz] ;
    Matrice ope(2*nr, 2*nr) ;
	
    for (int k=0 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    // quantic numbers and spectral bases
	    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q < 2)) {
		if (need_calculation) {
		ope.set_etat_qcq() ;
		Diff_dsdx od(base_r, nr) ; const Matrice& md = od.get_matrice() ;
		Diff_sx os(base_r, nr) ; const Matrice& ms = os.get_matrice() ;

		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col) = - md(lin,col) + 3*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin,col+nr) = -l_q*(l_q+1)*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col) = -0.5*ms(lin,col) ;
		for (int lin=0; lin<nr; lin++) 
		    for (int col=0; col<nr; col++) 
			ope.set(lin+nr,col+nr) = -md(lin,col) + 3*ms(lin, col) ;

		ope *= 1./alpha ;
		for (int col=0; col<2*nr; col++) {
		    ope.set(ind0+nr-2, col) = 0 ;
		    ope.set(ind0+nr-1, col) = 0 ;
		    ope.set(ind1+nr-2, col) = 0 ;
		    ope.set(ind1+nr-1, col) = 0 ;
		}
		for (int col=0; col<nr; col++) {
		    ope.set(ind0+nr-1, col+ind0) = 1 ;
		    ope.set(ind1+nr-1, col+ind1) = 1 ;
		}
		ope.set(ind0+nr-2, ind0+1) = 1 ;
		ope.set(ind1+nr-2, ind1+1) = 1 ;

		ope.set_lu() ;
		if ((par_mat != 0x0) && (mat_done(l_q) == 0)) {
		    Matrice* pope = new Matrice(ope) ;
		    par_mat->add_matrice_mod(*pope, lz*lmax + l_q) ;
		    mat_done.set(l_q) = 1 ;
		}
		} //End of case when a calculation is needed
		const Matrice& oper = (par_mat == 0x0 ? ope : 
				       par_mat->get_matrice_mod(lz*lmax + l_q) ) ;
		Tbl sec(2*nr) ;
		sec.set_etat_qcq() ;
		for (int lin=0; lin<nr; lin++)
		    sec.set(lin) = (*source.get_spectral_va().c_cf)
			(lz, k, j, lin) ;
		for (int lin=0; lin<nr; lin++)
		    sec.set(nr+lin) = -0.5*(*source.get_spectral_va().c_cf)
			(lz, k, j, lin) ;
		sec.set(ind0+nr-2) = 0 ;
		sec.set(ind0+nr-1) = 0 ;
		sec.set(ind1+nr-2) = 0 ;
		sec.set(ind1+nr-1) = 0 ;
 		Tbl sol = oper.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_part_hrr.set(lz, k, j, i) = sol(i) ;
		    sol_part_eta.set(lz, k, j, i) = sol(i+nr) ;
		}
		sec.annule_hard() ;
		sec.set(ind0+nr-2) = 1 ;
		sol = oper.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_hom1_hrr.set(lz, k, j, i) = sol(i) ;
		    sol_hom1_eta.set(lz, k, j, i) = sol(i+nr) ;
		}
		sec.set(ind0+nr-2) = 0 ;
		sec.set(ind1+nr-2) = 1 ;
		sol = oper.inverse(sec) ;
		for (int i=0; i<nr; i++) {
		    sol_hom2_hrr.set(lz, k, j, i) = sol(i) ;
		    sol_hom2_eta.set(lz, k, j, i) = sol(i+nr) ;
		}
	    }
	}
    }
    }

    int taille = 2*(nz-1) ;
    Mtbl_cf& mhrr = *hrr.set_spectral_va().c_cf ;
    Mtbl_cf& meta = *tilde_eta.set_spectral_va().c_cf ;
	
    Tbl sec_membre(taille) ; 
    Matrice systeme(taille, taille) ; 
    int ligne ;  int colonne ;
	
    // Loop on l and m
    //----------------
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	    if ((nullite_plm(j, nt, k, np, base) == 1) && (l_q < 2)) {
		ligne = 0 ;
		colonne = 0 ;
		systeme.annule_hard() ;
		sec_membre.annule_hard() ;

		//Nucleus 
		int nr = mgrid.get_nr(0) ;
		
		sec_membre.set(ligne) = -sol_part_hrr.val_out_bound_jk(0, j, k) ;
		ligne++ ;

		sec_membre.set(ligne) = -sol_part_eta.val_out_bound_jk(0, j, k) ;

		//shells
		for (int zone=1 ; zone<nz-1 ; zone++) {
		    nr = mgrid.get_nr(zone) ;
		    ligne-- ;

		    //Condition at x = -1
		    systeme.set(ligne, colonne) = 
			- sol_hom1_hrr.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			- sol_hom2_hrr.val_in_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) += sol_part_hrr.val_in_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			- sol_hom1_eta.val_in_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			- sol_hom2_eta.val_in_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) += sol_part_eta.val_in_bound_jk(zone, j, k) ;
		    ligne++ ;

		    // Condition at x=1
		    systeme.set(ligne, colonne) = 
			sol_hom1_hrr.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			sol_hom2_hrr.val_out_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) -= sol_part_hrr.val_out_bound_jk(zone, j, k) ;
		    ligne++ ;

		    systeme.set(ligne, colonne) = 
			sol_hom1_eta.val_out_bound_jk(zone, j, k) ;
		    systeme.set(ligne, colonne+1) = 
			sol_hom2_eta.val_out_bound_jk(zone, j, k) ;

		    sec_membre.set(ligne) -= sol_part_eta.val_out_bound_jk(zone, j, k) ;
		    
		    colonne += 2 ;
		}
    
		//Compactified external domain
		nr = mgrid.get_nr(nz-1) ;

		ligne-- ;

		systeme.set(ligne, colonne) = 
		    - sol_hom1_hrr.val_in_bound_jk(nz-1, j, k) ;
		systeme.set(ligne, colonne+1) = 
		    - sol_hom2_hrr.val_in_bound_jk(nz-1, j, k) ;
		
		sec_membre.set(ligne) += sol_part_hrr.val_in_bound_jk(nz-1, j, k) ;
		ligne++ ;
		
		systeme.set(ligne, colonne) = 
		    - sol_hom1_eta.val_in_bound_jk(nz-1, j, k) ;
		systeme.set(ligne, colonne+1) = 
		    - sol_hom2_eta.val_in_bound_jk(nz-1, j, k) ;
		
		sec_membre.set(ligne) += sol_part_eta.val_in_bound_jk(nz-1, j, k) ;
			
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
		    mhrr.set(0, k, j, i) = sol_part_hrr(0, k, j, i) ;
		    meta.set(0, k, j, i) = sol_part_eta(0, k, j, i) ;
 		}
 		for (int zone=1 ; zone<nz-1 ; zone++) { //shells
 		    nr = mgrid.get_nr(zone) ;
 		    for (int i=0 ; i<nr ; i++) {
		    mhrr.set(zone, k, j, i) = sol_part_hrr(zone, k, j, i)
			+ facteur(conte)*sol_hom1_hrr(zone, k, j, i) 
			+ facteur(conte+1)*sol_hom2_hrr(zone, k, j, i) ;
			
		    meta.set(zone, k, j, i) = sol_part_eta(zone, k, j, i)
			+ facteur(conte)*sol_hom1_eta(zone, k, j, i) 
			+ facteur(conte+1)*sol_hom2_eta(zone, k, j, i) ;
 		    }
 		    conte+=2 ;
 		}
 		nr = mgrid.get_nr(nz-1) ; //compactified external domain
 		for (int i=0 ; i<nr ; i++) {
		    mhrr.set(nz-1, k, j, i) = sol_part_hrr(nz-1, k, j, i)
			+ facteur(conte)*sol_hom1_hrr(nz-1, k, j, i) 
			+ facteur(conte+1)*sol_hom2_hrr(nz-1, k, j, i) ;
			
		    meta.set(nz-1, k, j, i) = sol_part_eta(nz-1, k, j, i)
			+ facteur(conte)*sol_hom1_eta(nz-1, k, j, i) 
			+ facteur(conte+1)*sol_hom2_eta(nz-1, k, j, i) ;
		}

	    } // End of nullite_plm  
	} //End of loop on theta
}

}
