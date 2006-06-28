/*
 *  Solution of the Poisson the 2nd-order PDE for hrr/eta (2 transverse conditions)
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005  Jerome Novak
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

char sym_tensor_trans_pde_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2006/06/28 07:48:26  j_novak
 * Better treatment of some null cases.
 *
 * Revision 1.9  2006/06/21 15:42:47  j_novak
 * Minor changes.
 *
 * Revision 1.8  2006/06/20 12:07:15  j_novak
 * Improved execution speed for sol_Dirac_tildeB...
 *
 * Revision 1.7  2006/06/14 10:04:21  j_novak
 * New methods sol_Dirac_l01, set_AtB_det_one and set_AtB_trace_zero.
 *
 * Revision 1.6  2006/06/13 13:30:12  j_novak
 * New members sol_Dirac_A and sol_Dirac_tildeB (see documentation).
 *
 * Revision 1.5  2006/06/12 13:37:23  j_novak
 * Added bounds in l (multipolar momentum) for Sym_tensor_trans::solve_hrr.
 *
 * Revision 1.4  2005/11/28 14:45:17  j_novak
 * Improved solution of the Poisson tensor equation in the case of a transverse
 * tensor.
 *
 * Revision 1.3  2005/11/24 14:07:54  j_novak
 * Use of Matrice::annule_hard()
 *
 * Revision 1.2  2005/11/24 09:24:25  j_novak
 * Corrected some missing references.
 *
 * Revision 1.1  2005/09/16 13:58:11  j_novak
 * New Poisson solver for a Sym_tensor_trans.
 *
 *
 * $Heade$
 *
 */

// C headers
#include <assert.h>
#include <math.h>

// Lorene headers
#include "tensor.h"
#include "diff.h"
#include "proto.h"
#include "param.h"

Sym_tensor_trans Sym_tensor_trans::poisson(const Scalar* h_guess) const {

    // All this has a meaning only for spherical components...
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    //## ... and affine mapping, for the moment!
    const Map_af* mpaff = dynamic_cast<const Map_af*>(mp) ;
    assert( mpaff!= 0x0) ;
    Sym_tensor_trans resu(*mp, *triad, *met_div) ;

    const Mg3d& gri = *mp->get_mg() ;
    int np = gri.get_np(0) ;
    int nt = gri.get_nt(0) ;
    assert (nt > 4) ;
    if (np == 1) {
	int nz = gri.get_nzone() ;
	double* bornes = new double[nz+1] ;
	const double* alp = mpaff->get_alpha() ;
	const double* bet = mpaff->get_beta() ; 
	for (int lz=0; lz<nz; lz++) {
	    assert (gri.get_np(lz) == np) ;
	    assert (gri.get_nt(lz) == nt) ;
	    switch (gri.get_type_r(lz)) {
		case RARE: {
		    bornes[lz] = bet[lz] ;
		    break ;
		}
		case FIN: {
		    bornes[lz] = bet[lz] - alp[lz] ;
		    break ;
		}
		case UNSURR: {
		    bornes[lz] = double(1) / ( bet[lz] - alp[lz] ) ;
		    break ;
		}
		default: {
		    cout << "Sym_tensor_trans::poisson() : problem with the grid!" 
			 << endl ;
		    abort() ;
		    break ;
		}
	    }
	}
	if (gri.get_type_r(nz-1) == UNSURR) 
	    bornes[nz] = 1./(alp[nz-1] + bet[nz-1]) ;
	else
	    bornes[nz] = alp[nz-1] + bet[nz-1] ;
	
	const Mg3d& gr2 = *gri.get_non_axi() ;
	Map_af mp2(gr2, bornes) ;
	int np2 = ( np > 3 ? np : 4 ) ;
	
	Sym_tensor sou_cart(mp2, CON, mp2.get_bvect_spher()) ;
	for (int l=1; l<=3; l++)
	    for (int c=l; c<=3; c++) {
		switch (this->operator()(l,c).get_etat() ) {
		    case ETATZERO: {
			sou_cart.set(l,c).set_etat_zero() ;
			break ;
		    }
		    case ETATUN: {
			sou_cart.set(l,c).set_etat_one() ;
			break ;
		    }
		    case ETATQCQ : {
			sou_cart.set(l,c).allocate_all() ;
			for (int lz=0; lz<nz; lz++) 			
			    for (int k=0; k<np2; k++)
				for (int j=0; j<nt; j++)
				    for(int i=0; i<gr2.get_nr(lz); i++)
					sou_cart.set(l,c).set_grid_point(lz, k, j, i)
				   = this->operator()(l,c).val_grid_point(lz, 0, j, i) ;
			break ;
		    }
		    default: {
			cout << 
			    "Sym_tensor_trans::poisson() : source in undefined state!" 
			     << endl ;
			abort() ;
			break ; 
		    }
		}
		sou_cart.set(l,c).set_dzpuis(this->operator()(l,c).get_dzpuis()) ;
	    }
	sou_cart.std_spectral_base() ;
	sou_cart.change_triad(mp2.get_bvect_cart()) ;
	Sym_tensor res_cart(mp2, CON, mp2.get_bvect_cart()) ;
	for (int i=1; i<=3; i++)
	    for(int j=i; j<=3; j++) 
		res_cart.set(i,j) = sou_cart(i,j).poisson() ;
	res_cart.change_triad(mp2.get_bvect_spher()) ;
	Scalar res_x(*mp) ;
	Scalar res_w(*mp) ;

 	switch (res_cart.xxx().get_etat() ) {
 	    case ETATZERO: {
		res_x.set_etat_zero() ;
 		break ;
 	    }
 	    case ETATUN : {
 		res_x.set_etat_one() ;
 		break ;
 	    }
 	    case ETATQCQ : {
 		res_x.allocate_all() ;
 		for (int lz=0; lz<nz; lz++) 			
 		    for (int k=0; k<np; k++)
 			for (int j=0; j<nt; j++)
 			    for(int i=0; i<gri.get_nr(lz); i++)
 				res_x.set_grid_point(lz, k, j, i)
 				    = res_cart.xxx().val_grid_point(lz, k, j, i) ;
 		break ;
 	    }
 	    default: {
 		cout << 
 		    "Sym_tensor_trans::poisson() : res_x in undefined state!" 
 		     << endl ;
 		abort() ;
 		break ; 
 	    }
 	}
	res_x.set_spectral_base(res_cart.xxx().get_spectral_base()) ;

	switch (res_cart.www().get_etat() ) {
	    case ETATZERO: {
		res_w.set_etat_zero() ;
		break ;
	    }
	    case ETATUN : {
		res_w.set_etat_one() ;
		break ;
	    }
	    case ETATQCQ : {
		res_w.allocate_all() ;
		for (int lz=0; lz<nz; lz++) 			
		    for (int k=0; k<np; k++)
			for (int j=0; j<nt; j++)
			    for(int i=0; i<gri.get_nr(lz); i++)
				res_w.set_grid_point(lz, k, j, i)
				    = res_cart.www().val_grid_point(lz, k, j, i) ;
		break ;
	    }
	    default: {
		cout << 
		    "Sym_tensor_trans::poisson() : res_w in undefined state!" 
		     << endl ;
		abort() ;
		break ; 
	    }
	}
	res_w.set_spectral_base(res_cart.www().get_spectral_base()) ;
	
	resu.set_WX_det_one(res_w, res_x, h_guess) ;
	
	delete [] bornes ;
    }
    else {
	assert (np >=4) ;
	Sym_tensor_trans sou_cart = *this ;
	sou_cart.change_triad(mp->get_bvect_cart()) ;
	
	Sym_tensor res_cart(*mp, CON, mp->get_bvect_cart()) ;
	for (int i=1; i<=3; i++)
	    for(int j=i; j<=3; j++) 
		res_cart.set(i,j) = sou_cart(i,j).poisson() ;

	res_cart.change_triad(*triad) ;
	
	resu.set_WX_det_one(res_cart.www(), res_cart.xxx(), h_guess) ;
	
    }
    Vector dive = resu.divergence(*met_div) ;
    dive.dec_dzpuis(2) ;

    maxabs(dive, "Sym_tensor_trans::poisson : divergence of the solution") ;
    
    return resu ;   
}

void Sym_tensor_trans::sol_Dirac_A(const Scalar& aaa, Scalar& tilde_mu, Scalar& x_new) 
    const {

    const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
    assert(mp_aff != 0x0) ; //Only affine mapping for the moment

    const Mg3d& mgrid = *mp_aff->get_mg() ;
    int nz = mgrid.get_nzone() ;
    assert(mgrid.get_type_r(0) == RARE)  ;
    assert(mgrid.get_type_r(nz-1) == UNSURR) ;
    if (aaa.get_etat() == ETATZERO) {
	tilde_mu = 0 ;
	x_new = 0 ;
	return ;
    }
    assert(aaa.get_etat() != ETATNONDEF) ;

    int nt = mgrid.get_nt(0) ;
    int np = mgrid.get_np(0) ;

    Scalar source = aaa ;
    Scalar source_coq = aaa ;
    source_coq.annule_domain(0) ;
    source_coq.annule_domain(nz-1) ;
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
	    donne_lm(nz, lz, j, k, base, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
	    {
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

    for (int lz=1; lz<nz-1; lz++) {
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
		donne_lm(nz, lz, j, k, base, m_q, l_q, base_r) ;
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
		{
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
    {int lz = nz-1 ;  
    int nr = mgrid.get_nr(lz) ;
    assert(mgrid.get_nt(lz) == nt) ;
    assert(mgrid.get_np(lz) == np) ;
    double alpha = mp_aff->get_alpha()[lz] ;
    Matrice ope(2*nr, 2*nr) ;
    ope.set_etat_qcq() ;
	
    for (int k=0 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    // quantic numbers and spectral bases
	    donne_lm(nz, lz, j, k, base, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
	    {
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

    int taille = 2*(nz-1) ;
    Mtbl_cf& mmu = *tilde_mu.set_spectral_va().c_cf ;
    Mtbl_cf& mw = *x_new.set_spectral_va().c_cf ;
	
    Tbl sec_membre(taille) ; 
    Matrice systeme(taille, taille) ; 
    int ligne ;  int colonne ;
	
    // Loop on l and m
    //----------------
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	    if ((nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)){
		ligne = 0 ;
		colonne = 0 ;
		systeme.annule_hard() ;
		sec_membre.annule_hard() ;

		//Nucleus 
		int nr = mgrid.get_nr(0) ;
		
		systeme.set(ligne, colonne) = sol_hom2_mu.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_mu.val_out_bound_jk(0, j, k) ;
		ligne++ ;

		systeme.set(ligne, colonne) = sol_hom2_x.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_x.val_out_bound_jk(0, j, k) ;
		colonne++ ;

		//shells
		for (int zone=1 ; zone<nz-1 ; zone++) {
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
    
		//Compactified external domain
		nr = mgrid.get_nr(nz-1) ;

		ligne-- ;

		systeme.set(ligne, colonne) = 
		    - sol_hom1_mu.val_in_bound_jk(nz-1, j, k) ;
		
		sec_membre.set(ligne) += sol_part_mu.val_in_bound_jk(nz-1, j, k) ;
		ligne++ ;
		
		systeme.set(ligne, colonne) = 
		    - sol_hom1_x.val_in_bound_jk(nz-1, j, k) ;
		
		sec_membre.set(ligne) += sol_part_x.val_in_bound_jk(nz-1, j, k) ;
			
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
 		for (int zone=1 ; zone<nz-1 ; zone++) { //shells
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
 		nr = mgrid.get_nr(nz-1) ; //compactified external domain
 		for (int i=0 ; i<nr ; i++) {
		    mmu.set(nz-1, k, j, i) = sol_part_mu(nz-1, k, j, i)
			+ facteur(conte)*sol_hom1_mu(nz-1, k, j, i) ;
			
		    mw.set(nz-1, k, j, i) = sol_part_x(nz-1, k, j, i)
			+ facteur(conte)*sol_hom1_x(nz-1, k, j, i) ;
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

void Sym_tensor_trans::sol_Dirac_tilde_B(const Scalar& tilde_b, const Scalar& hh, 
					 Scalar& hrr, Scalar& tilde_eta, Scalar& ww,
					 Param* par) 
    const {

    const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
    assert(mp_aff != 0x0) ; //Only affine mapping for the moment

    const Mg3d& mgrid = *mp_aff->get_mg() ;
    int nz = mgrid.get_nzone() ;
    assert(mgrid.get_type_r(0) == RARE)  ;
    assert(mgrid.get_type_r(nz-1) == UNSURR) ;
    if ( (tilde_b.get_etat() == ETATZERO) && (hh.get_etat() == ETATZERO) ) {
	hrr = 0 ;
	tilde_eta = 0 ;
	ww = 0 ;
	return ;
    }
    int nt = mgrid.get_nt(0) ;
    int np = mgrid.get_np(0) ;

    assert (tilde_b.get_etat() != ETATNONDEF) ;
    assert (hh.get_etat() != ETATNONDEF) ;

    Scalar source = tilde_b ;
    Scalar source_coq = tilde_b ;
    source_coq.annule_domain(0) ;
    source_coq.annule_domain(nz-1) ;
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
    if (par != 0x0) {
	bool param_new = false ;
	if ((par->get_n_int_mod() >= 4)
	    &&(par->get_n_tbl_mod()>=1)
	    &&(par->get_n_matrice_mod()>=1)
	    &&(par->get_n_itbl_mod()>=1)) {
	    if (par->get_int_mod(0) != nz) param_new = true ;
	    if (par->get_int_mod(1) != lmax) param_new = true ;
	    if (par->get_int_mod(2) != mgrid.get_type_t() ) param_new = true ;
	    if (par->get_int_mod(3) != mgrid.get_type_p() ) param_new = true ;
	    if (par->get_itbl_mod(0)(0) != mgrid.get_nr(0)) param_new = true ;		 
	    if (fabs(par->get_tbl_mod(0)(0) - mp_aff->get_alpha()[0]) > 2.e-15)
		param_new = true ; 
	    for (int l=1; l<nz-1; l++) {
		if (par->get_itbl_mod(0)(l) != mgrid.get_nr(l)) param_new = true ;
		if (fabs(par->get_tbl_mod(0)(l) - mp_aff->get_beta()[l] / 
		    mp_aff->get_alpha()[l]) > 2.e-15) param_new = true ;
	    }
	    if (par->get_itbl_mod(0)(nz-1) != mgrid.get_nr(nz-1)) param_new = true ;
	    if (fabs(par->get_tbl_mod(0)(nz-1) - mp_aff->get_alpha()[nz-1]) > 2.e-15)
		param_new = true ; 
	}
	else{
	    param_new = true ;
	}
	if (param_new) {
	    par->clean_all() ;
	    par->add_int_mod(*(new int(nz)), 0) ;
	    par->add_int_mod(*(new int(lmax)), 1) ;
	    par->add_int_mod(*(new int(mgrid.get_type_t())), 2) ;
	    par->add_int_mod(*(new int(mgrid.get_type_p())), 3) ;
	    Itbl* pnr = new Itbl(nz) ;
	    pnr->set_etat_qcq() ;
	    par->add_itbl_mod(*pnr) ;
	    for (int l=0; l<nz; l++)
		pnr->set(l) = mgrid.get_nr(l) ;
	    Tbl* palpha = new Tbl(nz) ;
	    palpha->set_etat_qcq() ;
	    par->add_tbl_mod(*palpha) ;
	    palpha->set(0) = mp_aff->get_alpha()[0] ;
	    for (int l=1; l<nz-1; l++)
		palpha->set(l) = mp_aff->get_beta()[l] / mp_aff->get_alpha()[l] ;
	    palpha->set(nz-1) = mp_aff->get_alpha()[nz-1] ;
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

    sol_Dirac_l01(hh, hrr, tilde_eta, par) ;
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
    if (need_calculation && (par != 0x0)) mat_done.annule_hard() ;
		
    for (int k=0 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    // quantic numbers and spectral bases
	    donne_lm(nz, lz, j, k, base, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
	    {
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
		    if ((par != 0x0) && (mat_done(l_q) == 0)) {
			Matrice* pope = new Matrice(ope) ;
			par->add_matrice_mod(*pope, lz*lmax + l_q) ;
			mat_done.set(l_q) = 1 ;
		    }
		} //End of case when a calculation is needed

		const Matrice& oper = (par == 0x0 ? ope : 
				       par->get_matrice_mod(lz*lmax + l_q) ) ;
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

    for (int lz=1; lz<nz-1; lz++) {
	if (need_calculation && (par != 0x0)) mat_done.annule_hard() ;
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
		donne_lm(nz, lz, j, k, base, m_q, l_q, base_r) ;
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
		{
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
		    if ((par != 0x0) && (mat_done(l_q) == 0)) {
			Matrice* pope = new Matrice(ope) ;
			par->add_matrice_mod(*pope, lz*lmax + l_q) ;
			mat_done.set(l_q) = 1 ;
		    }
		    } //End of case when a calculation is needed
		    const Matrice& oper = (par == 0x0 ? ope : 
				       par->get_matrice_mod(lz*lmax + l_q) ) ;
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
    {int lz = nz-1 ;  
    if (need_calculation && (par != 0x0)) mat_done.annule_hard() ;
    int nr = mgrid.get_nr(lz) ;
    int ind0 = 0 ;
    int ind1 = nr ;
    int ind2 = 2*nr ;
    double alpha = mp_aff->get_alpha()[lz] ;
    Matrice ope(3*nr, 3*nr) ;
	
    for (int k=0 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    // quantic numbers and spectral bases
	    donne_lm(nz, lz, j, k, base, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
	    {
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
		if ((par != 0x0) && (mat_done(l_q) == 0)) {
		    Matrice* pope = new Matrice(ope) ;
		    par->add_matrice_mod(*pope, lz*lmax + l_q) ;
		    mat_done.set(l_q) = 1 ;
		}
		} //End of case when a calculation is needed
		const Matrice& oper = (par == 0x0 ? ope : 
				       par->get_matrice_mod(lz*lmax + l_q) ) ;

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

    int taille = 3*(nz-1) ;
    Mtbl_cf& mhrr = *hrr.set_spectral_va().c_cf ;
    Mtbl_cf& meta = *tilde_eta.set_spectral_va().c_cf ;
    Mtbl_cf& mw = *ww.set_spectral_va().c_cf ;
	
    Tbl sec_membre(taille) ; 
    Matrice systeme(taille, taille) ; 
    int ligne ;  int colonne ;
	
    // Loop on l and m
    //----------------
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	    if ((nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1)){
		ligne = 0 ;
		colonne = 0 ;
		systeme.annule_hard() ;
		sec_membre.annule_hard() ;

		//Nucleus 
		int nr = mgrid.get_nr(0) ;
		
		systeme.set(ligne, colonne) = sol_hom3_hrr.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_hrr.val_out_bound_jk(0, j, k) ;
		ligne++ ;

		systeme.set(ligne, colonne) = sol_hom3_eta.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_eta.val_out_bound_jk(0, j, k) ;
		ligne++ ;

		systeme.set(ligne, colonne) = sol_hom3_w.val_out_bound_jk(0, j, k) ;
		sec_membre.set(ligne) = -sol_part_w.val_out_bound_jk(0, j, k) ;
		colonne++ ;

		//shells
		for (int zone=1 ; zone<nz-1 ; zone++) {
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
    
		//Compactified external domain
		nr = mgrid.get_nr(nz-1) ;

		ligne -= 2 ;

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
		ligne++ ;
		
		systeme.set(ligne, colonne) = 
		    - sol_hom1_w.val_in_bound_jk(nz-1, j, k) ;
		systeme.set(ligne, colonne+1) = 
		    - sol_hom2_w.val_in_bound_jk(nz-1, j, k) ;
		
		sec_membre.set(ligne) += sol_part_w.val_in_bound_jk(nz-1, j, k) ;
			
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
 		for (int zone=1 ; zone<nz-1 ; zone++) { //shells
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
 		nr = mgrid.get_nr(nz-1) ; //compactified external domain
 		for (int i=0 ; i<nr ; i++) {
		    mhrr.set(nz-1, k, j, i) = sol_part_hrr(nz-1, k, j, i)
			+ facteur(conte)*sol_hom1_hrr(nz-1, k, j, i) 
			+ facteur(conte+1)*sol_hom2_hrr(nz-1, k, j, i) ;

		    meta.set(nz-1, k, j, i) = sol_part_eta(nz-1, k, j, i)
			+ facteur(conte)*sol_hom1_eta(nz-1, k, j, i) 
			+ facteur(conte+1)*sol_hom2_eta(nz-1, k, j, i) ; 
			
		    mw.set(nz-1, k, j, i) = sol_part_w(nz-1, k, j, i)
			+ facteur(conte)*sol_hom1_w(nz-1, k, j, i) 
			+ facteur(conte+1)*sol_hom2_w(nz-1, k, j, i) ;
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

namespace {
    Map* mp_ref = 0x0 ;
    int l_max_ref = -1 ;
    Matrice** t_mat ;
}

void Sym_tensor_trans::solve_hrr(const Scalar& sou_hrr, Scalar& hrr_new, int l_in_min,
				 int l_in_max) const {

    assert(*mp == sou_hrr.get_mp() ) ;
    assert(*mp == hrr_new.get_mp() ) ;

    const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
    assert(mp_aff != 0x0) ;

    assert (sou_hrr.get_etat() != ETATNONDEF) ;
    assert (sou_hrr.check_dzpuis(0)) ;
    if (sou_hrr.get_etat() == ETATZERO) {
	hrr_new.set_etat_zero() ;
	return ;
    }

    Scalar source = sou_hrr ;
    source.set_spectral_va().ylm() ;
    assert (source.get_spectral_va().c_cf != 0x0) ;
    assert (source.get_spectral_va().c_cf->get_etat() != ETATNONDEF) ;
    if (source.get_spectral_va().c_cf->get_etat() == ETATZERO) {
	hrr_new.set_etat_zero() ;
	return ;
    }
    const Base_val& base = source.get_spectral_base() ;
    int l_max = base.give_lmax(*mp_aff->get_mg(), 0) ;
    if (l_in_max < 0) l_in_max = l_max ;

    hrr_new.annule_hard() ;
    hrr_new.set_spectral_base(base) ;
    hrr_new.set_spectral_va().ylm() ;

    int nz = mp_aff->get_mg()->get_nzone() ;
    int nzm1 = nz - 1 ;
    assert (mp_aff->get_mg()->get_type_r(0) == RARE) ;
    assert (mp_aff->get_mg()->get_type_r(nzm1) == UNSURR) ;
    int np = mp_aff->get_mg()->get_np(0) ;
    int nt = mp_aff->get_mg()->get_nt(0) ;
    int taille_ope = 0 ;
    for (int lz=0; lz<nz; lz++) {
	taille_ope += mp_aff->get_mg()->get_nr(lz) ;
	assert ( mp_aff->get_mg()->get_np(lz) == np ) ;
	assert ( mp_aff->get_mg()->get_nt(lz) == nt ) ;
    }

    bool need_calculation = false ;
    if (mp_ref == 0x0) {
	need_calculation = true ;
	mp_ref = new Map_af(*mp_aff) ;
	t_mat = new Matrice*[l_max+1] ; //+1 to be safe...
	for (int ll=0; ll<l_max+1; ll++)
	    t_mat[ll] = new Matrice(taille_ope, taille_ope) ;
	l_max_ref = l_max ;
    }
    else {
	if (!(*mp_ref == *mp_aff)) {
	    need_calculation = true ;
	    delete mp_ref ;
	    mp_ref = new Map_af(*mp_aff) ;
	    for (int ll=0; ll<l_max_ref+1; ll++)
		delete t_mat[ll] ;
	    if (l_max != l_max_ref) {
		delete [] t_mat ;
		t_mat = new Matrice*[l_max+1] ;
		l_max_ref = l_max ;
	    }
	    for (int ll=0; ll<l_max+1; ll++)
		t_mat[ll] = new Matrice(taille_ope, taille_ope) ;
	}
    }
    assert (l_max == l_max_ref) ;
    assert (t_mat != 0x0) ;
    int l_q, m_q, base_r ;

    //-----------------
    // Loop on l and m
    //-----------------
    for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	    if ((nullite_plm(j, nt, k, np, base) == 1) 
		&& (l_q>=l_in_min) && (l_q<=l_in_max)){
	    assert (t_mat[l_q] != 0x0) ;
	    Matrice& ope = *(t_mat[l_q]) ;
	    if (need_calculation) 
		ope.annule_hard() ;

	    Tbl ty(taille_ope) ;
	    ty.set_etat_qcq() ;
	    int parite = ( ( (l_q % 2) == 0) ? R_CHEBP : R_CHEBI );
	    bool big_l = ( l_q > 3) ;
	    int lin_ref = 0 ;
	    int col_ref = 0 ;
    {
	//Nucleus
	int nr = mp_aff->get_mg()->get_nr(0) ;
	assert (parite == base_r) ;
	if (need_calculation) {
	    Diff_xdsdx dx(parite, nr) ; const Matrice& mdx = dx.get_matrice() ;
	    Diff_x2dsdx2 dx2(parite, nr) ; const Matrice& md2 = dx2.get_matrice() ;
	    Diff_id sx2(parite, nr) ; const Matrice& ms2 = sx2.get_matrice() ;
	    int nr_nuc = ( big_l ? nr - 1 : nr ) ;
	
	    for (int lin=0; lin<nr_nuc; lin++) {
		for (int col=0; col<nr; col++)
		    ope.set(lin,col) = md2(lin, col)  + 7*mdx(lin, col) 
			+ (9-double(l_q)*double(l_q+1)/double(2))*ms2(lin, col) ;
		ty.set(lin) = (*source.get_spectral_va().c_cf)(0, k, j, lin) ;
	    }
	for (int col=0; col<nr; col++)
	    ope.set(nr_nuc, col) = 1 ;
	double alpha = mp_aff->get_alpha()[0] ;
	if  ( parite == R_CHEBP ) 
	    for (int col = 0; col<nr; col++) 
		ope.set(nr_nuc+1, col) = 4*col*col / alpha ;
	else
	    for (int col=0; col<nr; col++)
		ope.set(nr_nuc+1, col) = (2*col+1)*(2*col+1)/alpha ;
	lin_ref += nr_nuc ;
	col_ref += nr ;
	}
	else {
	    int nr_nuc = ( big_l ? nr - 1 : nr ) ;
	    for (int lin=0; lin<nr_nuc; lin++) {
		ty.set(lin) = (*source.get_spectral_va().c_cf)(0, k, j, lin) ;
	    }
	lin_ref += nr_nuc ;
	}
    }

    //Shells
    for (int lz=1; lz<nzm1; lz++) {
	int nr = mp_aff->get_mg()->get_nr(lz) ;
	if (need_calculation) {
	    double alpha = mp_aff->get_alpha()[lz] ;
	    double ech = mp_aff->get_beta()[lz] / mp_aff->get_alpha()[lz] ;
	    int sign = -1 ;
	    for (int col=0; col<nr; col++) {
		ope.set(lin_ref, col_ref+col) = sign ;
		sign *= -1 ;
	    }
	    ty.set(lin_ref) = 0 ;
	    sign = 1 ;
	    for (int col=0; col<nr; col++) {
		ope.set(lin_ref+1, col_ref+col) = sign*col*col/alpha  ;
		sign *= -1 ;
	    }
	    ty.set(lin_ref+1) = 0 ;
	    lin_ref += 2 ;

	    Diff_x2dsdx2 x22(R_CHEB, nr) ; const Matrice& m22 = x22.get_matrice() ;
	    Diff_xdsdx2 x12(R_CHEB, nr) ; const Matrice& m12 = x12.get_matrice() ;
	    Diff_dsdx2 x02(R_CHEB, nr) ; const Matrice& m02 = x02.get_matrice() ;
	    Diff_xdsdx x11(R_CHEB, nr) ; const Matrice& m11 = x11.get_matrice() ;
	    Diff_dsdx x01(R_CHEB, nr) ; const Matrice& m01 = x01.get_matrice() ;
	    Diff_id xid(R_CHEB, nr) ; const Matrice& mid = xid.get_matrice() ;
	    
	    for (int lin=0; lin<nr-2; lin++) {
		for (int col=0; col<nr; col++) { 
		    ope.set(lin+lin_ref,col+col_ref) = 
			m22(lin,col) + 2*ech*m12(lin,col) + ech*ech*m02(lin,col) 
			+ 7*(m11(lin,col) + ech*m01(lin,col))
			+ (9 - 0.5*l_q*(l_q+1))*mid(lin,col) ;
		}
		ty.set(lin+lin_ref) 
		    = (*source.get_spectral_va().c_cf)(lz, k, j, lin) ;
	    }
	    
	    lin_ref += nr - 2 ;
	    for (int col=col_ref; col<col_ref+nr; col++)
		ope.set(lin_ref,col) = 1 ; 
	    for (int col=0; col<nr; col++)
		ope.set(lin_ref+1,col+col_ref) = col*col/alpha ;
	    col_ref += nr ;
	}
	else {
	    ty.set(lin_ref) = 0 ;
	    ty.set(lin_ref+1) = 0 ;
	    lin_ref += 2 ;
	    for (int lin=0; lin<nr-2; lin++) 
		ty.set(lin+lin_ref) 
		    = (*source.get_spectral_va().c_cf)(lz, k, j, lin) ;
	    lin_ref += nr - 2 ;
	}	    
    }
    
    {    // Compactified external domain  
	int nr =  mp_aff->get_mg()->get_nr(nzm1) ;
	int n_hom = (big_l ? 1 : 2) ;
	if (need_calculation) {
	    double alpha = mp_aff->get_alpha()[nzm1] ;
	    int sign = -1 ;
	    for (int col=0; col<nr; col++) {
		ope.set(lin_ref, col_ref+col) = sign ;
		sign *= -1 ;
	    }
	    ty.set(lin_ref) = 0 ;
	    sign = -4 ;
	    for (int col=0; col<nr; col++) {
		ope.set(lin_ref+1, col_ref+col) = sign*col*col*alpha  ;
		sign *= -1 ;
	    }
	    ty.set(lin_ref+1) = 0 ;
	    lin_ref += 2 ;
	    
	    Diff_x2dsdx2 x22(R_CHEBU, nr) ; const Matrice& m22 = x22.get_matrice();
	    Diff_xdsdx x11(R_CHEBU, nr) ; const Matrice& m11 = x11.get_matrice() ;
	    Diff_id xid(R_CHEBU, nr) ; const Matrice& mid = xid.get_matrice() ;
	    
	    for (int lin=0; lin<nr-n_hom; lin++) {
		for (int col=0; col<nr; col++) 
		    ope.set(lin+lin_ref, col+col_ref) = 
			m22(lin,col) - 5*m11(lin,col) 
			+ (9 - 0.5*l_q*(l_q+1))*mid(lin,col) ;
		ty.set(lin+lin_ref) 
		    = (*source.get_spectral_va().c_cf)(nzm1, k, j, lin) ;
	    }
	}
	else {
	    ty.set(lin_ref) = 0 ;
	    ty.set(lin_ref+1) = 0 ;
	    lin_ref += 2 ;
	    for (int lin=0; lin<nr-n_hom; lin++) {
		ty.set(lin+lin_ref) 
		    = (*source.get_spectral_va().c_cf)(nzm1, k, j, lin) ;
	    }
	}
    }
    ope.set_lu() ;
    Tbl tx = ope.inverse(ty) ;

    int compte = 0 ;
    for (int lz=0; lz<nz; lz++) {
	int nr = mp_aff->get_mg()->get_nr(lz) ;
	for (int i=0; i<nr; i++) {
	    hrr_new.set_spectral_va().c_cf->set(lz, k, j, i)
		= tx(compte) ;
	    compte++ ;
	}
    }
	    } // End of nullite_plm (=> l,m loop)
	}// Theta loop
    hrr_new.set_spectral_va().ylm_i() ;
    if (hrr_new.set_spectral_va().c != 0x0) 
	delete hrr_new.set_spectral_va().c ;
    hrr_new.set_spectral_va().c = 0x0 ;
    
    return ;
}

void Sym_tensor_trans::sol_Dirac_l01(const Scalar& hh, Scalar& hrr, Scalar& tilde_eta,
				     Param* par) const {

    const Map_af* mp_aff = dynamic_cast<const Map_af*>(mp) ;
    assert(mp_aff != 0x0) ; //Only affine mapping for the moment

    const Mg3d& mgrid = *mp_aff->get_mg() ;
    int nz = mgrid.get_nzone() ;
    assert(mgrid.get_type_r(0) == RARE)  ;
    assert(mgrid.get_type_r(nz-1) == UNSURR) ;

    if (hh.get_etat() == ETATZERO) {
	hrr.annule_hard() ;
	tilde_eta.annule_hard() ;
	return ;
    }

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
    if (par != 0x0)
	if (par->get_n_matrice_mod() > 0) 
	    if (&par->get_matrice_mod(0) != 0x0) need_calculation = false ;

    int l_q, m_q, base_r ;
    Itbl mat_done(lmax) ;

    //---------------
    //--  NUCLEUS ---
    //---------------
    {int lz = 0 ;  
    int nr = mgrid.get_nr(lz) ;
    double alpha = mp_aff->get_alpha()[lz] ;
    Matrice ope(2*nr, 2*nr) ;
    if (need_calculation && (par != 0x0)) mat_done.annule_hard() ;

    for (int k=0 ; k<np+1 ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    // quantic numbers and spectral bases
	    donne_lm(nz, lz, j, k, base, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q < 2))
	    {
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
		    if ((par != 0x0) && (mat_done(l_q) == 0)) {
			Matrice* pope = new Matrice(ope) ;
			par->add_matrice_mod(*pope, lz*lmax + l_q) ;
			mat_done.set(l_q) = 1 ;
		    }
		} //End of case when a calculation is needed

		const Matrice& oper = (par == 0x0 ? ope : 
				       par->get_matrice_mod(lz*lmax + l_q) ) ;
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
	if (need_calculation && (par != 0x0)) mat_done.annule_hard() ;
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
		donne_lm(nz, lz, j, k, base, m_q, l_q, base_r) ;
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q < 2))
		{
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
		    if ((par != 0x0) && (mat_done(l_q) == 0)) {
			Matrice* pope = new Matrice(ope) ;
			par->add_matrice_mod(*pope, lz*lmax + l_q) ;
			mat_done.set(l_q) = 1 ;
		    }
		    } //End of case when a calculation is needed
		    const Matrice& oper = (par == 0x0 ? ope : 
				       par->get_matrice_mod(lz*lmax + l_q) ) ;
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
    if (need_calculation && (par != 0x0)) mat_done.annule_hard() ;
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
	    donne_lm(nz, lz, j, k, base, m_q, l_q, base_r) ;
	    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q < 2))
	    {
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
		if ((par != 0x0) && (mat_done(l_q) == 0)) {
		    Matrice* pope = new Matrice(ope) ;
		    par->add_matrice_mod(*pope, lz*lmax + l_q) ;
		    mat_done.set(l_q) = 1 ;
		}
		} //End of case when a calculation is needed
		const Matrice& oper = (par == 0x0 ? ope : 
				       par->get_matrice_mod(lz*lmax + l_q) ) ;
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
	    if ((nullite_plm(j, nt, k, np, base) == 1) && (l_q < 2)){
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

