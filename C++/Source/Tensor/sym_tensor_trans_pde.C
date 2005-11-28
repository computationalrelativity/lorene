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

// Lorene headers
#include "tensor.h"
#include "diff.h"
#include "proto.h"

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

namespace {
    Map* mp_ref = 0x0 ;
    int l_max_ref = -1 ;
    Matrice** t_mat ;
}

void Sym_tensor_trans::solve_hrr(const Scalar& sou_hrr, Scalar& hrr_new) const {

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

    for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) 
	    if (nullite_plm(j, nt, k, np, base) == 1) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
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
//    cout << "l_q: " << l_q << ", determinant: " << ope.determinant() << endl ;
    ope.set_lu() ;
    
//    Matrice myy(ty) ;

    Tbl tx = ope.inverse(ty) ;

//      Matrice mxx(tx) ;
//      cout << "difference: " << 
//  	diffrelmax((ope*mxx).get_array(), myy.get_array()) << endl ;
    
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
    hrr_new.set_spectral_va().ylm_i() ;
    if (hrr_new.set_spectral_va().c != 0x0) 
	delete hrr_new.set_spectral_va().c ;
    hrr_new.set_spectral_va().c = 0x0 ;
//     double Rdes = 2*mp_aff->val_r_jk(nz-2, 1., 0, 0) ;
//     des_profile(hrr_new, 0., Rdes, 1, 1) ;

//     Scalar verif = 0.5*hrr_new.lapang() ;
//     verif.set_spectral_va().ylm_i() ;
//     verif += 9*hrr_new ;
//     Scalar tmp = hrr_new.dsdr().dsdr() ;
//     tmp.mult_r_dzpuis(2) ;
//     tmp += 7*hrr_new.dsdr() ;
//     tmp.mult_r_dzpuis(0) ;
//     verif += tmp ;
//     verif.set_spectral_va().ylm() ;
//     verif -= source ;
//     verif.spectral_display("sym_tensor_trans::solve_hrr: check of accuracy") ;
    
    return ;
}
