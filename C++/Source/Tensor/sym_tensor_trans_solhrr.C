/*
 *  Solution of the 2nd-order PDE for hrr (2 transverse conditions)
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

char sym_tensor_trans_solhrr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2005/09/07 16:47:43  j_novak
 * Removed method Sym_tensor_trans::T_from_det_one
 * Modified Sym_tensor::set_auxiliary, so that it takes eta/r and mu/r as
 * arguments.
 * Modified Sym_tensor_trans::set_hrr_mu.
 * Added new protected method Sym_tensor_trans::solve_hrr
 *
 *
 * $Header$
 *
 */

// C headers
#include <assert.h>

// Lorene headers
#include "tensor.h"
#include "diff.h"
#include "proto.h"
#include "graphique.h"

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

    Matrice ope(taille_ope, taille_ope) ;
    ope.set_etat_qcq() ;
    int l_q, m_q, base_r ;

    for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) 
	    if (nullite_plm(j, nt, k, np, base) == 1) {
	    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
	    for (int lin=0; lin<taille_ope; lin++)
		for (int col=0; col<taille_ope; col++)
		    ope.set(lin,col) = 0 ;

	    Tbl ty(taille_ope) ;
	    ty.set_etat_qcq() ;
	    int parite = ( ( (l_q % 2) == 0) ? R_CHEBP : R_CHEBI );
	    bool big_l = ( l_q > 3) ;
	    int lin_ref = 0 ;
	    int col_ref = 0 ;
    {
	//Nucleus
	int nr = mp_aff->get_mg()->get_nr(0) ;
	Diff_xdsdx dx(parite, nr) ; const Matrice& mdx = dx.get_matrice() ;
	Diff_x2dsdx2 dx2(parite, nr) ; const Matrice& md2 = dx2.get_matrice() ;
	Diff_id sx2(parite, nr) ; const Matrice ms2 = sx2.get_matrice() ;
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

    //Shells
    for (int lz=1; lz<nzm1; lz++) {
	int nr = mp_aff->get_mg()->get_nr(lz) ;
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
	Diff_id xid(R_CHEB, nr) ; const Matrice mid = xid.get_matrice() ;

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

    {    // Compactified external domain  
	int nr =  mp_aff->get_mg()->get_nr(nzm1) ;
	int n_hom = (big_l ? 1 : 2) ;
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
	Diff_id xid(R_CHEBU, nr) ; const Matrice mid = xid.get_matrice() ;
   
	for (int lin=0; lin<nr-n_hom; lin++) {
	    for (int col=0; col<nr; col++) 
		ope.set(lin+lin_ref, col+col_ref) = 
		    m22(lin,col) - 5*m11(lin,col) 
		    + (9 - 0.5*l_q*(l_q+1))*mid(lin,col) ;
	    ty.set(lin+lin_ref) 
		= (*source.get_spectral_va().c_cf)(nzm1, k, j, lin) ;
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
