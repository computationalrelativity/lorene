/*
 *  Resolution of the divergence ODE: df/df + n*f/r = source (source must have dzpuis =2)
 *
 *    (see file scalar.h for documentation).
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

char scalar_sol_div_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:53:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:16:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2005/09/16 14:33:00  j_novak
 * Added #include <math.h>.
 *
 * Revision 1.2  2005/09/16 12:49:52  j_novak
 * The case with dzpuis=1 is added.
 *
 * Revision 1.1  2005/06/08 12:35:22  j_novak
 * New method for solving divergence-like ODEs.
 *
 *
 * $Header$
 *
 */

// C headers
#include <cassert>
#include <cmath>

//Lorene headers
#include "tensor.h"
#include "diff.h"
#include "proto.h"

// Local prototypes 
namespace Lorene {
void _sx_r_chebp(Tbl* , int& ) ;
void _sx_r_chebi(Tbl* , int& ) ;


Scalar Scalar::sol_divergence(int n_factor) const {

    assert(etat != ETATNONDEF) ;
    const Map_af* mpaff = dynamic_cast<const Map_af*>(mp) ;
    assert( mpaff != 0x0) ;

    Scalar result(*mp) ;

    if ( etat == ETATZERO )
	result.set_etat_zero() ;
    else {                                             //source not zero
	Base_val base_resu = get_spectral_base() ;
	base_resu.mult_x() ;
	const Mg3d* mg = mp->get_mg() ;
	result.set_etat_qcq() ; result.set_spectral_base(base_resu) ;
	result.set_spectral_va().set_etat_cf_qcq() ;
	Valeur sigma(va) ;
 	sigma.ylm_i() ;                               // work on Fourier basis
 	const Mtbl_cf& source = *sigma.c_cf ;

	// Checks on the type of domains
	int nz = mg->get_nzone() ;
	bool ced = (mg->get_type_r(nz-1) == UNSURR ) ;
	assert ( (!ced) || (check_dzpuis(2)) || (check_dzpuis(1)) ) ;
	assert (mg->get_type_r(0) == RARE) ;
	int nt = mg->get_nt(0) ;
	int np = mg->get_np(0) ;
#ifndef NDEBUG
	for (int lz = 0; lz<nz; lz++)
	    assert( (mg->get_nt(lz) == nt) && (mg->get_np(lz) == np) ) ;
#endif
	int nr, base_r,l_quant, m_quant;
	Tbl *so ;
	Tbl *s_hom ;
	Tbl *s_part ;
  
	// Working objects and initialization
	Mtbl_cf sol_part(mg, base_resu) ;
	Mtbl_cf sol_hom(mg, base_resu) ;
	Mtbl_cf& resu = *result.set_spectral_va().c_cf ;
  	sol_part.annule_hard();
	sol_hom.annule_hard() ;
	resu.annule_hard() ;

	//---------------
	//--  NUCLEUS ---
	//---------------
	int lz = 0 ;  
	nr = mg->get_nr(lz) ;

        int dege = 1 ; // the operator is degenerate
	int nr0 = nr - dege ;
	Tbl vect1(3, 1, nr) ; 
	Tbl vect2(3, 1, nr) ;
	int base_pipo = 0 ;
	double alpha = mpaff->get_alpha()[lz] ;
	double beta =  0. ;
	Matrice ope_even(nr0, nr0) ; //when the *result* is decomposed on R_CHEBP
	ope_even.set_etat_qcq() ;
	for (int i=dege; i<nr; i++) {
	    vect1.annule_hard() ;
	    vect2.annule_hard() ;
	    vect1.set(0,0,i) = 1. ; vect2.set(0,0,i) = 1. ;
	    _dsdx_r_chebp(&vect1, base_pipo) ;
	    _sx_r_chebp(&vect2, base_pipo) ;
	    for (int j=0; j<nr0; j++)
		ope_even.set(j,i-dege) = (vect1(0,0,j) + n_factor*vect2(0,0,j)) / alpha ;
	}
	ope_even.set_lu() ;
	Matrice ope_odd(nr0, nr0) ; //when the *result* is decomposed on R_CHEBI
	ope_odd.set_etat_qcq() ;
	for (int i=0; i<nr0; i++) {
	    vect1.annule_hard() ;
	    vect2.annule_hard() ;
	    vect1.set(0,0,i) = 1. ; vect2.set(0,0,i) = 1. ;
	    _dsdx_r_chebi(&vect1, base_pipo) ;
	    _sx_r_chebi(&vect2, base_pipo) ;
	    for (int j=0; j<nr0; j++)
		ope_odd.set(j,i) = (vect1(0,0,j) + n_factor*vect2(0,0,j)) / alpha ;
	}
	ope_odd.set_lu() ;
  
	for (int k=0 ; k<np+1 ; k++) 
	    for (int j=0 ; j<nt ; j++) {
		// to get the spectral base
		base_resu.give_quant_numbers(lz, k, j, m_quant, l_quant, base_r) ;
		assert ( (base_r == R_CHEBP) || (base_r == R_CHEBI) ) ;
		const Matrice& operateur = (( base_r == R_CHEBP ) ? 
					    ope_even : ope_odd ) ;
		// particular solution 
		so = new Tbl(nr0) ;
		so->set_etat_qcq() ;
		for (int i=0 ; i<nr0 ; i++)
		    so->set(i) = source(lz, k, j, i) ;

		s_part = new Tbl(operateur.inverse(*so)) ;
		
		// Putting to Mtbl_cf	  
		double somme = 0 ;
		for (int i=0 ; i<nr0 ; i++) {
		    if (base_r == R_CHEBP) {
			resu.set(lz, k, j, i+dege) = (*s_part)(i) ;
			somme += ((i+dege)%2 == 0 ? 1 : -1)*(*s_part)(i) ;
		    }
		    else 
			resu.set(lz,k,j,i) = (*s_part)(i) ;
		}
		if (base_r == R_CHEBI) 
		    for (int i=nr0; i<nr; i++) 
			resu.set(lz,k,j,i) = 0 ; 
		if (base_r == R_CHEBP) //result must vanish at r=0
		    resu.set(lz, k, j, 0) -= somme ;

		delete so ;
		delete s_part ;
		
	    }
	
	//---------------------
	//--      SHELLS     --
	//---------------------
	int nz0 = (ced ? nz - 1 : nz) ;
	for (lz=1 ; lz<nz0 ; lz++) {
	    nr = mg->get_nr(lz) ;    
	    alpha = mpaff->get_alpha()[lz] ;
	    beta =  mpaff->get_beta()[lz];
	    double ech = beta / alpha ;
	    Diff_id id(R_CHEB, nr) ; const Matrice& mid = id.get_matrice() ; 
	    Diff_xdsdx xd(R_CHEB, nr) ; const Matrice& mxd = xd.get_matrice() ;
	    Diff_dsdx dx(R_CHEB, nr) ; const Matrice& mdx = dx.get_matrice() ;
	    Matrice operateur = mxd + ech*mdx + n_factor*mid ;
	    operateur.set_lu() ;
	    // homogeneous solution
	    s_hom = new Tbl(solh(nr, n_factor-1, ech, R_CHEB)) ;
	    
	    for (int k=0 ; k<np+1 ; k++)
		for (int j=0 ; j<nt ; j++) {
		    // to get the spectral base
		    base_resu.give_quant_numbers(lz, k, j, m_quant, l_quant, base_r) ;
		    assert (base_r == R_CHEB) ;
		    
		    so = new Tbl(nr) ;
		    so->set_etat_qcq() ;
		    // particular solution
		    Tbl tmp(nr) ;
		    tmp.set_etat_qcq() ;
		    for (int i=0 ; i<nr ; i++)
			tmp.set(i) = source(lz, k, j, i) ;
		    for (int i=0; i<nr; i++) so->set(i) = beta*tmp(i) ;
		    multx_1d(nr, &tmp.t, R_CHEB) ;
		    for (int i=0; i<nr; i++) so->set(i) += alpha*tmp(i)  ;
		    
		    s_part = new Tbl (operateur.inverse(*so)) ;
		    
		    // cleaning things...
		    for (int i=0 ; i<nr ; i++) {
			sol_part.set(lz, k, j, i) = (*s_part)(i) ;
			sol_hom.set(lz, k, j, i) = (*s_hom)(1,i) ;
		    }
		    
		    delete so ;
		    delete s_part ;
		}
	    delete s_hom ;
	}
	if (ced) {
	//---------------
	//--  CED   -----
	//---------------
	    int dzp = ( check_dzpuis(2) ? 2 : 1) ;
	    nr = source.get_mg()->get_nr(nz-1) ;
	    alpha = mpaff->get_alpha()[nz-1] ;
	    beta = mpaff->get_beta()[nz-1] ;
	    dege = dzp  ;
	    nr0 = nr - dege ;
	    Diff_dsdx dx(R_CHEBU, nr) ; const Matrice& mdx = dx.get_matrice() ;
	    Diff_sx sx(R_CHEBU, nr) ; const Matrice& msx = sx.get_matrice() ;
	    Diff_xdsdx xdx(R_CHEBU, nr) ; const Matrice& mxdx = xdx.get_matrice() ;
	    Diff_id id(R_CHEBU, nr) ; const Matrice& mid = id.get_matrice() ;
	    Matrice operateur(nr0, nr0) ;
	    operateur.set_etat_qcq() ;
	    if (dzp == 2)
		for (int lin=0; lin<nr0; lin++) 
		    for (int col=dege; col<nr; col++)
			operateur.set(lin,col-dege) = (-mdx(lin,col) 
					+ n_factor*msx(lin, col)) / alpha ;
	    else {
 		for (int lin=0; lin<nr0; lin++) {
 		    for (int col=dege; col<nr; col++)
 			operateur.set(lin,col-dege) = (-mxdx(lin,col) 
 					+ n_factor*mid(lin, col)) ;
		}
	    }
	    operateur.set_lu() ;
	    // homogeneous solution
	    s_hom = new Tbl(solh(nr, n_factor-1, 0., R_CHEBU)) ;
	    for (int k=0 ; k<np+1 ; k++)
		for (int j=0 ; j<nt ; j++) {
		    base_resu.give_quant_numbers(lz, k, j, m_quant, l_quant, base_r) ;
		    assert(base_r == R_CHEBU) ;	    
		    
		    // particular solution
		    so = new Tbl(nr0) ;
		    so->set_etat_qcq() ;
		    for (int i=0 ; i<nr0 ; i++)
			so->set(i) = source(nz-1, k, j, i) ;
		    s_part = new Tbl(operateur.inverse(*so)) ;

		    // cleaning
		    double somme = 0 ;
		    for (int i=0 ; i<nr0 ; i++) {
			sol_part.set(nz-1, k, j, i+dege) = (*s_part)(i) ;
			somme += (*s_part)(i) ;
			sol_hom.set(nz-1, k, j, i) = (*s_hom)(i) ;
		    }
		    for (int i=nr0; i<nr; i++)
			sol_hom.set(nz-1, k, j, i) = (*s_hom)(i) ;
		    //result must vanish at infinity
		    sol_part.set(nz-1, k, j, 0) = -somme ; 
		    delete so ;
		    delete s_part ;
		}
	    delete s_hom ;
	}
	
	//-------------------------
	//-- matching solutions ---
	//-------------------------
	if (nz > 1) { 
	    Tbl echelles(nz-1) ;
	    echelles.set_etat_qcq() ;
	    for (lz=1; lz<nz; lz++) 
		echelles.set(lz-1) 
		    = pow ( (mpaff->get_beta()[lz]/mpaff->get_alpha()[lz] -1), 
			    n_factor) ;
	    if (ced) echelles.set(nz-2) = 1./pow(-2., n_factor) ;
	    
	    for (int k=0 ; k<np+1 ; k++)
		for (int j=0 ; j<nt ; j++) {
		    for (lz=1; lz<nz; lz++) {
			double val1 = 0 ;
			double valm1 = 0 ;
			double valhom1 = 0 ;
			int nr_prec = mg->get_nr(lz-1) ;
			nr = mg->get_nr(lz) ;
			for (int i=0; i<nr_prec; i++)
			    val1 += resu(lz-1, k, j, i) ;
			for (int i=0; i<nr; i++) {
			    valm1 += ( i%2 == 0 ? 1 : -1)*sol_part(lz, k, j, i) ;
			    valhom1 += ( i%2 == 0 ? 1 : -1)*sol_hom(lz, k, j, i) ;
			}
			double lambda = (val1 - valm1) * echelles(lz-1) ;
			for (int i=0; i<nr; i++)
			    resu.set(lz, k, j, i) = sol_part(lz, k, j, i) 
				+ lambda*sol_hom(lz, k, j, i) ;

		    }
		}
	}
    }
    return result ;
}

}
