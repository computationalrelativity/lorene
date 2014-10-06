/*
 *  Functions to solve various PDEs for a divergence-free symmetric tensor.
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2006  Jerome Novak
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
 * Revision 1.15  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.14  2010/10/11 10:38:34  j_novak
 * *** empty log message ***
 *
 * Revision 1.13  2010/10/11 10:23:03  j_novak
 * Removed methods Sym_tensor_trans::solve_hrr() and Sym_tensor_trans::set_WX_det_one(), as they are no longer relevant.
 *
 * Revision 1.12  2006/09/05 15:38:45  j_novak
 * The fuctions sol_Dirac... are in a seperate file, with new parameters to
 * control the boundary conditions.
 *
 * Revision 1.11  2006/08/31 12:13:22  j_novak
 * Added an argument of type Param to Sym_tensor_trans::sol_Dirac_A().
 *
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
	Scalar res_A(*mp) ; Scalar big_A = res_cart.compute_A() ;
	Scalar res_B(*mp) ; Scalar big_B = res_cart.compute_tilde_B_tt() ;

 	switch (big_A.get_etat() ) {
 	    case ETATZERO: {
		res_A.set_etat_zero() ;
 		break ;
 	    }
 	    case ETATUN : {
 		res_A.set_etat_one() ;
 		break ;
 	    }
 	    case ETATQCQ : {
 		res_A.allocate_all() ;
 		for (int lz=0; lz<nz; lz++) 			
 		    for (int k=0; k<np; k++)
 			for (int j=0; j<nt; j++)
 			    for(int i=0; i<gri.get_nr(lz); i++)
 				res_A.set_grid_point(lz, k, j, i)
 				    = big_A.val_grid_point(lz, k, j, i) ;
 		break ;
 	    }
 	    default: {
 		cout << 
 		    "Sym_tensor_trans::poisson() : res_A in undefined state!" 
 		     << endl ;
 		abort() ;
 		break ; 
 	    }
 	}
	res_A.set_spectral_base(big_A.get_spectral_base()) ;
	int dzA = big_A.get_dzpuis() ;
	res_A.set_dzpuis(dzA) ;

	switch (big_B.get_etat() ) {
	    case ETATZERO: {
		res_B.set_etat_zero() ;
		break ;
	    }
	    case ETATUN : {
		res_B.set_etat_one() ;
		break ;
	    }
	    case ETATQCQ : {
		res_B.allocate_all() ;
		for (int lz=0; lz<nz; lz++) 			
		    for (int k=0; k<np; k++)
			for (int j=0; j<nt; j++)
			    for(int i=0; i<gri.get_nr(lz); i++)
				res_B.set_grid_point(lz, k, j, i)
				    = big_B.val_grid_point(lz, k, j, i) ;
		break ;
	    }
	    default: {
		cout << 
		    "Sym_tensor_trans::poisson() : res_B in undefined state!" 
		     << endl ;
		abort() ;
		break ; 
	    }
	}
	res_B.set_spectral_base(big_B.get_spectral_base()) ;
	int dzB = big_B.get_dzpuis() ;
	res_B.set_dzpuis(dzB) ;
	
	resu.set_AtBtt_det_one(res_A, res_B, h_guess) ;
	
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
	
	resu.set_AtBtt_det_one(res_cart.compute_A(), res_cart.compute_tilde_B_tt(), h_guess) ;
	
    }
#ifndef NDEBUG
    Vector dive = resu.divergence(*met_div) ;
    dive.dec_dzpuis(2) ;
    maxabs(dive, "Sym_tensor_trans::poisson : divergence of the solution") ;
#endif    
    return resu ;   
}


