/*
 * Method of the class Map_radial for resolution of the equation
 *
 *  a \Delta\psi + {\bf b} \cdot \nabla \psi = \sigma
 *
 * Computes the unique solution such that psi(0) = 0.
 * bb must be given in spherical coordinates.
 *
 * (see file map.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000-2007 Eric Gourgoulhon
 *   Copyright (c) 2007 Michal Bejger
 *   Copyright (c) 2000-2001 Philippe Grandclement
 *   Copyright (c) 2001 Keisuke Taniguchi
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char map_radial_poisson_cpt_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/06 15:13:14  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2007/10/18 08:19:32  e_gourgoulhon
 * Suppression of the abort for nzet > 2 : the function should be able
 * to treat an arbitrary number of domains inside the star.
 * NB: tested only for nzet = 2.
 *
 * Revision 1.4  2007/10/16 21:53:13  e_gourgoulhon
 * Added new function poisson_compact (multi-domain version)
 *
 * Revision 1.3  2003/10/03 15:58:48  j_novak
 * Cleaning of some headers
 *
 * Revision 1.2  2002/12/11 13:17:07  k_taniguchi
 * Change the multiplication "*" to "%".
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.17  2001/08/28  15:08:06  keisuke
 *  Uncomment "sour_j.std_base_scal()" and "sour_j.ylm()".
 *
 * Revision 2.16  2001/08/28  14:45:06  keisuke
 * Uncomment "sour_j.coef()".
 *
 * Revision 2.15  2001/08/28  14:10:42  keisuke
 * Comment out "sour_j.std_base_scal()", "sour_j.coef()", and "sour_j.ylm()".
 *
 * Revision 2.14  2000/03/30  09:18:53  eric
 * Modifs affichage.
 *
 * Revision 2.13  2000/03/17  15:24:01  phil
 * suppression du nettoyage brutal
 *
 * Revision 2.12  2000/03/16  16:26:07  phil
 * Ajout du nettoyage des hautes frequences
 *
 * Revision 2.11  2000/03/10  09:18:36  eric
 * Annulation du terme l=0 de la source effective.
 *
 * Revision 2.10  2000/03/09  13:52:19  phil
 * *** empty log message ***
 *
 * Revision 2.9  2000/02/28  14:34:31  eric
 * *** empty log message ***
 *
 * Revision 2.8  2000/02/28  14:31:32  eric
 * Suppression de mpaff: remplacement de pre_lap = (1-r^2/R^2) par
 *  Id - .mult_x().mult_x().
 * Introduction de dpsi.
 * Modif affichages.
 *
 * Revision 2.7  2000/02/25  13:48:00  eric
 * Suppression des appels a Mtbl_cf::nettoie().
 *
 * Revision 2.6  2000/02/22  11:43:10  eric
 * Appel de ylm() sur d2psi.
 * Modif affichage.
 *
 * Revision 2.5  2000/02/21  12:58:12  eric
 * Modif affichage.
 *
 * Revision 2.4  2000/01/27  15:58:36  eric
 * Utilisation du nouveau constructeur Map_af::Map_af(const Map&)
 * pour la construction du mapping auxiliaire mpaff.
 * Suppression des affichages intermediaires.
 *
 * Revision 2.3  2000/01/14  17:33:44  eric
 * Amelioration du calcul (le Cmp intermediaire psi0 est supprime).
 *
 * Revision 2.2  2000/01/13  16:43:59  eric
 * Premiere version testee : OK !
 *
 * Revision 2.1  2000/01/12  16:03:23  eric
 * Premiere version complete.
 *
 * Revision 2.0  2000/01/10  09:14:52  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C++

// Headers C
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "tenseur.h"
#include "param.h"
#include "proto.h"
#include "graphique.h"
#include "utilitaires.h"

// Local prototypes
Mtbl_cf sol_poisson_compact(const Mtbl_cf&, double, double, bool) ;
Mtbl_cf sol_poisson_compact(const Map_af&, const Mtbl_cf&, const Tbl&, 
	const Tbl&, bool) ;  


	/////////////////////////////
	//	Single domain version  //
	/////////////////////////////

void Map_radial::poisson_compact(const Cmp& source, const Cmp& aa, 
				 const Tenseur& bb, const Param& par, 
				 Cmp& psi) const {


    // Protections
    // -----------
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(aa.get_etat() != ETATNONDEF) ; 
    assert(bb.get_etat() != ETATNONDEF) ; 
    assert(aa.get_mp() == source.get_mp()) ; 
    assert(bb.get_mp() == source.get_mp()) ; 
    assert(psi.get_mp() == source.get_mp()) ; 
    
    
    // The components of vector b must be given in the spherical basis
    //  associated with the mapping : 
    assert(*(bb.get_triad()) == bvect_spher) ; 

    // Computation parameters
    // ----------------------
    int nitermax = par.get_int() ; 
    int& niter = par.get_int_mod() ; 
    double precis = par.get_double(0) ; 
    double relax = par.get_double(1) ; 
    double unmrelax = 1. - relax ; 
        
    // Maybe nothing to do ?
    // ---------------------
    
    if ( source.get_etat() == ETATZERO ) {	
	psi.set_etat_zero() ; 
	return ; 
    }	
	
    // Factors alpha ( in front of (1-xi^2) Lap_xi(psi) )
    //  and    beta  ( in front of xi dpsi/dxi )
    // --------------------------------------------------
    
    Mtbl tmp = dxdr ; 
    double dxdr_c = tmp(0, 0, 0, 0) ;  // central value of 1/(dR/dxi)
    
    double alph = dxdr_c * dxdr_c * aa(0, 0, 0, 0) ; 

    const Valeur& b_r = bb(0).va ;  // radial component b^r of vector b
    
    Valeur vatmp = dxdr*dxdr*b_r ;
    
    double bet = min(vatmp)(0) ;  // Minimal value in domain no. 0
    
    cout << "Map_radial::poisson_compact : alph, bet : " << alph << "  "
	 << bet << endl ; 


    // Everything is set to zero except in the innermost domain (nucleus) :
    // ------------------------------------------------------------------

    int nz = mg->get_nzone() ;

    psi.annule(1, nz-1) ; 

    // Auxilary quantities
    // -------------------
    Cmp psi_jm1 = psi ; 
    Cmp b_grad_psi(this) ; 
    Valeur sour_j(*mg) ; 
    Valeur aux_psi(*mg) ;		
    Valeur lap_xi_psi(*mg) ; 
    Valeur oper_psi(*mg) ; 
    Valeur dpsi(*mg) ; 
    Valeur d2psi(*mg) ; 

    Valeur& vpsi = psi.va ; 

//==========================================================================
//			    Start of iteration 
//==========================================================================

    Tbl tdiff(nz) ; 
    double diff ; 
    niter = 0 ; 
 
    
    do {
    
	// Computation of the source for sol_poisson_compact
	// -------------------------------------------------

	b_grad_psi = bb(0) % psi.dsdr() + 
			 bb(1) % psi.srdsdt() +
			 bb(2) % psi.srstdsdp() ; 
    

	vpsi.ylm() ;	// Expansion of psi onto spherical harmonics

	// Lap_xi(psi) :

	dpsi = vpsi.dsdx() ; 
	dpsi.ylm() ;		    // necessary because vpsi.p_dsdx
				    // has been already computed (in 
				    // non-spherical harmonics bases) by
				    // the call psi.dsdr() above

	aux_psi = 2.*dpsi + (vpsi.lapang()).sx() ;

	d2psi = vpsi.d2sdx2() ; 
	d2psi.ylm() ; 
	
	lap_xi_psi = d2psi + aux_psi.sx() ; 

	// Effective source : 

	sour_j = source.va 
		    + alph * ( lap_xi_psi - (lap_xi_psi.mult_x()).mult_x() )
		    - aa.va   * (psi.laplacien()).va
		    + bet * dpsi.mult_x() 
		    - b_grad_psi.va ; 
			   
	sour_j.std_base_scal() ; 
	sour_j.coef() ; 
	sour_j.ylm() ;	    // Spherical harmonics expansion 

	// The term l=0 of the effective source is set to zero : 
	// ---------------------------------------------------
	double somlzero = 0 ; 
    	for (int i=0; i<mg->get_nr(0); i++) {
	    somlzero += fabs( (*(sour_j.c_cf))(0, 0, 0, i) ) ; 
	    (sour_j.c_cf)->set(0, 0, 0, i) = 0 ;  
	}
	
	if (somlzero > 1.e-10) {
	    cout << "### WARNING : Map_radial::poisson_compact : " << endl
		 << " the l=0 part of the effective source is > 1.e-10  : "
		 << somlzero << endl ; 
	}

	
	// Resolution of the equation 
	//   a (1-xi^2) Lap_xi(psi) + b xi dpsi/dxi = sour_j 
	//---------------------------------------------------
	
	bool reamorce = (niter == 0) ;
	
	assert(sour_j.c_cf != 0x0) ; 
	
	psi.set_etat_zero() ;  // to call Cmp::del_deriv().
	psi.set_etat_qcq() ;
	vpsi = sol_poisson_compact( *(sour_j.c_cf), alph, bet, reamorce ) ;


	// Test: multiplication by the operator matrix 
	// -------------------------------------------
	
//	Mtbl_cf cvpsi = *(vpsi.c_cf) ; 
//	Mtbl_cf csour = *(sour_j.c_cf) ; 
//	
//	int np = mg->get_np(0) ; 
//	int nt = mg->get_nt(0) ; 
//	int nr = mg->get_nr(0) ;
//	
//	for (int k=0 ; k<np+1 ; k++) {
//	    for (int j=0 ; j<nt ; j++) {
//		if (nullite_plm(j, nt, k, np, cvpsi.base) == 1) {
//	    int l_quant, m_quant, base_r ; 
//	    donne_lm(nz, 0, j, k, cvpsi.base, m_quant, l_quant, base_r) ;
//	    
//	    cout << "### k, j, l, m : " << k << "  " << j << "    "
//		    <<  l_quant << "  " << m_quant << endl ; 
//	    Matrice operateur = alph * laplacien_nul_mat(nr, l_quant, base_r)
//			         + bet * xdsdx_mat(nr, l_quant, base_r) ;
//	    Tbl coef(nr) ; 	
//	    operateur = combinaison_nul(operateur, l_quant, coef, base_r) ;
//	    
//		    Tbl so(nr) ;
//		    so.set_etat_qcq() ;
//		    for (int i=0 ; i<nr ; i++)
//			so.set(i) = csour(0, k, j, i) ;
//		    so = combinaison_nul(so, l_quant, coef, base_r) ;
//
//		    Tbl psi1(nr) ;
//		    Tbl opi1(nr) ;
//		    psi1.set_etat_qcq() ; 
//		    opi1.set_etat_qcq() ; 
//		    
//	    bool discrep = false ; 
//
//	    for (int i=0; i<nr; i++) {
//		
//		double op = 0 ; 
//		for (int i1=0; i1<nr; i1++) {
//		    op += operateur(i, i1) * cvpsi(0, k, j, i1) ; 
//		}
//		
//		psi1.set(i) = cvpsi(0, k, j, i) ; 
//		opi1.set(i) = op ; 
//		
//		double x1 = so(i) ; 
//		double x2 = op - so(i) ; 
//		double seuil = 1e-11 ; 
//		if ( fabs(x2) > seuil )
//		    {
//		    discrep = true ; 
//		    cout << "i : op , sou,  diff : " << i << " : " << op << "  " 
//		     << x1 << "  " << x2 << endl ; 
//		     } 
//		
//	    }
//
//	    if ( discrep ) {
//		cout << "Matrice pour k,  j = " << k << "  " << j << endl ; 
//		cout << operateur << endl ; 
//		cout << "psi : " << psi1 << endl ; 
//		cout << "op(psi) : " << opi1 << endl ; 
//		cout << "so : " << so << endl ; 
//	    }
//
//		}   // fin du cas ou m<=l 
//	    }	// fin de boucle sur j	
//	}   // fin de boucle sur k 


	// Test: has the equation been correctly solved ?
	// ---------------------------------------------
	
	// Lap_xi(psi) :
	aux_psi = vpsi.dsdx() ; 
	
	aux_psi = 2.*aux_psi + (vpsi.lapang()).sx() ;
		
	lap_xi_psi = vpsi.d2sdx2() + aux_psi.sx() ; 
			
	oper_psi = alph * ( lap_xi_psi - (lap_xi_psi.mult_x()).mult_x() )
		   + bet * (vpsi.dsdx()).mult_x() ; 
			

	double maxc = fabs( max(*(vpsi.c_cf))(0) ) ; 
	double minc = fabs( min(*(vpsi.c_cf))(0) ) ; 
	double max_abs_psi = ( maxc > minc ) ? maxc : minc ; 	

	maxc = fabs( max(*(sour_j.c_cf))(0) ) ; 
	minc = fabs( min(*(sour_j.c_cf))(0) ) ; 
	double max_abs_sou = ( maxc > minc ) ? maxc : minc ; 	

	Mtbl_cf diff_opsou = *oper_psi.c_cf - *sour_j.c_cf ; 
	maxc = fabs( max(diff_opsou)(0) ) ; 
	minc = fabs( min(diff_opsou)(0) ) ; 
	double max_abs_diff = ( maxc > minc ) ? maxc : minc ; 	


//	cout << " Coef of oper_psi - sour_j : " << endl ; 
//	diff_opsou.affiche_seuil(cout, 4, 1.e-11) ; 

	cout << "  Step " << niter 
	     << "  : test (1-xi^2) Lap_xi(psi) + b xi dpsi/dxi : "  << endl ; 
	cout << "   max |psi| |sou| |oper(psi)-sou|: " 
	     << max_abs_psi << "  " << max_abs_sou << "  "  
	     << max_abs_diff << endl ; 
	
	// Relaxation : 
	// -----------
	
	vpsi.ylm_i() ;    // Inverse spherical harmonics transform

	psi = relax * psi + unmrelax * psi_jm1 ; 

	tdiff = diffrel(psi, psi_jm1) ; 

	diff = max(tdiff) ;

	cout << 
	"   Relative difference psi^J <-> psi^{J-1} :                       " 
	     << tdiff(0) << endl ; 

//	arrete() ; 

	// Updates for the next iteration
	// -------------------------------    

	psi_jm1 = psi ;
	niter++ ; 

    }	// End of iteration 
    while ( (diff > precis) && (niter < nitermax) ) ;
    
//==========================================================================
//			    End of iteration 
//==========================================================================        
					 
}



	/////////////////////////////
	//	Multi domain version  //
	/////////////////////////////


void Map_radial::poisson_compact(int nzet, const Cmp& source, const Cmp& aa, 
				 const Tenseur& bb, const Param& par, 
				 Cmp& psi) const {
	
	if (nzet == 1) {
		poisson_compact(source, aa, bb, par, psi) ;
		return ; 
	}


    // Protections
    // -----------
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(aa.get_etat() != ETATNONDEF) ; 
    assert(bb.get_etat() != ETATNONDEF) ; 
    assert(aa.get_mp() == source.get_mp()) ; 
    assert(bb.get_mp() == source.get_mp()) ; 
    assert(psi.get_mp() == source.get_mp()) ; 
    
    
    // The components of vector b must be given in the spherical basis
    //  associated with the mapping : 
    assert(*(bb.get_triad()) == bvect_spher) ; 

    // Maybe nothing to do ?
    // ---------------------
    
	if ( source.get_etat() == ETATZERO ) {	
		psi.set_etat_zero() ; 
		return ; 
	}	
	
    // Computation parameters
    // ----------------------
    int nitermax = par.get_int() ; 
    int& niter = par.get_int_mod() ; 
    double precis = par.get_double(0) ; 
    double relax = par.get_double(1) ; 
    double unmrelax = 1. - relax ; 

	// Auxiliary affine mapping
    Map_af mpaff(*this) ; 

	// Coefficients to fit the profiles of aa and bb in each domain
	// ------------------------------------------------------------
	Tbl ac(nzet,3) ; 
	ac.annule_hard() ;	// initialization to zero 
	Tbl bc(nzet,3) ; 
	bc.annule_hard() ; // initialization to zero 

	Valeur ap = aa.va ;
	Valeur bp = bb(0).va ;

	// Coefficients in the nucleus
	int nrm1 = mg->get_nr(0) - 1 ; 
	ac.set(0,0) = ap(0,0,0,0) ;
	ac.set(0,2) = ap(0,0,0,nrm1) - ac(0,0) ; 
	
	bc.set(0,1) = bp(0,0,0,nrm1) ; 

	// Coefficients in the intermediate shells
	for (int lz=1; lz<nzet-1; lz++) {
		nrm1 = mg->get_nr(lz) - 1 ; 
		ac.set(lz,0) = 0.5 * ( ap(lz,0,0,nrm1) + ap(lz,0,0,0) ) ; 
		ac.set(lz,1) = 0.5 * ( ap(lz,0,0,nrm1) - ap(lz,0,0,0) ) ; 

		bc.set(lz,0) = 0.5 * ( bp(lz,0,0,nrm1) + bp(lz,0,0,0) ) ; 
		bc.set(lz,1) = 0.5 * ( bp(lz,0,0,nrm1) - bp(lz,0,0,0) ) ; 
	}

	// Coefficients in the external shell
	int lext = nzet - 1 ; 
	nrm1 = mg->get_nr(lext) - 1 ; 
	ac.set(lext,0) = 0.5 * ap(lext,0,0,0) ; 
	ac.set(lext,1) = - ac(lext,0) ; 

	bc.set(lext,0) = 0.5 * ( bp(lext,0,0,nrm1) + bp(lext,0,0,0) ) ; 
	bc.set(lext,1) = 0.5 * ( bp(lext,0,0,nrm1) - bp(lext,0,0,0) ) ;   

	cout << "ac : " << ac << endl ; 
	cout << "bc : " << bc << endl ; 

	// Prefactor of Lap_xi(Psi) and dPsi/dxi
	// -------------------------------------

	Mtbl ta(mg) ;
	Mtbl tb(mg) ;
	ta.annule_hard() ; 
	tb.annule_hard() ; 
	for (int lz=0; lz<nzet; lz++) {
		const double* xi = mg->get_grille3d(lz)->x ;
		double* tta = ta.set(lz).t ; 
		double* ttb = tb.set(lz).t ; 
		int np = mg->get_np(lz) ; 
		int nt = mg->get_nt(lz) ; 
		int nr = mg->get_nr(lz) ; 
		int pt = 0 ; 
		for (int k=0; k<np; k++) {
			for (int j=0; j<nt; j++) {
				for (int i=0; i<nr; i++) {
					tta[pt] = ac(lz,0) + xi[i] * (ac(lz,1) + ac(lz,2) * xi[i]) ;
					ttb[pt] = bc(lz,0) + xi[i] * (bc(lz,1) + bc(lz,2) * xi[i]) ;
					pt++ ;  
				}
			}
		}
	}

// Verification
// 	cout << "Map :" << *(aa.get_mp()) << endl ; 
// 	Cmp tverif(*this) ; 
// 	tverif = ap ;
// 	tverif.std_base_scal() ;  
// 	des_profile(tverif,0., 4., 0., 0.) ; 
// 	tverif = ta ; 
// 	tverif.std_base_scal() ;  
// 	des_profile(tverif,0., 4., 0., 0.) ; 
// 
// 	tverif = bp ;
// 	tverif.std_base_scal() ;  
// 	des_profile(tverif,0., 4., 0., 0.) ; 
// 	tverif = tb ; 
// 	tverif.std_base_scal() ;  
// 	des_profile(tverif,0., 4., 0., 0.) ; 


    // Everything is set to zero except inside the star
    // -------------------------------------------------

    int nz = mg->get_nzone() ;

    psi.annule(nzet, nz-1) ; 

    // Auxilary quantities
    // -------------------
    Cmp psi_jm1 = psi ; 
    Cmp b_grad_psi(this) ; 
    Valeur sour_j(*mg) ; 
    Valeur aux_psi(*mg) ;		
    Valeur lap_xi_psi(*mg) ; 
    Valeur oper_psi(*mg) ; 
    Valeur dpsi(*mg) ; 
    Valeur d2psi(*mg) ; 

    Valeur& vpsi = psi.va ; 

//==========================================================================
//			    Start of iteration 
//==========================================================================

    Tbl tdiff(nz) ; 
    double diff ; 
    niter = 0 ; 
 
    
    do {
    
	// Computation of the source for sol_poisson_compact
	// -------------------------------------------------

 	b_grad_psi = bb(0) % psi.dsdr() + bb(1) % psi.srdsdt() + bb(2) % psi.srstdsdp() ; 

//?? 
 	vpsi.ylm() ;	// Expansion of psi onto spherical harmonics

	// Effective source : 

	Cmp lap_zeta(mpaff) ; 
	mpaff.laplacien(psi, 0, lap_zeta) ; 

	Cmp grad_zeta(mpaff) ;
	mpaff.dsdr(psi, grad_zeta) ; 

	sour_j = source.va 
			+ ta * lap_zeta.va - aa.va * (psi.laplacien()).va
		    + tb * grad_zeta.va - b_grad_psi.va ; 
			   
	sour_j.std_base_scal() ; 
	sour_j.coef() ; 
	sour_j.ylm() ;	    // Spherical harmonics expansion 


	// The term l=0 of the effective source is set to zero : 
	// ---------------------------------------------------

	for (int lz=0; lz<nzet; lz++) {
		double somlzero = 0 ; 
    	for (int i=0; i<mg->get_nr(lz); i++) {
	    	somlzero += fabs( (*(sour_j.c_cf))(lz, 0, 0, i) ) ; 
	    	(sour_j.c_cf)->set(lz, 0, 0, i) = 0 ;  
		}
		if (somlzero > 1.e-10) {
	    	cout << "### WARNING : Map_radial::poisson_compact : " << endl
		 	<< " domain no. " << lz << " : the l=0 part of the effective source is > 1.e-10  : "
		 	<< somlzero << endl ; 
		}
	}

	// Resolution of the equation 
	//----------------------------
	
	bool reamorce = (niter == 0) ;
	
	assert(sour_j.c_cf != 0x0) ; 
	
	psi.set_etat_zero() ;  // to call Cmp::del_deriv().
	psi.set_etat_qcq() ;
	
	vpsi = sol_poisson_compact(mpaff, *(sour_j.c_cf), ac, bc, reamorce) ;

	// Test: has the equation been correctly solved ?
	// ---------------------------------------------
	
	mpaff.laplacien(psi, 0, lap_zeta) ; 
	mpaff.dsdr(psi, grad_zeta) ; 

	oper_psi = ta * lap_zeta.va + tb * grad_zeta.va ; 
	oper_psi.std_base_scal() ; 
	oper_psi.coef() ; 
	oper_psi.ylm() ;			

	Mtbl_cf diff_opsou = *oper_psi.c_cf - *sour_j.c_cf ; 
		//cout << " Coef of oper_psi - sour_j : " << endl ; 
		// diff_opsou.affiche_seuil(cout, 4, 1.e-11) ; 
	
	cout << "poisson_compact:  step " << niter << " : " << endl ;  
	for (int lz=0; lz<nzet; lz++) {
		double maxc = fabs( max(*(vpsi.c_cf))(lz) ) ; 
		double minc = fabs( min(*(vpsi.c_cf))(lz) ) ; 
		double max_abs_psi = ( maxc > minc ) ? maxc : minc ; 	

		maxc = fabs( max(*(sour_j.c_cf))(lz) ) ; 
		minc = fabs( min(*(sour_j.c_cf))(lz) ) ; 
		double max_abs_sou = ( maxc > minc ) ? maxc : minc ; 	

		maxc = fabs( max(diff_opsou)(lz) ) ; 
		minc = fabs( min(diff_opsou)(lz) ) ; 
		double max_abs_diff = ( maxc > minc ) ? maxc : minc ; 	

		cout << "  lz = " << lz << " : max |psi| |sou| |oper(psi)-sou|: " 
	     << max_abs_psi << "  " << max_abs_sou << "  "  
	     << max_abs_diff << endl ; 
	}

	// Relaxation : 
	// -----------
	
	vpsi.ylm_i() ;    // Inverse spherical harmonics transform

	psi = relax * psi + unmrelax * psi_jm1 ; 

	tdiff = diffrel(psi, psi_jm1) ; 

	diff = max(tdiff) ;

	cout << "   Relative difference psi^J <-> psi^{J-1} :    " 
	     << tdiff << endl ;

	// Updates for the next iteration
	// -------------------------------    

	psi_jm1 = psi ;
	niter++ ; 

    }	// End of iteration 
    while ( (diff > precis) && (niter < nitermax) ) ;
    
//==========================================================================
//			    End of iteration 
//==========================================================================        

	psi.annule(nzet, nz-1) ; 

}


















