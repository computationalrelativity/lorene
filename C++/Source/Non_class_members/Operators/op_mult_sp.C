/*
 * Sets of routines for multiplication by sin(phi)
 * (internal use)
 *
 *	void _mult_sp_XXXX(Tbl * t, int & b)
 *	t	pointer on the Tbl containing the spectral coefficients
 *	b	spectral base
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char op_mult_sp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:53:25  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2007/12/14 10:19:33  jl_cornou
 * *** empty log message ***
 *
 * Revision 1.3  2007/10/05 12:37:20  j_novak
 * Corrected a few errors in the theta-nonsymmetric case (bases T_COSSIN_C and
 * T_COSSIN_S).
 *
 * Revision 1.2  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/11/14  15:12:01  eric
 * Traitement du cas np=1
 *
 * Revision 2.3  2000/09/18  10:14:59  eric
 * Ajout des bases P_COSSIN_P et P_COSSIN_I
 *
 * Revision 2.2  2000/09/12  13:36:38  eric
 * Met les bonnes bases meme dans le cas ETATZERO
 *
 * Revision 2.1  2000/09/12  08:29:11  eric
 * Traitement des bases qui dependent de la parite de m
 *
 * Revision 2.0  2000/09/11  13:53:42  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tbl.h"

		//-----------------------------------
		//    Routine for unknown cases
		//-----------------------------------

namespace Lorene {
void _mult_sp_pas_prevu(Tbl* , int& base) {
    cout << "mult_sp() is not not implemented for the basis " << base << " !"
       << endl ;
    abort() ;
}

		    //----------------
		    // basis P_COSSIN
		    //----------------

void _mult_sp_p_cossin(Tbl* tb, int& base) {
    
    assert(tb->get_etat() != ETATNONDEF) ;	// Protection

    // The spectral bases in r and theta which depend on the parity of m 
    // are changed
    // -----------------------------------------------------------------
    
    int base_r = base & MSQ_R ; 
    int base_t = base & MSQ_T ; 
    const int base_p = base & MSQ_P ; 
    
    switch (base_r) {

	case R_CHEBPIM_P : {
	    base_r = R_CHEBPIM_I ; 
	    break ; 
	}
	
	case R_CHEBPIM_I : {
	    base_r = R_CHEBPIM_P ; 
	    break ; 	    
	}

        case R_CHEBPI_P : {
	    break ; 
	}
	
	case R_CHEBPI_I : {
	    break ; 	    
	}  
	  
	case R_CHEB : {
	    break ; 
	}
	
	case R_JACO02 : {
	    break ;
	}
	
	case R_CHEBU : {
	    break ; 
	}
	
	default : {
	    cout << "_mult_cp_p_cossin : unknown basis in r !" << endl ; 
	    cout << "  base_r = " << hex << base_r << endl ; 
	    abort() ; 
	}	
    }
    
    switch (base_t) {

	case T_COSSIN_CP : {
	    base_t = T_COSSIN_SI ; 
	    break ; 
	}
	
	case T_COSSIN_SI : {
	    base_t = T_COSSIN_CP ; 
	    break ; 
	}
	
	case T_COSSIN_CI : {
	    base_t = T_COSSIN_SP ; 
	    break ; 
	}
	
	case T_COSSIN_SP : {
	    base_t = T_COSSIN_CI ; 
	    break ; 
	}
	  
        case T_COSSIN_S : {
	    base_t = T_COSSIN_C ; 
	    break ; 
	}

	case T_COSSIN_C : {
	    base_t = T_COSSIN_S ; 
	    break ; 
	}  
	  
	default : {
	    cout << "_mult_cp_p_cossin : unknown basis in theta !" << endl ; 
	    cout << "  base_t = " << hex << base_t << endl ; 
	    abort() ; 
	}	
    }
    
    base = base_r | base_t | base_p ; 
    

    //----------------------------------
    //  Start of computation 
    //----------------------------------

    // Nothing to do ? 
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    assert(tb->get_etat() == ETATQCQ) ;	// Protection
  
    // Number of degrees of freedom
    int nr = tb->get_dim(0) ;
    int nt = tb->get_dim(1) ;
    int np = tb->get_dim(2) - 2 ;
    

    // Special case np = 1 (axisymmetry)  --> zero is returned
    // ---------------------------------
    
    if (np==1) {
	tb->set_etat_zero() ; 
	return ; 
    }


    assert(np >= 4) ; 
    
    int ntnr = nt*nr ;
    
    double* const cf = tb->t ;				// input coefficients
    double* const resu = new double[ tb->get_taille() ] ;	// final result
    double* co = resu ;		// initial position 
        
    // Case k=0 (m=0)
    // --------------
    
    int q = 3 * ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = 0.5 * cf[q + i] ; 
    }
    co += ntnr ; 
    
    // Case k = 1 
    // ----------
    
    for (int i=0; i<ntnr; i++) {
	co[i] = 0 ; 
    }
    co += ntnr ; 
    
    if (np==4) {

	// Case k = 2  for np=4 
	// ----------
    
	for (int i=0; i<ntnr; i++) {
	    co[i] = 0 ; 
	}
	co += ntnr ; 
	
	// Case k = 3  for np=4 
	// ----------
    
	q = 4*ntnr ; 
	for (int i=0; i<ntnr; i++) {
	    co[i] = cf[i] - 0.5 * cf[q+i] ; 
	}
	co += ntnr ; 
		
    }
    else{

	// Case k = 2  for np>4 
	// ----------
    
	q = 5*ntnr ; 
	for (int i=0; i<ntnr; i++) {
	    co[i] = 0.5 * cf[q+i] ; 
	}
	co += ntnr ; 
	
	// Case k = 3 for np>4
	// ----------

	q = 4*ntnr ; 
	for (int i=0; i<ntnr; i++) {
	    co[i] = cf[i] - 0.5 * cf[q+i] ; 
	}
	co += ntnr ; 
	
	
	// Cases 4 <= k <= np-3
	// --------------------
    
	for (int k=4; k<np-2; k+=2) {

	    int k1 = k ;	// k1 even
	    
	    int q1 = (k1-1)*ntnr ;
	    int q2 = (k1+3)*ntnr ; 
	    for (int i=0; i<ntnr; i++) {
		co[i] = 0.5 * ( cf[q2+i] - cf[q1+i] ) ; 
	    }
	    co += ntnr ; 
	    k1++ ;		// k1 odd
    
	    q1 = (k1-3)*ntnr ; 
	    q2 = (k1+1)*ntnr ; 
	    for (int i=0; i<ntnr; i++) {
		co[i] = 0.5 * ( cf[q1+i] - cf[q2+i] ) ; 
	    }
	    co += ntnr ; 
	}	    
	
	// Case k = np - 2  for np>4
	// ---------------
	
	q = (np-3)*ntnr ; 
	for (int i=0; i<ntnr; i++) {
	    co[i] = - 0.5 * cf[q + i] ; 
	}
	co += ntnr ; 
	
	// Case k = np - 1  for np>4
	// ---------------
	
	int q1 = (np-4)*ntnr ; 
	int q2 = np*ntnr ; 
	for (int i=0; i<ntnr; i++) {
	    co[i] = 0.5 * ( cf[q1+i] - cf[q2+i] ) ; 
	}
	co += ntnr ; 
	 	
    }  // End of case np > 4
    
    
    // Case k = np
    // -----------
    
    q = (np-1)*ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = - 0.5 * cf[q+i] ; 
    }
    co += ntnr ; 
    
    // Case k = np+1
    // -------------

    for (int i=0; i<ntnr; i++) {
	co[i] = 0 ; 
    }
    
    //## verif : 
    co += ntnr ;
    assert( co == resu + (np+2)*ntnr ) ; 
    
    
    // The result is set to tb :
    // ----------------------- 
    delete [] tb->t ;
    tb->t = resu ;
        
    return ; 
}

		    //------------------
		    // basis P_COSSIN_P
		    //------------------

void _mult_sp_p_cossin_p(Tbl* tb, int& base) {
    
    assert(tb->get_etat() != ETATNONDEF) ;	// Protection

    // New spectral bases
    // ------------------        
    int base_r = base & MSQ_R ; 
    int base_t = base & MSQ_T ; 
    base = base_r | base_t | P_COSSIN_I ; 
    
    //----------------------------------
    //  Start of computation 
    //----------------------------------

    // Nothing to do ? 
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    assert(tb->get_etat() == ETATQCQ) ;	// Protection
  
    // Number of degrees of freedom
    int nr = tb->get_dim(0) ;
    int nt = tb->get_dim(1) ;
    int np = tb->get_dim(2) - 2 ;
    
    // Special case np = 1 (axisymmetry)  --> zero is returned
    // ---------------------------------
    
    if (np==1) {
	tb->set_etat_zero() ; 
	return ; 
    }

    assert(np >= 4) ; 
    
    int ntnr = nt*nr ;
    
    double* const cf = tb->t ;				// input coefficients
    double* const resu = new double[ tb->get_taille() ] ;	// final result
    double* co = resu ;		// initial position 
        
    // Case k=0 
    // --------
    
    int q = 3 * ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = 0.5 * cf[q + i] ; 
    }
    co += ntnr ; 
    
    // Case k = 1 
    // ----------
    
    for (int i=0; i<ntnr; i++) {
	co[i] = 0 ; 
    }
    co += ntnr ; 
    
    // Case k = 2 
    // ----------
    
    q = 2*ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = cf[i] - 0.5 * cf[q+i] ; 
    }
    co += ntnr ; 
    
    // Cases 3 <= k <= np-2
    // --------------------

    for (int k=3; k<np-1; k+=2) {

        int k1 = k ;	// k1 odd
	    
        int q1 = k1*ntnr ;
        int q2 = (k1+2)*ntnr ; 
        for (int i=0; i<ntnr; i++) {
	    co[i] = 0.5 * ( cf[q2+i] - cf[q1+i] ) ; 
	}
	co += ntnr ; 
	k1++ ;		// k1 even
    
	q1 = (k1-2)*ntnr ; 
	q2 = k1*ntnr ; 
	for (int i=0; i<ntnr; i++) {
	    co[i] = 0.5 * ( cf[q1+i] - cf[q2+i] ) ; 
	}
	co += ntnr ; 
    }
	
    // Case k = np - 1  
    // ---------------
        
    q = (np-1)*ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = - 0.5 * cf[q+i] ; 
    }
    co += ntnr ; 
        
    // Case k = np
    // -----------
    
    int q1 = (np-2)*ntnr ; 
    int q2 = np*ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = 0.5 * ( cf[q1+i] - cf[q2+i] ) ;
    }
    co += ntnr ; 
    
    // Case k = np+1
    // -------------

    for (int i=0; i<ntnr; i++) {
	co[i] = 0 ; 
    }
    
    //## verif : 
    co += ntnr ;
    assert( co == resu + (np+2)*ntnr ) ; 
    
    
    // The result is set to tb :
    // ----------------------- 
    delete [] tb->t ;
    tb->t = resu ;
        
    return ; 
}


		    //------------------
		    // basis P_COSSIN_I
		    //------------------

void _mult_sp_p_cossin_i(Tbl* tb, int& base) {
    
    assert(tb->get_etat() != ETATNONDEF) ;	// Protection

    // New spectral bases
    // ------------------        
    int base_r = base & MSQ_R ; 
    int base_t = base & MSQ_T ; 
    base = base_r | base_t | P_COSSIN_P ; 
    
    //----------------------------------
    //  Start of computation 
    //----------------------------------

    // Nothing to do ? 
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    assert(tb->get_etat() == ETATQCQ) ;	// Protection
  
    // Number of degrees of freedom
    int nr = tb->get_dim(0) ;
    int nt = tb->get_dim(1) ;
    int np = tb->get_dim(2) - 2 ;
    
    // Special case np = 1 (axisymmetry)  --> zero is returned
    // ---------------------------------
    
    if (np==1) {
	tb->set_etat_zero() ; 
	return ; 
    }

    assert(np >= 4) ; 
    
    int ntnr = nt*nr ;
    
    double* const cf = tb->t ;				// input coefficients
    double* const resu = new double[ tb->get_taille() ] ;	// final result
    double* co = resu ;		// initial position 
        
    // Case k=0 
    // --------
    
    int q = 2 * ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = 0.5 * cf[q + i] ; 
    }
    co += ntnr ; 
    
    // Case k = 1 
    // ----------
    
    for (int i=0; i<ntnr; i++) {
	co[i] = 0 ; 
    }
    co += ntnr ; 
    
    // Case k = 2 
    // ----------
    
    int q1 = 2*ntnr ; 
    int q2 = 4*ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = 0.5 * ( cf[q2+i] - cf[q1+i] ) ; 
    }
    co += ntnr ; 
    
    // Case k = 3 
    // ----------
    
    q = 3*ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = 0.5 * ( cf[i] - cf[q+i] ) ; 
    }
    co += ntnr ; 
    
    // Cases 4 <= k <= np-1
    // --------------------

    for (int k=4; k<np; k+=2) {

        int k1 = k ;	// k1 even
	    
        q1 = k1*ntnr ;
        q2 = (k1+2)*ntnr ; 
        for (int i=0; i<ntnr; i++) {
	    co[i] = 0.5 * ( cf[q2+i] - cf[q1+i] ) ; 
	}
	co += ntnr ; 
	k1++ ;		// k1 odd
    
	q1 = (k1-2)*ntnr ; 
	q2 = k1*ntnr ; 
	for (int i=0; i<ntnr; i++) {
	    co[i] = 0.5 * ( cf[q1+i] - cf[q2+i] ) ; 
	}
	co += ntnr ; 
    }
	
    // Case k = np
    // -----------
        
    q = np*ntnr ; 
    for (int i=0; i<ntnr; i++) {
	co[i] = - 0.5 * cf[q+i] ; 
    }
    co += ntnr ; 
            
    // Case k = np+1
    // -------------

    for (int i=0; i<ntnr; i++) {
	co[i] = 0 ; 
    }
    
    //## verif : 
    co += ntnr ;
    assert( co == resu + (np+2)*ntnr ) ; 
    
    
    // The result is set to tb :
    // ----------------------- 
    delete [] tb->t ;
    tb->t = resu ;
        
    return ; 
}

}
