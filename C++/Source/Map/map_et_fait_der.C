/*
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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


char map_et_fait_der_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/13 08:53:04  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:13:13  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2013/06/05 15:10:42  j_novak
 * Suppression of FINJAC sampling in r. This Jacobi(0,2) base is now
 * available by setting colloc_r to BASE_JAC02 in the Mg3d constructor.
 *
 * Revision 1.4  2008/08/27 08:49:01  jl_cornou
 * Added R_JACO02 case
 *
 * Revision 1.3  2003/10/15 10:40:26  e_gourgoulhon
 * Added new Coord's drdt and stdrdp.
 * Changed cast (const Map_et *) into static_cast<const Map_et*>.
 *
 * Revision 1.2  2002/10/16 14:36:41  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  1999/12/20  17:27:37  eric
 * Modif lapr_tp.
 *
 * Revision 1.1  1999/11/24  11:23:06  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Includes
#include <cassert>
#include <cstdlib>
#include <cmath>

#include "map.h"


//////////////////////////////////////////////////////////////////////////
//		Derivees d'ordre 1 du changement de variables		//
//////////////////////////////////////////////////////////////////////////

/*
 ************************************************************************
 *			    1/(dR/dx)	    ( -1/(dU/dx) ds la ZEC )
 ************************************************************************
 */

namespace Lorene {
Mtbl* map_et_fait_dxdr(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& da = *((cv->daa)[l]) ; 
	const Tbl& db = *((cv->dbb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case FIN: case RARE: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = 1. / 
				  ( alpha[l] * ( 1. + da(i) * ff(l, k, j, 0)
						    + db(i) * gg(l, k, j, 0) ) 
				  ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }
	        
	    case UNSURR: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = - 1. / 
				    ( alpha[l] * ( 1. + da(i) * ff(l, k, j, 0) )
				    ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_dxdr: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
}

/*
 *******************************************************************************
 *			     dR/dtheta	    ( dans la ZEC: - dU/dth )
 *******************************************************************************
 */

Mtbl* map_et_fait_drdt(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    const Valeur& dffdt = ff.dsdt() ; 
    const Valeur& dggdt = gg.dsdt() ; 

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE : case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = alpha[l] * (   aa(i) * dffdt(l, k, j, 0)
					    + bb(i) * dggdt(l, k, j, 0) ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = - aa(i) * dffdt(l, k, j, 0)	;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_drdt: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 


/*
 *******************************************************************************
 *			  1/sin(theta)   dR/dphi	    ( dans la ZEC: - 1/sin(th) dU/dth )
 *******************************************************************************
 */

Mtbl* map_et_fait_stdrdp(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    const Valeur& stdffdp = ff.stdsdp() ; 
    const Valeur& stdggdp = gg.stdsdp() ; 

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE : case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = alpha[l] * (   aa(i) * stdffdp(l, k, j, 0)
					    + bb(i) * stdggdp(l, k, j, 0) ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = - aa(i) * stdffdp(l, k, j, 0)	;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_stdrdp: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 


/*
 *******************************************************************************
 *			    1/R dR/dtheta	    ( dans la ZEC: - 1/U dU/dth )
 *******************************************************************************
 */

Mtbl* map_et_fait_srdrdt(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const double* beta = cv->beta ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    const Valeur& dffdt = ff.dsdt() ; 
    const Valeur& dggdt = gg.dsdt() ; 

    
    for (int l=0 ; l<nz ; l++) {

	const Grille3d* g = mg->get_grille3d(l) ;

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    
	    const Tbl& asx = cv->aasx ; 
	    const Tbl& bsx = cv->bbsx ; 
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r =	(   asx(i) * dffdt(l, k, j, 0) 
				  + bsx(i) * dggdt(l, k, j, 0)  ) /
				( 1. + asx(i) * ff(l, k, j, 0)
				     + bsx(i) * gg(l, k, j, 0) ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }

	    case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = alpha[l] * (   aa(i) * dffdt(l, k, j, 0)
					    + bb(i) * dggdt(l, k, j, 0) )
				/ ( alpha[l] * (
					(g->x)[i] + aa(i) * ff(l, k, j, 0)
						  + bb(i) * gg(l, k, j, 0) ) 
				      + beta[l] ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    const Tbl& asxm1 = cv->zaasx ; 

	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = - asxm1(i) * dffdt(l, k, j, 0)
				/  ( 1. + asxm1(i) * ff(l, k, j, 0) )	;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_srdrdt: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 

/*
 *****************************************************************************
 *	1/(R sin(theta)) dR/dphi	( ds la ZEC: - 1/(U sin(th)) dU/dphi )
 *****************************************************************************
 */

Mtbl* map_et_fait_srstdrdp(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const double* beta = cv->beta ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    const Valeur& stdffdp = ff.stdsdp() ; 
    const Valeur& stdggdp = gg.stdsdp() ; 

    
    for (int l=0 ; l<nz ; l++) {

	const Grille3d* g = mg->get_grille3d(l) ;

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    
	    const Tbl& asx = cv->aasx ; 
	    const Tbl& bsx = cv->bbsx ; 
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r =	(   asx(i) * stdffdp(l, k, j, 0) 
				  + bsx(i) * stdggdp(l, k, j, 0)  ) /
				( 1. + asx(i) * ff(l, k, j, 0)
				     + bsx(i) * gg(l, k, j, 0) ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }

	    case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = alpha[l] * (   aa(i) * stdffdp(l, k, j, 0)
					    + bb(i) * stdggdp(l, k, j, 0) )
				/ ( alpha[l] * (
					(g->x)[i] + aa(i) * ff(l, k, j, 0)
						  + bb(i) * gg(l, k, j, 0) ) 
				      + beta[l] ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    const Tbl& asxm1 = cv->zaasx ; 

	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = - asxm1(i) * stdffdp(l, k, j, 0)
				/  ( 1. + asxm1(i) * ff(l, k, j, 0) )	;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_srstdrdp: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    return mti ; 
} 

/*
 *******************************************************************************
 *			    1/R^2 dR/dtheta	( dans la ZEC: - 1/U^2 dU/dth )
 *******************************************************************************
 */

Mtbl* map_et_fait_sr2drdt(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const double* beta = cv->beta ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    const Valeur& dffdt = ff.dsdt() ; 
    const Valeur& dggdt = gg.dsdt() ; 

    
    for (int l=0 ; l<nz ; l++) {

	const Grille3d* g = mg->get_grille3d(l) ;

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    
	    const Tbl& asx = cv->aasx ; 
	    const Tbl& bsx = cv->bbsx ; 
	    const Tbl& asx2 = cv->aasx2 ; 
	    const Tbl& bsx2 = cv->bbsx2 ; 
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			double ww =  1. + asx(i) * ff(l, k, j, 0)
				        + bsx(i) * gg(l, k, j, 0) ;
		    
			*p_r =	(   asx2(i) * dffdt(l, k, j, 0) 
				  + bsx2(i) * dggdt(l, k, j, 0)  ) /
				(alpha[l] * ww*ww) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }

		    
	    case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			double ww = alpha[l] * (
					(g->x)[i] + aa(i) * ff(l, k, j, 0)
						  + bb(i) * gg(l, k, j, 0) ) 
				      + beta[l] ; 

			*p_r = alpha[l] * (   aa(i) * dffdt(l, k, j, 0)
					    + bb(i) * dggdt(l, k, j, 0) )
				/ ( ww*ww ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    const Tbl& asxm1 = cv->zaasx ; 
	    const Tbl& asxm1car = cv->zaasx2 ; 

	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			double ww = 1. + asxm1(i) * ff(l, k, j, 0) ;
			
			*p_r = - asxm1car(i) * dffdt(l, k, j, 0)
				/  ( alpha[l] * ww*ww )	;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_sr2drdt: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 

/*
 *******************************************************************************
 *	  1/(R^2 sin(theta)) dR/dphi	( ds la ZEC: - 1/(U^2 sin(th)) dU/dphi )
 *******************************************************************************
 */

Mtbl* map_et_fait_sr2stdrdp(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const double* beta = cv->beta ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    const Valeur& stdffdp = ff.stdsdp() ; 
    const Valeur& stdggdp = gg.stdsdp() ; 

    
    for (int l=0 ; l<nz ; l++) {

	const Grille3d* g = mg->get_grille3d(l) ;

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    
	    const Tbl& asx = cv->aasx ; 
	    const Tbl& bsx = cv->bbsx ; 
	    const Tbl& asx2 = cv->aasx2 ; 
	    const Tbl& bsx2 = cv->bbsx2 ; 
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			double ww =  1. + asx(i) * ff(l, k, j, 0)
				        + bsx(i) * gg(l, k, j, 0) ;
		    
			*p_r =	(   asx2(i) * stdffdp(l, k, j, 0) 
				  + bsx2(i) * stdggdp(l, k, j, 0)  ) /
				(alpha[l] * ww*ww) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }

	    case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			double ww = alpha[l] * (
					(g->x)[i] + aa(i) * ff(l, k, j, 0)
						  + bb(i) * gg(l, k, j, 0) ) 
				      + beta[l] ; 

			*p_r = alpha[l] * (   aa(i) * stdffdp(l, k, j, 0)
					    + bb(i) * stdggdp(l, k, j, 0) )
				/ ( ww*ww ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    const Tbl& asxm1 = cv->zaasx ; 
	    const Tbl& asxm1car = cv->zaasx2 ; 

	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			double ww = 1. + asxm1(i) * ff(l, k, j, 0) ;
			
			*p_r = - asxm1car(i) * stdffdp(l, k, j, 0)
				/  ( alpha[l] * ww*ww )	;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_sr2stdrdp: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 

/*
 ******************************************************************************
 *		1/(dR/dx) R/x			    ds. le noyau 
 *		1/(dR/dx) R/(x + beta_l/alpha_l)    ds. les coquilles
 *		1/(dU/dx) U/(x-1)		    ds. la ZEC
 ******************************************************************************
 */

Mtbl* map_et_fait_rsxdxdr(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const double* beta = cv->beta ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    
    for (int l=0 ; l<nz ; l++) {

	const Grille3d* g = mg->get_grille3d(l) ;

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 
	const Tbl& da = *((cv->daa)[l]) ; 
	const Tbl& db = *((cv->dbb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    
	    const Tbl& asx = cv->aasx ; 
	    const Tbl& bsx = cv->bbsx ; 
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			*p_r = ( 1. + asx(i) * ff(l, k, j, 0) 
				    + bsx(i) * gg(l, k, j, 0) ) / 
			       ( 1. + da(i) * ff(l, k, j, 0) 
				    + db(i) * gg(l, k, j, 0) ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }

		    
	    case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			*p_r = ( (g->x)[i] + aa(i) * ff(l, k, j, 0) 
					   + bb(i) * gg(l, k, j, 0) 
					   + beta[l]/alpha[l]	     )  /
			       ( ( 1. + da(i) * ff(l, k, j, 0) 
				      + db(i) * gg(l, k, j, 0) ) *
				 ( (g->x)[i] + beta[l]/alpha[l] ) ) ;			
			
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    const Tbl& asxm1 = cv->zaasx ; 

	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			*p_r = ( 1. + asxm1(i) * ff(l, k, j, 0) ) /
			       ( 1. + da(i) * ff(l, k, j, 0)	) ; 
			
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_rsxdxdr: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 

/*
 ******************************************************************************
 * [ R/ (alpha_l x + beta_l) ]^2 (dR/dx) /alpha_l   ds. le noyau et les coquilles
 * dU/dx /alpha_l				    ds. la ZEC
 ******************************************************************************
 */

Mtbl* map_et_fait_rsx2drdx(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const double* beta = cv->beta ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    
    for (int l=0 ; l<nz ; l++) {

	const Grille3d* g = mg->get_grille3d(l) ;

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 
	const Tbl& da = *((cv->daa)[l]) ; 
	const Tbl& db = *((cv->dbb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    
	    const Tbl& asx = cv->aasx ; 
	    const Tbl& bsx = cv->bbsx ; 
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			double ww =  1. + asx(i) * ff(l, k, j, 0)
				        + bsx(i) * gg(l, k, j, 0) ;

			*p_r = ww * ww * 
			       ( 1. + da(i) * ff(l, k, j, 0) 
				    + db(i) * gg(l, k, j, 0) ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }

		    
	    case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			double ww = ( (g->x)[i] + aa(i) * ff(l, k, j, 0)
						+ bb(i) * gg(l, k, j, 0) 
						+ beta[l]/alpha[l]	  )  /
				    ( (g->x)[i] + beta[l]/alpha[l] ) ; 

			*p_r = ww * ww * 
			       ( 1. + da(i) * ff(l, k, j, 0) 
				    + db(i) * gg(l, k, j, 0) ) ;			
			
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			*p_r =  1. + da(i) * ff(l, k, j, 0) ;
			
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_rsx2drdx: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 



//////////////////////////////////////////////////////////////////////////
//		Derivees d'ordre 2 du changement de variables		//
//////////////////////////////////////////////////////////////////////////


/*
 *******************************************************************************
 *			    d^2R/dx^2	    ( dans la ZEC: - d^2U/dx^2 )
 *******************************************************************************
 */

Mtbl* map_et_fait_d2rdx2(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& dda = *((cv->ddaa)[l]) ; 
	const Tbl& ddb = *((cv->ddbb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case FIN: case RARE: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			*p_r = alpha[l] * (   dda(i) * ff(l, k, j, 0)
					    + ddb(i) * gg(l, k, j, 0) ) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }
	        
	    case UNSURR: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			*p_r = -  alpha[l] * dda(i) * ff(l, k, j, 0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    } 

	    default: {
	    cout << "map_et_fait_d2rdx2: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 

/*
 *****************************************************************************
 *  1/R^2 (  1/sin(th) d/dth( sin(th) dR/dth ) + 1/sin(th)^2 d^2R/dphi^2  )		    
 *
 * [ dans la ZEC : 
 *  - 1/U^2 (  1/sin(th) d/dth( sin(th) dU/dth ) + 1/sin(th)^2 d^2U/dphi^2  )  ]		    
 *****************************************************************************
 */

Mtbl* map_et_fait_lapr_tp(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const double* beta = cv->beta ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 

    Valeur ff_tmp = ff ; 
    Valeur gg_tmp = gg ; 
    ff_tmp.ylm() ; 
    gg_tmp.ylm() ; 
    const Valeur& lapff = ff_tmp.lapang() ; 
    const Valeur& lapgg = gg_tmp.lapang() ; 

    
    for (int l=0 ; l<nz ; l++) {

	const Grille3d* g = mg->get_grille3d(l) ;

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    
	    const Tbl& asx = cv->aasx ; 
	    const Tbl& bsx = cv->bbsx ; 
	    const Tbl& asx2 = cv->aasx2 ; 
	    const Tbl& bsx2 = cv->bbsx2 ; 
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			double ww =  1. + asx(i) * ff(l, k, j, 0)
				        + bsx(i) * gg(l, k, j, 0) ;
		    
			*p_r =	(   asx2(i) * lapff(l, k, j, 0) 
				  + bsx2(i) * lapgg(l, k, j, 0)  ) /
				(alpha[l] * ww*ww) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }

	    case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			double ww = alpha[l] * (
					(g->x)[i] + aa(i) * ff(l, k, j, 0)
						  + bb(i) * gg(l, k, j, 0) ) 
				      + beta[l] ; 

			*p_r = alpha[l] * (   aa(i) * lapff(l, k, j, 0)
					    + bb(i) * lapgg(l, k, j, 0) )
				/ ( ww*ww ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    const Tbl& asxm1 = cv->zaasx ; 
	    const Tbl& asxm1car = cv->zaasx2 ; 

	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			double ww = 1. + asxm1(i) * ff(l, k, j, 0) ;
			
			*p_r = - asxm1car(i) * lapff(l, k, j, 0)
				/  ( alpha[l] * ww*ww )	;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_lapr_tp: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 

/*
 *******************************************************************************
 *			    d^2R/dthdx		( dans la ZEC: - d^2U/dthdx )
 *******************************************************************************
 */

Mtbl* map_et_fait_d2rdtdx(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    const Valeur& dffdt = ff.dsdt() ; 
    const Valeur& dggdt = gg.dsdt() ; 
    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& da = *((cv->daa)[l]) ; 
	const Tbl& db = *((cv->dbb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case FIN: case RARE: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			*p_r = alpha[l] * (   da(i) * dffdt(l, k, j, 0)
					    + db(i) * dggdt(l, k, j, 0) ) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }
	        
	    case UNSURR: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			*p_r = -  alpha[l] * da(i) * dffdt(l, k, j, 0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    } 

	    default: {
	    cout << "map_et_fait_d2rdtdx: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 

/*
 *******************************************************************************
 *		1/sin(th) d^2R/dphidx   ( dans la ZEC: - 1/sin(th) d^2U/dphidx )
 *******************************************************************************
 */

Mtbl* map_et_fait_sstd2rdpdx(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    const Valeur& stdffdp = ff.stdsdp() ; 
    const Valeur& stdggdp = gg.stdsdp() ; 
    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& da = *((cv->daa)[l]) ; 
	const Tbl& db = *((cv->dbb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case FIN: case RARE: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			*p_r = alpha[l] * (   da(i) * stdffdp(l, k, j, 0)
					    + db(i) * stdggdp(l, k, j, 0) ) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }
	        
	    case UNSURR: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			*p_r = -  alpha[l] * da(i) * stdffdp(l, k, j, 0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    } 

	    default: {
	    cout << "map_et_fait_sstd2rdpdx: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 

/*
 *******************************************************************************
 *	    1/R^2 d^2R/dtheta^2	( dans la ZEC: - 1/U^2 d^2U/dth^2 )
 *******************************************************************************
 */

Mtbl* map_et_fait_sr2d2rdt2(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = static_cast<const Map_et*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const double* alpha = cv->alpha ;
    const double* beta = cv->beta ;
    const Valeur& ff = cv->ff ; 
    const Valeur& gg = cv->gg ; 
    const Valeur& d2ffdt2 = ff.d2sdt2() ; 
    const Valeur& d2ggdt2 = gg.d2sdt2() ; 

    
    for (int l=0 ; l<nz ; l++) {

	const Grille3d* g = mg->get_grille3d(l) ;

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    
	    const Tbl& asx = cv->aasx ; 
	    const Tbl& bsx = cv->bbsx ; 
	    const Tbl& asx2 = cv->aasx2 ; 
	    const Tbl& bsx2 = cv->bbsx2 ; 
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			double ww =  1. + asx(i) * ff(l, k, j, 0)
				        + bsx(i) * gg(l, k, j, 0) ;
		    
			*p_r =	(   asx2(i) * d2ffdt2(l, k, j, 0) 
				  + bsx2(i) * d2ggdt2(l, k, j, 0)  ) /
				(alpha[l] * ww*ww) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }

	    case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {

			double ww = alpha[l] * (
					(g->x)[i] + aa(i) * ff(l, k, j, 0)
						  + bb(i) * gg(l, k, j, 0) ) 
				      + beta[l] ; 

			*p_r = alpha[l] * (   aa(i) * d2ffdt2(l, k, j, 0)
					    + bb(i) * d2ggdt2(l, k, j, 0) )
				/ ( ww*ww ) ; 
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

	    case UNSURR: {
	    
	    const Tbl& asxm1 = cv->zaasx ; 
	    const Tbl& asxm1car = cv->zaasx2 ; 

	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
		    
			double ww = 1. + asxm1(i) * ff(l, k, j, 0) ;
			
			*p_r = - asxm1car(i) * d2ffdt2(l, k, j, 0)
				/  ( alpha[l] * ww*ww )	;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }    

	    default: {
	    cout << "map_et_fait_sr2d2rdt2: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
} 

}
