/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char map_et_fait_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.1  1999/11/24  11:23:00  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Includes
#include <assert.h>
#include <stdlib.h>
#include <iostream.h>
#include <math.h>

#include "map.h"

		    //----------------//
		    // Coord. radiale //
		    //----------------//

Mtbl* map_et_fait_r(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = (const Map_et*) cvi ;
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

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	int np = mg->get_np(l) ;
	int nt = mg->get_nt(l) ;
	int nr = mg->get_nr(l) ;
	
	switch(mg->get_type_r(l)) {

	    case FIN: case RARE: {

	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = alpha[l] * ( (g->x)[i] 
					    + aa(i) * ff(l, k, j, 0)
					    + bb(i) * gg(l, k, j, 0)
					   ) + beta[l] ;
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
			*p_r = 1./( alpha[l]  * (
					    (g->x)[i] + aa(i) * ff(l, k, j, 0)
						)    
				    + beta[l] ) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }
	    
	    default: {
	    cout << "map_et_fait_r: Unknown type_r !" << endl ;
	    abort () ;
	    }
	    
	}	    // Fin du switch
    }			// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

		    //--------------//
		    // Coord. Theta //
		    //--------------//

Mtbl* map_et_fait_tet(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = (const Map_et *) cvi ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
        
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    for (int l=0 ; l<nz ; l++) {
	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l);
	int np = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (int k=0 ; k<np ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		for (int i=0 ; i<nr ; i++) {
		    *p_r = (g->tet)[j] ;
		    p_r++ ;
		}   // Fin de boucle sur r
	    }	// Fin de boucle sur theta
	}   // Fin de boucle sur phi
    }	// Fin de boucle sur zone

    // Termine
    return mti ;    
}

			//------------//
			// Coord. Phi //
			//------------//

Mtbl* map_et_fait_phi(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = (const Map_et *) cvi ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    for (int l=0 ; l<nz ; l++) {
	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l);
	int np = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (int k=0 ; k<np ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		for (int i=0 ; i<nr ; i++) {
		    *p_r = (g->phi)[k] ;
		    p_r++ ;
		}   // Fin de boucle sur r
	    }	// Fin de boucle sur theta
	}   // Fin de boucle sur phi
    }	// Fin de boucle sur zone
   
    // Termine
    return mti ; 
}

			//----------//
			// Coord. X //
			//----------//

Mtbl* map_et_fait_x(const Map* cvi) {

    // recup de la grille
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    
    *mti = (cvi->r) * (cvi->sint) * (cvi->cosp) ;

    // Termine
    return mti ;
}

			//----------//
			// Coord. Y //
			//----------//

Mtbl* map_et_fait_y(const Map* cvi) {

    // recup de la grille
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    
    *mti = (cvi->r) * (cvi->sint) * (cvi->sinp) ;
   
    // Termine
    return mti ; 
}

			//----------//
			// Coord. Z //
			//----------//

Mtbl* map_et_fait_z(const Map* cvi) {

    // recup de la grille
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    
    *mti = (cvi->r) * (cvi->cost) ;
    
    // Termine
    return mti ;
}

			//--------------------//
			// Coord. X "absolue" //
			//--------------------//

Mtbl* map_et_fait_xa(const Map* cvi) {

    // recup de la grille
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    
    double r_phi = cvi->get_rot_phi() ; 
    double t_x = cvi->get_ori_x() ; 

    if ( fabs(r_phi) < 1.e-14 ) {	
	*mti = (cvi->x) + t_x ; 
    }
    else if ( fabs(r_phi - M_PI) < 1.e-14 ) {
	*mti = - (cvi->x) + t_x ; 	
    }
    else {
	Mtbl phi = cvi->phi + r_phi ;
	*mti = (cvi->r) * (cvi->sint) * cos(phi) + t_x ;
    }

    // Termine
    return mti ;
}

			//--------------------//
			// Coord. Y "absolue" //
			//--------------------//

Mtbl* map_et_fait_ya(const Map* cvi) {

    // recup de la grille
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    
    double r_phi = cvi->get_rot_phi() ; 
    double t_y = cvi->get_ori_y() ; 

    if ( fabs(r_phi) < 1.e-14 ) {	
	*mti = (cvi->y) + t_y ; 
    }
    else if ( fabs(r_phi - M_PI) < 1.e-14 ) {
	*mti = - (cvi->y) + t_y ; 	
    }
    else {
	Mtbl phi = cvi->phi + r_phi ;
	*mti = (cvi->r) * (cvi->sint) * sin(phi) + t_y ;
    }

    // Termine
    return mti ;
}

			//--------------------//
			// Coord. Z "absolue" //
			//--------------------//

Mtbl* map_et_fait_za(const Map* cvi) {

    // recup de la grille
    const Mg3d* mg = cvi->get_mg() ;
    
    double t_z = cvi->get_ori_z() ; 

    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    
    *mti = (cvi->r) * (cvi->cost) + t_z ; 

    // Termine
    return mti ;       
}

			//---------------//
			// Trigonometrie //
			//---------------//

Mtbl* map_et_fait_sint(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = (const Map_et *) cvi ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    for (int l=0 ; l<nz ; l++) {
	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l);
	int np = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (int k=0 ; k<np ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		for (int i=0 ; i<nr ; i++) {
		    *p_r = sin(g->tet[j]) ;
		    p_r++ ;
		}   // Fin de boucle sur r
	    }	// Fin de boucle sur theta
	}   // Fin de boucle sur phi
    }	// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

Mtbl* map_et_fait_cost(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = (const Map_et *) cvi ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    for (int l=0 ; l<nz ; l++) {
	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l);
	int np = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (int k=0 ; k<np ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		for (int i=0 ; i<nr ; i++) {
		    *p_r = cos(g->tet[j]) ;
		    p_r++ ;
		}   // Fin de boucle sur r
	    }	// Fin de boucle sur theta
	}   // Fin de boucle sur phi
    }	// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

Mtbl* map_et_fait_sinp(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = (const Map_et *) cvi ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    for (int l=0 ; l<nz ; l++) {
	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l);
	int np = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (int k=0 ; k<np ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		for (int i=0 ; i<nr ; i++) {
		    *p_r = sin(g->phi[k]) ;
		    p_r++ ;
		}   // Fin de boucle sur r
	    }	// Fin de boucle sur theta
	}   // Fin de boucle sur phi
    }	// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

Mtbl* map_et_fait_cosp(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = (const Map_et *) cvi ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    for (int l=0 ; l<nz ; l++) {
	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l);
	int np = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (int k=0 ; k<np ; k++) {
	    for (int j=0 ; j<nt ; j++) {
		for (int i=0 ; i<nr ; i++) {
		    *p_r = cos(g->phi[k]) ;
		    p_r++ ;
		}   // Fin de boucle sur r
	    }	// Fin de boucle sur theta
	}   // Fin de boucle sur phi
    }	// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

/*
 ************************************************************************
 *	x/R dans le noyau,  1/R dans les coquilles,  (x-1)/U dans la ZEC
 ************************************************************************
 */

Mtbl* map_et_fait_xsr(const Map* cvi) {

    // recup du changement de variable
    const Map_et* cv = (const Map_et *) cvi ;
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
    const Tbl& asx = cv->aasx ; 
    const Tbl& bsx = cv->bbsx ; 
    const Tbl& asxm1 = cv->zaasx ; 
    
    for (int l=0 ; l<nz ; l++) {
	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;

	const Tbl& aa = *((cv->aa)[l]) ; 
	const Tbl& bb = *((cv->bb)[l]) ; 

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    assert(beta[l]==0) ;
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = 1. / ( alpha[l] * ( 1. + asx(i) * ff(l, k, j, 0)
						      + bsx(i) * gg(l, k, j, 0)
						  ) )	;
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
			*p_r = 1. / ( alpha[l] * ( (g->x)[i] 
						    + aa(i) * ff(l, k, j, 0)
						    + bb(i) * gg(l, k, j, 0)
						) + beta[l] );
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }
	    
	    case UNSURR: {
	    assert(beta[l] == - alpha[l]) ;
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = 1. / ( alpha[l] * ( 1. 
						    + asxm1(i) * ff(l, k, j, 0)
						  ) )	;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }
	    
	    default: {
	    cout << "map_et_fait_xsr: unknown type_r !" << endl ;
	    abort() ;
	    }
	    
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;

}

