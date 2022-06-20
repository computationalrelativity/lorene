/*
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


// Includes
#include <cassert>
#include <cstdlib>
#include <cmath>

#include "mtbl.h"
#include "map.h"
#include "proto.h"
#include "valeur.h"








			//----------------//
		    // Coord. radiale //
		    //----------------//

namespace Lorene {
Mtbl* map_star_fait_r(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& beta = cv->get_beta() ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l) ;
	int ip = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	    case RARE: case FIN:
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = alpha(l,k,j,0)  * (g->x)[i] + beta(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;

		// case FIN:{
		// cout << "map_star_fait_r: Shells not implemented yet..." << endl;
		// abort() ; 
		// break ;
		// }
	    
	    case UNSURR:
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			cout << "map_star_fait_r: Warning ! No compactified zone allowed !" << endl;
			abort() ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	    default:
	    cout << "map_star_fait_r: unknown type_r !\n" ;
	    abort () ;
	    exit(-1) ;
	    
	}	    // Fin du switch
    }			// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

		    //--------------//
		    // Coord. Theta //
		    //--------------//

Mtbl* map_star_fait_tet(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
        
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l);
	int ip = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (k=0 ; k<ip ; k++) {
	    for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
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

Mtbl* map_star_fait_phi(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l);
	int ip = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (k=0 ; k<ip ; k++) {
	    for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
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

Mtbl* map_star_fait_x(const Map* cvi) {

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

Mtbl* map_star_fait_y(const Map* cvi) {

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

Mtbl* map_star_fait_z(const Map* cvi) {

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

Mtbl* map_star_fait_xa(const Map* cvi) {

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

Mtbl* map_star_fait_ya(const Map* cvi) {

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

Mtbl* map_star_fait_za(const Map* cvi) {

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

Mtbl* map_star_fait_sint(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l);
	int ip = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (k=0 ; k<ip ; k++) {
	    for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		    *p_r = sin(g->tet[j]) ;
		    p_r++ ;
		}   // Fin de boucle sur r
	    }	// Fin de boucle sur theta
	}   // Fin de boucle sur phi
    }	// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

Mtbl* map_star_fait_cost(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l);
	int ip = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (k=0 ; k<ip ; k++) {
	    for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		    *p_r = cos(g->tet[j]) ;
		    p_r++ ;
		}   // Fin de boucle sur r
	    }	// Fin de boucle sur theta
	}   // Fin de boucle sur phi
    }	// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

Mtbl* map_star_fait_sinp(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l);
	int ip = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (k=0 ; k<ip ; k++) {
	    for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		    *p_r = sin(g->phi[k]) ;
		    p_r++ ;
		}   // Fin de boucle sur r
	    }	// Fin de boucle sur theta
	}   // Fin de boucle sur phi
    }	// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

Mtbl* map_star_fait_cosp(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l);
	int ip = mg->get_np(l);
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	for (k=0 ; k<ip ; k++) {
	    for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
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

Mtbl* map_star_fait_xsr(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
        
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& beta = cv->get_beta() ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l) ;
	int ip = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE:
		assert(beta(l).get_etat()==ETATZERO) ;
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = 1. / alpha(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 

	    case FIN: 
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    if (ir == 1) { //Some hack for angular grid case...
			*p_r = 1. / beta(l,k,j,0) ;
			p_r++ ;
		    }
		    else 
			for (i=0 ; i<ir ; i++) {
			    *p_r = 1. / ( alpha(l,k,j,0) * (g->x)[i] + beta(l,k,j,0) ) ;
			    p_r++ ;
			}	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	    case UNSURR:
		cout << "map_star_fait_xsr: Compactified domain not allowed !" << endl;
		abort() ;
	    break ;
	    
	    default:
	    cout << "map_star_fait_xsr: unknown type_r !" << endl ;
	    abort() ;
	    
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;

}

/*
 ************************************************************************
 *			    1/(dR/dx)	    ( -1/(dU/dx) ds la ZEC )
 ************************************************************************
 */

Mtbl* map_star_fait_dxdr(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l) ;
	int ip = mg->get_np(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: case FIN:
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = 1. / alpha(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 

		// case FIN:{
		// cout << "map_star_fait_dxdr: Shells not implemented yet..." << endl;
		// abort() ; 
		// break ;
		// }
		    
	    case UNSURR:
	    cout << "map_star_fait_dxdr: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
		    

	    default:
	    cout << "map_star_fait_dxdr: unknown type_r !" << endl ;
	    abort() ;
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine
    return mti ;
}

/*
 ************************************************************************
 *			    dR/dtheta
 ************************************************************************
 */

Mtbl* map_star_fait_drdt(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& dalphadt = alpha.dsdt() ;
	const Valeur& beta = cv->get_beta() ;
	const Valeur& dbetadt = beta.dsdt() ;

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE: case FIN:  {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = dalphadt(l,k,j,0) * (g->x)[i] + dbetadt(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }
		// case FIN:{
		// cout << "map_star_fait_drdt: Shells not implemented yet..." << endl;
		// abort() ; 
		// break ;
		// }

	    case UNSURR: {
	    cout << "map_star_fait_drdt: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_drdt: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/sin(theta) dR/dphi
 ************************************************************************
 */

Mtbl* map_star_fait_stdrdp(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& stalphadp = alpha.stdsdp() ;
	const Valeur& beta = cv->get_beta() ;
	const Valeur& stbetadp = beta.stdsdp() ;

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE: case FIN: {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = stalphadp(l,k,j,0) * (g->x)[i] + stbetadp(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

		// case FIN:{
		// cout << "map_star_fait_stdrdp: Shells not implemented yet..." << endl;
		// abort() ; 
		// break ;
		// }

	    case UNSURR: {
	    cout << "map_star_fait_stdrdp: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_stdrdp: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/R dR/dtheta
 ************************************************************************
 */

Mtbl* map_star_fait_srdrdt(const Map* cvi) {

  	// recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& dalphadt = alpha.dsdt() ;
	const Valeur& beta = cv->get_beta() ;
	const Valeur& dbetadt = beta.dsdt() ;

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE : {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = dalphadt(l,k,j,0) / alpha(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

		case FIN:{
		for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = (dalphadt(l,k,j,0)*(g->x)[i] + dbetadt(l,k,j,0)) / (alpha(l,k,j,0)*(g->x)[i] + beta(l,k,j,0)) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
		break ;
		}

	    case UNSURR: {
	    cout << "map_star_fait_srdrdt: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_srdrdt: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/(R sin(theta)) dR/dphi
 ************************************************************************
 */

Mtbl* map_star_fait_srstdrdp(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& stalphadp = alpha.stdsdp() ;
	const Valeur& beta = cv->get_beta() ;
	const Valeur& stbetadp = beta.stdsdp() ;

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE : {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = stalphadp(l,k,j,0) / alpha(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

		case FIN:{
		for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = (stalphadp(l,k,j,0)*(g->x)[i] + stbetadp(l,k,j,0)) / (alpha(l,k,j,0)*(g->x)[i] + beta(l,k,j,0)) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
		break ;
		}

	    case UNSURR: {
	    cout << "map_star_fait_srstdrdp: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_srstdrdp: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/R^2 dR/dtheta
 ************************************************************************
 */

Mtbl* map_star_fait_sr2drdt(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& dalphadt = alpha.dsdt() ;

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE : {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			double ww = alpha(l,k,j,0) ;
			*p_r = dalphadt(l,k,j,0) / (ww*ww * (g->x)[i])  ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }
		
		case FIN:{
		cout << "map_star_fait_sr2drdt: Shells not implemented yet..." << endl;
		abort() ; 
		break ;
		}

	    case UNSURR: {
	    cout << "map_star_fait_sr2drdt: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_sr2drdt: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/(R^2 sin(theta)) dR/dphi
 ************************************************************************
 */

Mtbl* map_star_fait_sr2stdrdp(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& stalphadp = alpha.stdsdp() ;

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE : {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			double ww = alpha(l,k,j,0) ;
			*p_r = stalphadp(l,k,j,0) / (ww*ww * (g->x)[i]) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

		case FIN:{
		cout << "map_star_fait_sr2stdrdp: Shells not implemented yet..." << endl;
		abort() ; 
		break ;
		}

	    case UNSURR: {
	    cout << "map_star_fait_sr2stdrdp: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_sr2stdrdp: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine

    return mti ; 
} 

/*
 ************************************************************************
 *			    d^2R/dx^2
 ************************************************************************
 */

Mtbl* map_star_fait_d2rdx2(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 *****************************************************************************
 *  1/R^2 (  1/sin(th) d/dth( sin(th) dR/dth ) + 1/sin(th)^2 d^2R/dphi^2  )		    
 *****************************************************************************
 */

Mtbl* map_star_fait_lapr_tp(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;

	Valeur alpha_tmp = alpha ;
	alpha_tmp.ylm() ;
    const Valeur& lapalpha = alpha_tmp.lapang() ; 

    
    for (int l=0 ; l<nz ; l++) {

	const Grille3d* g = mg->get_grille3d(l) ;

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;

	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case RARE: {
	    
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = (g->x)[i]*lapalpha(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    }

		case FIN:{
		cout << "map_star_fait_lapr_tp: Shells not implemented yet..." << endl;
		abort() ; 
		break ;
		}

	    case UNSURR: {
	    cout << "map_star_fait_lapr_tp: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_lapr_tp: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine

    return mti ; 
} 

/*
 ************************************************************************
 *			    d^2R/dthdx
 ************************************************************************
 */

Mtbl* map_star_fait_d2rdtdx(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& d2alphadxdt = alpha.dsdt().dsdx() ;

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE :  {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = d2alphadxdt(l,k,j,0)/alpha(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }
		case FIN:{
		cout << "map_star_fait_d2rdtdx: Shells not implemented yet..." << endl;
		abort() ; 
		break ;
		}

	    case UNSURR: {
	    cout << "map_star_fait_d2rdtdx: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_d2rdtdx: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/sin(th) d^2R/dphidx
 ************************************************************************
 */

Mtbl* map_star_fait_sstd2rdpdx(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& stalphadp = alpha.stdsdp() ;

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE : {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			*p_r = stalphadp(l,k,j,0) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }

		case FIN:{
		cout << "map_star_fait_sstd2rdpdx: Shells not implemented yet..." << endl;
		abort() ; 
		break ;
		}

	    case UNSURR: {
	    cout << "map_star_fait_sstd2rdpdx: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_sstd2rdpdx: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine


    return mti ; 
} 

/*
 ************************************************************************
 *			 1/R^2 d^2R/dtheta^2
 ************************************************************************
 */

Mtbl* map_star_fait_sr2d2rdt2(const Map* cvi) {

    // recup du changement de variable
    const Map_star* cv = static_cast<const Map_star*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    const Valeur& alpha = cv->get_alpha() ;
	const Valeur& d2alphadt2 = alpha.d2sdt2() ;

    
    for (int l=0 ; l<nz ; l++) {

	int nr = mg->get_nr(l);
	int nt = mg->get_nt(l) ;
	int np = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
		    
	    case RARE :  {
	    for (int k=0 ; k<np ; k++) {
		for (int j=0 ; j<nt ; j++) {
		    for (int i=0 ; i<nr ; i++) {
			double ww = alpha(l,k,j,0) ;
			*p_r = d2alphadt2(l,k,j,0) / (ww*ww * (g->x)[i]) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    }
		case FIN:{
		cout << "map_star_fait_drdt: Shells not implemented yet..." << endl;
		abort() ; 
		break ;
		}

	    case UNSURR: {
	    cout << "map_star_fait_drdt: Compactified domain not allowed !" << endl;
		abort() ;
	    break ; 
	    }    

	    default: {
	    cout << "map_star_fait_drdt: unknown type_r !" << endl ;
	    abort() ;
	    }
	}	    // Fin du switch
    }			// Fin de boucle sur zone

    // Termine

    return mti ; 
} 

}
