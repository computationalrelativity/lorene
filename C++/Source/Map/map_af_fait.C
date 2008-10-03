/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
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


char map_af_fait_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2008/10/03 09:05:29  j_novak
 * Improved the treatment of angular mapping in the computation of xsr
 *
 * Revision 1.7  2007/12/20 09:11:04  jl_cornou
 * Correction of an error in op_sxpun about Jacobi(0,2) polynomials
 *
 * Revision 1.6  2007/12/14 10:19:30  jl_cornou
 * *** empty log message ***
 *
 * Revision 1.5  2007/12/11 15:28:14  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.4  2006/06/09 14:58:48  j_novak
 * Added a hack in the case of pure angular grid for the Coord xsr in the shells.
 *
 * Revision 1.3  2003/10/15 10:34:46  e_gourgoulhon
 * Added new Coord's: drdt and stdrdp.
 *
 * Revision 1.2  2002/10/16 14:36:41  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  1999/10/15  09:16:38  eric
 * Changement prototypes: const.
 *
 * Revision 2.5  1999/10/14  14:27:26  eric
 * const.
 *
 * Revision 2.4  1999/09/30  12:56:10  eric
 * const Grille3d*
 *
 * Revision 2.3  1999/04/09  12:59:21  phil
 * correction de map_af_fait_dxdr
 *
 * Revision 2.2  1999/03/04  15:22:53  eric
 * *** empty log message ***
 *
 * Revision 2.1  1999/03/04  13:12:03  eric
 * Ajout des derivees du changement de variable
 *
 * Revision 2.0  1999/02/15  10:42:45  hyc
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/15  09:59:50  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/01/15  09:10:39  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Includes
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "mtbl.h"
#include "map.h"
#include "proto.h"

		    //----------------//
		    // Coord. radiale //
		    //----------------//

Mtbl* map_af_fait_r(const Map* cvi) {

    // recup du changement de variable
    const Map_af* cv = static_cast<const Map_af*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    double* alpha = cv->alpha ;
    double* beta = cv->beta ;
    
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
	    case FIN: case RARE: case FINJAC :
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = alpha[l] * (g->x)[i] + beta[l] ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	    case UNSURR:
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = 1./(alpha[l] * (g->x)[i] + beta[l]) ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	    default:
	    cout << "Map_af_fait_r: unknown type_r !\n" ;
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

Mtbl* map_af_fait_tet(const Map* cvi) {

    // recup du changement de variable
    const Map_af* cv = static_cast<const Map_af*>(cvi) ;
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

Mtbl* map_af_fait_phi(const Map* cvi) {

    // recup du changement de variable
    const Map_af* cv = static_cast<const Map_af*>(cvi) ;
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

Mtbl* map_af_fait_x(const Map* cvi) {

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

Mtbl* map_af_fait_y(const Map* cvi) {

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

Mtbl* map_af_fait_z(const Map* cvi) {

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

Mtbl* map_af_fait_xa(const Map* cvi) {

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

Mtbl* map_af_fait_ya(const Map* cvi) {

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

Mtbl* map_af_fait_za(const Map* cvi) {

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

Mtbl* map_af_fait_sint(const Map* cvi) {

    // recup du changement de variable
    const Map_af* cv = static_cast<const Map_af*>(cvi) ;
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

Mtbl* map_af_fait_cost(const Map* cvi) {

    // recup du changement de variable
    const Map_af* cv = static_cast<const Map_af*>(cvi) ;
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

Mtbl* map_af_fait_sinp(const Map* cvi) {

    // recup du changement de variable
    const Map_af* cv = static_cast<const Map_af*>(cvi) ;
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

Mtbl* map_af_fait_cosp(const Map* cvi) {

    // recup du changement de variable
    const Map_af* cv = static_cast<const Map_af*>(cvi) ;
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

Mtbl* map_af_fait_xsr(const Map* cvi) {

    // recup du changement de variable
    const Map_af* cv = static_cast<const Map_af*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
        
    // Pour le confort
    double* alpha = cv->alpha ;
    double* beta = cv->beta ;
    
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
	    assert(beta[l]==0) ;
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = 1. / alpha[l] ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 

	    case FINJAC:
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = 1. / alpha[l] ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	
	    case FIN: 
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    if (ir == 1) { //Some hack for angular grid case...
			*p_r = 1. / beta[l] ;
			p_r++ ;
		    }
		    else 
			for (i=0 ; i<ir ; i++) {
			    *p_r = 1. / ( alpha[l] * (g->x)[i] + beta[l] ) ;
			    p_r++ ;
			}	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	    case UNSURR:
	    assert(beta[l] == - alpha[l]) ;
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = 1. / alpha[l] ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	    default:
	    cout << "map_af_fait_xsr: unknown type_r !" << endl ;
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

Mtbl* map_af_fait_dxdr(const Map* cvi) {

    // recup du changement de variable
    const Map_af* cv = static_cast<const Map_af*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    double* alpha = cv->alpha ;
    
    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l) ;
	int ip = mg->get_np(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch(mg->get_type_r(l)) {
	
	    case FIN: case RARE: case FINJAC :
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = 1. / alpha[l] ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
		    
	    case UNSURR:
	    for (k=0 ; k<ip ; k++) {
		for (j=0 ; j<it ; j++) {
		    for (i=0 ; i<ir ; i++) {
			*p_r = - 1. / alpha[l] ;
			p_r++ ;
		    }	    // Fin de boucle sur r
		}	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
		    

	    default:
	    cout << "map_af_fait_dxdr: unknown type_r !" << endl ;
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

Mtbl* map_af_fait_drdt(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/sin(theta) dR/dphi
 ************************************************************************
 */

Mtbl* map_af_fait_stdrdp(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/R dR/dtheta
 ************************************************************************
 */

Mtbl* map_af_fait_srdrdt(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/(R sin(theta)) dR/dphi
 ************************************************************************
 */

Mtbl* map_af_fait_srstdrdp(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/R^2 dR/dtheta
 ************************************************************************
 */

Mtbl* map_af_fait_sr2drdt(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/(R^2 sin(theta)) dR/dphi
 ************************************************************************
 */

Mtbl* map_af_fait_sr2stdrdp(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    d^2R/dx^2
 ************************************************************************
 */

Mtbl* map_af_fait_d2rdx2(const Map* cvi) {

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

Mtbl* map_af_fait_lapr_tp(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    d^2R/dthdx
 ************************************************************************
 */

Mtbl* map_af_fait_d2rdtdx(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/sin(th) d^2R/dphidx
 ************************************************************************
 */

Mtbl* map_af_fait_sstd2rdpdx(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    d^2R/dtheta^2
 ************************************************************************
 */

Mtbl* map_af_fait_sr2d2rdt2(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

