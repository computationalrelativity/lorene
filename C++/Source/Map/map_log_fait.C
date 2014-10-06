/*
 *   Copyright (c) 2004 Philippe Granclement
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


char map_log_fait_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:13:13  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/06/22 08:49:58  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 *
 *
 * $Header$
 *
 */

// Includes
#include <cassert>
#include <cstdlib>
#include <cmath>

#include "mtbl.h"
#include "map.h"
#include "proto.h"

		    //----------------//
		    // Coord. radiale //
		    //----------------//

Mtbl* map_log_fait_r(const Map* cvi) {

    // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
    
    // Pour le confort
    Tbl alpha = cv->alpha ;
    Tbl beta = cv->beta ;
    Itbl type_var = cv->type_var ;

    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l) ;
	int ip = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch (type_var(l)) {
	case AFFINE : {

	  switch(mg->get_type_r(l)) {
	  case FIN: case RARE:
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = alpha(l) * (g->x)[i] + beta(l) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	  case UNSURR:
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 1./(alpha(l) * (g->x)[i] + beta(l)) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	  default:
	    cout << "Map_log_fait_r: unknown type_r !\n" ;
	    abort () ;
	    exit(-1) ;
	    
	  }	    // Fin du switch 1
	  break ;
	}

	case LOG : {
	  switch(mg->get_type_r(l)) {
	  case FIN:
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = exp(alpha(l) * (g->x)[i] + beta(l)) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	  default: {
	    cout << "Map_log_fait_r: unknown type_r !\n" ;
	    abort () ;
	    exit(-1) ;
	  }
	  }	    // Fin du switch 2
	  break ;
	}
	default: {
	  cout << "Map_log_fait_r: unknown type_r !\n" ;
	  abort () ;
	  exit(-1) ;
	}
	}
 
    }		// Fin de boucle sur zone
    
    // Termine
    return mti ;
}

		    //--------------//
		    // Coord. Theta //
		    //--------------//

Mtbl* map_log_fait_tet(const Map* cvi) {

    // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
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

Mtbl* map_log_fait_phi(const Map* cvi) {

    // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
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

Mtbl* map_log_fait_x(const Map* cvi) {

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

Mtbl* map_log_fait_y(const Map* cvi) {

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

Mtbl* map_log_fait_z(const Map* cvi) {

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

Mtbl* map_log_fait_xa(const Map* cvi) {

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

Mtbl* map_log_fait_ya(const Map* cvi) {

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

Mtbl* map_log_fait_za(const Map* cvi) {

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

Mtbl* map_log_fait_sint(const Map* cvi) {

    // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
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

Mtbl* map_log_fait_cost(const Map* cvi) {

    // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
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

Mtbl* map_log_fait_sinp(const Map* cvi) {

    // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
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

Mtbl* map_log_fait_cosp(const Map* cvi) {

    // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
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

Mtbl* map_log_fait_xsr(const Map* cvi) {

    // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
        
    // Pour le confort
    Tbl alpha = cv->alpha ;
    Tbl beta = cv->beta ;
    Itbl type_var = cv->type_var ;

    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l) ;
	int ip = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch (type_var(l)) {
	case AFFINE : {
	  
	  switch(mg->get_type_r(l)) {
	    
	  case RARE: {
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 1. / alpha(l) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	  }
	  case FIN: { 
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 1. / ( alpha(l) * (g->x)[i] + beta(l) ) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	  }
	  case UNSURR: {
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 1. / alpha(l) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	  }
	  default: {
	    cout << "map_log_fait_xsr: unknown type_r !" << endl ;
	    abort() ;
	  }
	  }	    // Fin du switch 1 
	  break ;
	}
	  
	case LOG: {
	  switch (mg->get_type_r(l)) {
	  case FIN: { 
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 1. / exp( alpha(l) * (g->x)[i] + beta(l) ) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	  }
	  default: {
	    cout << "map_log_fait_xsr: unknown type_r !" << endl ;
	    abort() ;
	  }
	  }
	  break ;
	}

	default:
	  cout << "map_log_fait_xsr: unknown type_r !" << endl ;
	  abort() ; 
	}			// Fin de boucle sur zone
    }
    // Termine
    return mti ;

}

/*
 ************************************************************************
 *			    1/(dR/dx)	    ( -1/(dU/dx) ds la ZEC )
 ************************************************************************
 */

Mtbl* map_log_fait_dxdr(const Map* cvi) {

 // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
        
    // Pour le confort
    Tbl alpha = cv->alpha ;
    Tbl beta = cv->beta ;
    Itbl type_var = cv->type_var ;

    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l) ;
	int ip = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch (type_var(l)) {
	case AFFINE : {
	  switch(mg->get_type_r(l)) {
	    
	  case RARE: case FIN:
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 1. / alpha(l) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    
	    
	  case UNSURR:
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = - 1. / alpha(l) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	  default:
	    cout << "map_log_fait_dxdr: unknown type_r !" << endl ;
	    abort() ;
	    
	  }	    // Fin du switch 1 
	  break ;
	}
	case LOG : {
	  switch(mg->get_type_r(l)) {
	  case FIN: 
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 1./ (alpha(l) * exp(alpha(l) * (g->x)[i] + beta(l))) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	  
	  default:
	    cout << "map_log_fait_dxdr: unknown type_r !" << endl ;
	    abort() ;
	  }
	  break ;
	}

	default:
	  cout << "map_log_fait_dxdr: unknown type_r !" << endl ;
	  abort() ; 
	}			// Fin de boucle sur zone
    }
    // Termine
    return mti ;

}
	  

  
/*
 ************************************************************************
 *			    dR/dtheta
 ************************************************************************
 */

Mtbl* map_log_fait_drdt(const Map* cvi) {

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

Mtbl* map_log_fait_stdrdp(const Map* cvi) {

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

Mtbl* map_log_fait_srdrdt(const Map* cvi) {

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

Mtbl* map_log_fait_srstdrdp(const Map* cvi) {

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

Mtbl* map_log_fait_sr2drdt(const Map* cvi) {

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

Mtbl* map_log_fait_sr2stdrdp(const Map* cvi) {

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

Mtbl* map_log_fait_d2rdx2(const Map* cvi) {
 // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
        
    // Pour le confort
    Tbl alpha = cv->alpha ;
    Tbl beta = cv->beta ;
    Itbl type_var = cv->type_var ;

    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l) ;
	int ip = mg->get_np(l) ;
	const Grille3d* g = mg->get_grille3d(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch (type_var(l)) {
	case AFFINE : {
	  switch(mg->get_type_r(l)) {
	    
	  case RARE: case FIN : case UNSURR:
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 0. ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	   
	  default:
	    cout << "map_log_fait_d2rdx2: unknown type_r !" << endl ;
	    abort() ;
	    
	  }	    // Fin du switch 1 
	  break ;
	}
	case LOG : {
	  switch(mg->get_type_r(l)) {
	  case FIN: 
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = exp (alpha(l) * (g->x)[i] + beta(l)) * 
		    alpha(l)*alpha(l) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	  
	  default:
	    cout << "map_log_fait_d2rdx2: unknown type_r !" << endl ;
	    abort() ;
	  }
	  break ;
	}
	default:
	  cout << "map_log_fait_d2rdx2: unknown type_r !" << endl ;
	  abort() ; 
	}			// Fin de boucle sur zone
    }
    // Termine
    return mti ;

} 

/*
 *****************************************************************************
 *  1/R^2 (  1/sin(th) d/dth( sin(th) dR/dth ) + 1/sin(th)^2 d^2R/dphi^2  )		    
 *****************************************************************************
 */

Mtbl* map_log_fait_lapr_tp(const Map* cvi) {

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

Mtbl* map_log_fait_d2rdtdx(const Map* cvi) {

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

Mtbl* map_log_fait_sstd2rdpdx(const Map* cvi) {

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

Mtbl* map_log_fait_sr2d2rdt2(const Map* cvi) {

    // recup de la grille 
    const Mg3d* mg = cvi->get_mg() ;
    
    // Le resultat est nul : 
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_zero() ; 

    return mti ; 
} 

/*
 ************************************************************************
 *			    1/(dR/dx)	    ( -1/(dU/dx) ds la ZEC )
 ************************************************************************
 */

Mtbl* map_log_fait_dxdlnr(const Map* cvi) {

 // recup du changement de variable
    const Map_log* cv = static_cast<const Map_log*>(cvi) ;
    const Mg3d* mg = cv->get_mg() ;
    int nz = mg->get_nzone() ;
    
    // Le resultat
    Mtbl* mti = new Mtbl(mg) ;
    mti->set_etat_qcq() ;
        
    // Pour le confort
    Tbl alpha = cv->alpha ;
    Tbl beta = cv->beta ;
    Itbl type_var = cv->type_var ;

    int i, j, k ;
    for (int l=0 ; l<nz ; l++) {
	int ir = mg->get_nr(l);
	int it = mg->get_nt(l) ;
	int ip = mg->get_np(l) ;
	Tbl* tb = (mti->t)[l] ;
	tb->set_etat_qcq() ;
	double* p_r = tb->t ;
	
	switch (type_var(l)) {
	case AFFINE : {
	  switch(mg->get_type_r(l)) {
	    
	  case RARE: case FIN:
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 1. / alpha(l) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ; 
	    
	    
	  case UNSURR:
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = - 1. / alpha(l) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	    
	  default:
	    cout << "map_log_fait_dxdr: unknown type_r !" << endl ;
	    abort() ;
	    
	  }	    // Fin du switch 1 
	  break ;
	}
	case LOG : {
	  switch(mg->get_type_r(l)) {
	  case FIN: 
	    for (k=0 ; k<ip ; k++) {
	      for (j=0 ; j<it ; j++) {
		for (i=0 ; i<ir ; i++) {
		  *p_r = 1./ alpha(l) ;
		  p_r++ ;
		}	    // Fin de boucle sur r
	      }	// Fin de boucle sur theta
	    }	    // Fin de boucle sur phi
	    break ;
	  
	  default:
	    cout << "map_log_fait_dxdr: unknown type_r !" << endl ;
	    abort() ;
	  }
	  break ;
	}

	default:
	  cout << "map_log_fait_dxdr: unknown type_r !" << endl ;
	  abort() ; 
	}			// Fin de boucle sur zone
    }
    // Termine
    return mti ;

}
