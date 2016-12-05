/*
 *  Resolution of the angular Poisson equation. 
 *
 *
 */

/*
 *   Copyright (c) 2003-2005 Eric Gourgoulhon & Jerome Novak
 *
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


 


/*
 * $Id$
 * $Log$
 * Revision 1.10  2016/12/05 16:18:09  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.9  2014/10/13 08:53:29  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:16:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2009/10/23 12:55:04  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.6  2009/10/13 19:45:01  j_novak
 * New base T_LEG_MP.
 *
 * Revision 1.5  2005/04/08 07:36:20  f_limousin
 * Add #include <math.h> to avoid error in the compilation with gcc 3.3.1
 * (problem with fabs).
 *
 * Revision 1.4  2005/04/05 08:34:10  e_gourgoulhon
 * Treatment case l(l+1) = lambda.
 *
 * Revision 1.3  2005/04/04 21:33:37  e_gourgoulhon
 * Solving now for the generalized angular Poisson equation
 *    Lap_ang u + lambda u = source
 * with new parameter lambda.
 *
 * Revision 1.2  2004/12/17 13:35:03  m_forot
 * Add the case T_LEG
 *
 * Revision 1.1  2003/10/15 21:13:28  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "mtbl_cf.h"
#include "grilles.h"
#include "type_parite.h"


		//--------------------------------------
		// Routine pour les cas non prevus ----
		//------------------------------------

namespace Lorene {
void _poisangu_pas_prevu(Mtbl_cf* mt, int l, double) {
    cout << "Unknwon theta basis in the operator Mtbl_cf::poisson_angu() !" << endl ;
    cout << " basis : " << hex << (mt->base).b[l] << endl ; 
    abort () ;
}

			//---------------
			// cas T_LEG --
			//---------------

void _poisangu_t_leg(Mtbl_cf* mt, int l, double lambda)
{

    Tbl* tb = mt->t[l] ;	    // pt. sur tbl de travail
    
    // Peut-etre rien a faire ?
    	if (tb->get_etat() == ETATZERO) {
		return ;
    }
    
    int k, j, i ; 
    // Pour le confort
    int nr = mt->get_mg()->get_nr(l) ;   // Nombre
    int nt = mt->get_mg()->get_nt(l) ;   //	de points
    int np = mt->get_mg()->get_np(l) ;   //	    physiques
    
    int np1 = ( np == 1 ) ? 1 : np+1 ; 
	
    double* tuu = tb->t ; 

    // k = 0  :
     
    for (j=0 ; j<nt ; j++) {
		int ll = j ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
		tuu  += nr ;
    }     // Fin de boucle sur theta

    // On saute k = 1 : 
    tuu += nt*nr ; 
	
    // k=2,...
    for (k=2 ; k<np1 ; k++) {
	int m = k/2 ;
	tuu  += (m/2)*nr ;
	for (j=m/2 ; j<nt ; j++) {
	    int ll = j ;
	    double xl = - ll*(ll+1) + lambda ;
		
		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	    tuu  += nr ;
	}     // Fin de boucle sur theta
    }	// Fin de boucle sur phi	

    // base de developpement inchangee 
}

			//---------------
			// cas T_LEG_P --
			//---------------

void _poisangu_t_leg_p(Mtbl_cf* mt, int l, double lambda)
{

    Tbl* tb = mt->t[l] ;	    // pt. sur tbl de travail
    
    // Peut-etre rien a faire ?
    	if (tb->get_etat() == ETATZERO) {
		return ;
    }
    
    int k, j, i ; 
    // Pour le confort
    int nr = mt->get_mg()->get_nr(l) ;   // Nombre
    int nt = mt->get_mg()->get_nt(l) ;   //	de points
    int np = mt->get_mg()->get_np(l) ;   //	    physiques
    
    int np1 = ( np == 1 ) ? 1 : np+1 ; 
	
    double* tuu = tb->t ; 

    // k = 0  :
     
    for (j=0 ; j<nt ; j++) {
		int ll = 2*j ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
		tuu  += nr ;
    }     // Fin de boucle sur theta

    // On saute k = 1 : 
    tuu += nt*nr ; 
	
    // k=2,...
    for (k=2 ; k<np1 ; k++) {
	int m = k/2 ;
	tuu  += (m/2)*nr ;
	for (j=m/2 ; j<nt ; j++) {
	    int ll = 2*j + (m%2) ;
	    double xl = - ll*(ll+1) + lambda ;
		
		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	    tuu  += nr ;
	}     // Fin de boucle sur theta
    }	// Fin de boucle sur phi	
	    
    // base de developpement inchangee 
}

			//------------------
			// cas T_LEG_PP --
			//----------------

void _poisangu_t_leg_pp(Mtbl_cf* mt, int l, double lambda)
{

    Tbl* tb = mt->t[l] ;	    // pt. sur tbl de travail
    
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    int k, j, i ; 
    // Pour le confort
    int nr = mt->get_mg()->get_nr(l) ;   // Nombre
    int nt = mt->get_mg()->get_nt(l) ;   //	de points
    int np = mt->get_mg()->get_np(l) ;   //	    physiques
    
    int np1 = ( np == 1 ) ? 1 : np+1 ; 
	
    double* tuu = tb->t ; 

    // k = 0  :
     
    for (j=0 ; j<nt ; j++) {
		int ll = 2*j ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	tuu  += nr ;
    }     // Fin de boucle sur theta

    // On saute k = 1 : 
    tuu += nt*nr ; 
	
    // k=2,...
    for (k=2 ; k<np1 ; k++) {
	int m = 2*(k/2);
	tuu  += (m/2)*nr ;
	for (j=m/2 ; j<nt ; j++) {
	    int ll = 2*j ;
	    double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	    tuu  += nr ;
	}     // Fin de boucle sur theta
    }	// Fin de boucle sur phi	
	    
    // base de developpement inchangee 
}

			//----------------
			// cas T_LEG_I --
			//---------------

void _poisangu_t_leg_i(Mtbl_cf* mt, int l, double lambda)
{

    Tbl* tb = mt->t[l] ;	    // pt. sur tbl de travail
    
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    int k, j, i ; 
    // Pour le confort
    int nr = mt->get_mg()->get_nr(l) ;   // Nombre
    int nt = mt->get_mg()->get_nt(l) ;   //	de points
    int np = mt->get_mg()->get_np(l) ;   //	    physiques
    
    int np1 = ( np == 1 ) ? 1 : np+1 ; 
	
    double* tuu = tb->t ; 

    // k = 0  :
     
    for (j=0 ; j<nt-1 ; j++) {
		int ll = 2*j+1 ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	tuu  += nr ;
    }     // Fin de boucle sur theta
    tuu  += nr ; // On saute j=nt-1

    // On saute k = 1 : 
    tuu += nt*nr ; 
	
    // k=2,...
    for (k=2 ; k<np1 ; k++) {
	int m = k/2 ;
	tuu  += ((m+1)/2)*nr ;
	for (j=(m+1)/2 ; j<nt-1 ; j++) {
	    int ll = 2*j + ((m+1)%2) ;
	    double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	    tuu  += nr ;
	}     // Fin de boucle sur theta
	tuu  += nr ; // On saute j=nt-1
    }	// Fin de boucle sur phi	
	    
    // base de developpement inchangee 
}

			//------------------
			// cas T_LEG_IP --
			//----------------

void _poisangu_t_leg_ip(Mtbl_cf* mt, int l, double lambda)
{

    Tbl* tb = mt->t[l] ;	    // pt. sur tbl de travail
    
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    int k, j, i ; 
    // Pour le confort
    int nr = mt->get_mg()->get_nr(l) ;   // Nombre
    int nt = mt->get_mg()->get_nt(l) ;   //	de points
    int np = mt->get_mg()->get_np(l) ;   //	    physiques
    
    int np1 = ( np == 1 ) ? 1 : np+1 ; 
	
    double* tuu = tb->t ; 

    // k = 0  :
     
    for (j=0 ; j<nt-1 ; j++) {
		int ll = 2*j+1 ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
		tuu  += nr ;
    }     // Fin de boucle sur theta
    tuu  += nr ; // On saute j=nt-1

    // On saute k = 1 : 
    tuu += nt*nr ; 
	
    // k=2,...
    for (k=2 ; k<np1 ; k++) {
	int m = 2*(k/2);
	tuu  += (m/2)*nr ;
	for (j=m/2 ; j<nt-1 ; j++) {
	    int ll = 2*j+1 ;
	    double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	    tuu  += nr ;
	}     // Fin de boucle sur theta
	tuu  += nr ; // On saute j=nt-1
    }	// Fin de boucle sur phi	

//## Verif
    assert (tuu == tb->t + (np+1)*nt*nr) ;
	    
    // base de developpement inchangee 
}

			//------------------
			// cas T_LEG_PI --
			//----------------

void _poisangu_t_leg_pi(Mtbl_cf* mt, int l, double lambda)
{

    Tbl* tb = mt->t[l] ;	    // pt. sur tbl de travail
    
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    int k, j, i ; 
    // Pour le confort
    int nr = mt->get_mg()->get_nr(l) ;   // Nombre
    int nt = mt->get_mg()->get_nt(l) ;   //	de points
    int np = mt->get_mg()->get_np(l) ;   //	    physiques
    
    double* tuu = tb->t ; 

    // k = 0  :	    cos(phi)
    // -----
     
    for (j=0 ; j<nt-1 ; j++) {
		int ll = 2*j+1 ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
		tuu  += nr ;
    }     // Fin de boucle sur theta
    tuu  += nr ; // On saute j=nt-1

    if (np==1) {
	return ; 
    }

    // k = 1 : on saute
    // -----
    tuu += nt*nr ; 
	
    // k = 2 :	sin(phi)
    // ------
    for (j=0 ; j<nt-1 ; j++) {
		int ll = 2*j+1 ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
		tuu  += nr ;
    }     // Fin de boucle sur theta
    tuu  += nr ; // On saute j=nt-1

    // 3 <= k <= np
    // ------------
    for (k=3 ; k<np+1 ; k++) {
	int m = (k%2 == 0) ? k-1 : k ;
	tuu  += (m-1)/2*nr ;
	for (j=(m-1)/2 ; j<nt-1 ; j++) {
	    int ll = 2*j+1 ;
	    double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	    tuu  += nr ;
	}     // Fin de boucle sur theta
	tuu  += nr ; // On saute j=nt-1
    }	// Fin de boucle sur phi	

//## Verif
    assert (tuu == tb->t + (np+1)*nt*nr) ;
	    
    // base de developpement inchangee 
}

			//------------------
			// cas T_LEG_II --
			//----------------

void _poisangu_t_leg_ii(Mtbl_cf* mt, int l, double lambda)
{

    Tbl* tb = mt->t[l] ;	    // pt. sur tbl de travail
    
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    int k, j, i ; 
    // Pour le confort
    int nr = mt->get_mg()->get_nr(l) ;   // Nombre
    int nt = mt->get_mg()->get_nt(l) ;   //	de points
    int np = mt->get_mg()->get_np(l) ;   //	    physiques
    
    double* tuu = tb->t ; 

    // k = 0  :	    cos(phi)
    // -----
     
    for (j=0 ; j<nt-1 ; j++) {
		int ll = 2*j ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
		tuu  += nr ;
    }     // Fin de boucle sur theta
    tuu  += nr ; // On saute j=nt-1

    if (np==1) {
	return ; 
    }

    // k = 1 : on saute
    // -----
    tuu += nt*nr ; 
	
    // k = 2 :	sin(phi)
    // ------
    for (j=0 ; j<nt-1 ; j++) {
		int ll = 2*j+1 ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
		tuu  += nr ;
    }     // Fin de boucle sur theta
    tuu  += nr ; // On saute j=nt-1

    // 3 <= k <= np
    // ------------
    for (k=3 ; k<np+1 ; k++) {
	int m = (k%2 == 0) ? k-1 : k ;
	tuu  += (m+1)/2*nr ;
	for (j=(m+1)/2 ; j<nt-1 ; j++) {
	    int ll = 2*j ;
	    double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	    tuu  += nr ;
	}     // Fin de boucle sur theta
	tuu  += nr ; // On saute j=nt-1
    }	// Fin de boucle sur phi	

//## Verif
    assert (tuu == tb->t + (np+1)*nt*nr) ;
	    
    // base de developpement inchangee 
}

			//------------------
			// cas T_LEG_MP --
			//----------------

void _poisangu_t_leg_mp(Mtbl_cf* mt, int l, double lambda)
{

    Tbl* tb = mt->t[l] ;	    // pt. sur tbl de travail
    
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    int k, j, i ; 
    // Pour le confort
    int nr = mt->get_mg()->get_nr(l) ;   // Nombre
    int nt = mt->get_mg()->get_nt(l) ;   //	de points
    int np = mt->get_mg()->get_np(l) ;   //	    physiques
    
    int np1 = ( np == 1 ) ? 1 : np+1 ; 
	
    double* tuu = tb->t ; 

    // k = 0  :
     
    for (j=0 ; j<nt ; j++) {
		int ll = j ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	tuu  += nr ;
    }     // Fin de boucle sur theta

    // On saute k = 1 : 
    tuu += nt*nr ; 
	
    // k=2,...
    for (k=2 ; k<np1 ; k++) {
	int m = 2*(k/2);
	tuu  += m*nr ;
	for (j=m ; j<nt ; j++) {
	    int ll = j ;
	    double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	    tuu  += nr ;
	}     // Fin de boucle sur theta
    }	// Fin de boucle sur phi	
	    
    // base de developpement inchangee 
}

			//----------------
			// cas T_LEG_MI --
			//----------------

void _poisangu_t_leg_mi(Mtbl_cf* mt, int l, double lambda)
{

    Tbl* tb = mt->t[l] ;	    // pt. sur tbl de travail
    
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    int k, j, i ; 
    // Pour le confort
    int nr = mt->get_mg()->get_nr(l) ;   // Nombre
    int nt = mt->get_mg()->get_nt(l) ;   //	de points
    int np = mt->get_mg()->get_np(l) ;   //	    physiques
    
    int np1 = ( np == 1 ) ? 1 : np+1 ; 
	
    double* tuu = tb->t ; 

    // k = 0  :
     
    for (j=0 ; j<nt-1 ; j++) {
		int ll = j ;
		double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	tuu  += nr ;
    }     // Fin de boucle sur theta
    tuu  += nr ; // On saute j=nt-1

    // On saute k = 1 : 
    tuu += nt*nr ; 
	
    // k=2,...
    for (k=2 ; k<np1 ; k++) {
	int m = 2*((k-1)/2) + 1 ;
	tuu  += m*nr ;
	for (j=m ; j<nt-1 ; j++) {
	    int ll = j ;
	    double xl = - ll*(ll+1) + lambda ;

		if (fabs(xl) < 1.e-14) {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] = 0 ;
			}	
		}
		else {
			for (i=0 ; i<nr ; i++) {
	    		tuu[i] /= xl ;
			}	
		}
	    tuu  += nr ;
	}     // Fin de boucle sur theta
    tuu  += nr ; // On saute j=nt-1
    }	// Fin de boucle sur phi	
	    
    // base de developpement inchangee 
}

}
