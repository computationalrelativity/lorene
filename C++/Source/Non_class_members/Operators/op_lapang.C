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


char op_lapang_C[] = "$Header$" ;

/* 
 * Ensemble des routines de base pour le calcul du laplacien angulaire,
 * c'est-a-dire de l'operateur
 * 
 *  d^2/dtheta^2 + cos(theta)/sin(theta) d/dtheta + 1/sin(theta) d^2/dphi^2
 * 
 * (Utilisation interne)
 * 
 *	void _lapang_XXXX(Mtbl_cf * mt, int l)
 *	mt	pointeur sur le Mtbl_cf d'entree-sortie
 *	l	indice de la zone ou l'on doit effectuer le calcul
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:58  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2000/11/14  15:09:08  eric
 * Traitement du cas np=1 dans T_LEG_PI
 *
 * Revision 2.1  2000/10/04  14:54:59  eric
 * Ajout des bases T_LEG_IP et T_LEG_PI.
 *
 * Revision 2.0  1999/04/26  16:42:04  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>

// Headers Lorene
#include "mtbl_cf.h"
#include "grilles.h"
#include "type_parite.h"


		//--------------------------------------
		// Routine pour les cas non prevus ----
		//------------------------------------

void _lapang_pas_prevu(Mtbl_cf* mt, int l) {
    cout << "Unknwon theta basis in the operator Mtbl_cf::lapang() !" << endl ;
    cout << " basis : " << hex << (mt->base).b[l] << endl ; 
    abort () ;
}

			//---------------
			// cas T_LEG_P --
			//---------------

void _lapang_t_leg_p(Mtbl_cf* mt, int l)
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
	double xl = - ll*(ll+1) ;
	for (i=0 ; i<nr ; i++) {
	    tuu[i] *= xl ;
	}	// Fin de boucle sur r
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
	    double xl = - ll*(ll+1) ;
	    for (i=0 ; i<nr ; i++) {
		tuu[i] *= xl ;
	    }	// Fin de boucle sur r
	    tuu  += nr ;
	}     // Fin de boucle sur theta
    }	// Fin de boucle sur phi	
	    
    // base de developpement inchangee 
}

			//------------------
			// cas T_LEG_PP --
			//----------------

void _lapang_t_leg_pp(Mtbl_cf* mt, int l)
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
	double xl = - ll*(ll+1) ;
	for (i=0 ; i<nr ; i++) {
	    tuu[i] *= xl ;
	}	// Fin de boucle sur r
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
	    double xl = - ll*(ll+1) ;
	    for (i=0 ; i<nr ; i++) {
		tuu[i] *= xl ;
	    }	// Fin de boucle sur r
	    tuu  += nr ;
	}     // Fin de boucle sur theta
    }	// Fin de boucle sur phi	
	    
    // base de developpement inchangee 
}

			//----------------
			// cas T_LEG_I --
			//---------------

void _lapang_t_leg_i(Mtbl_cf* mt, int l)
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
	double xl = - ll*(ll+1) ;
	for (i=0 ; i<nr ; i++) {
	    tuu[i] *= xl ;
	}	// Fin de boucle sur r
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
	    double xl = - ll*(ll+1) ;
	    for (i=0 ; i<nr ; i++) {
		tuu[i] *= xl ;
	    }	// Fin de boucle sur r
	    tuu  += nr ;
	}     // Fin de boucle sur theta
	tuu  += nr ; // On saute j=nt-1
    }	// Fin de boucle sur phi	
	    
    // base de developpement inchangee 
}

			//------------------
			// cas T_LEG_IP --
			//----------------

void _lapang_t_leg_ip(Mtbl_cf* mt, int l)
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
	double xl = - ll*(ll+1) ;
	for (i=0 ; i<nr ; i++) {
	    tuu[i] *= xl ;
	}	// Fin de boucle sur r
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
	    double xl = - ll*(ll+1) ;
	    for (i=0 ; i<nr ; i++) {
		tuu[i] *= xl ;
	    }	// Fin de boucle sur r
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

void _lapang_t_leg_pi(Mtbl_cf* mt, int l)
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
	double xl = - ll*(ll+1) ;
	for (i=0 ; i<nr ; i++) {
	    tuu[i] *= xl ;
	}	// Fin de boucle sur r
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
	double xl = - ll*(ll+1) ;
	for (i=0 ; i<nr ; i++) {
	    tuu[i] *= xl ;
	}	// Fin de boucle sur r
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
	    double xl = - ll*(ll+1) ;
	    for (i=0 ; i<nr ; i++) {
		tuu[i] *= xl ;
	    }	// Fin de boucle sur r
	    tuu  += nr ;
	}     // Fin de boucle sur theta
	tuu  += nr ; // On saute j=nt-1
    }	// Fin de boucle sur phi	

//## Verif
    assert (tuu == tb->t + (np+1)*nt*nr) ;
	    
    // base de developpement inchangee 
}


