/*
 * Member functions of the Mtbl_cf class for computing the value of a field
 *  at an arbitrary point, when the field is antisymmetric with respect to the
 *  y=0 plane.
 *
 * (see file mtbl_cf.h for the documentation).
 */

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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:00  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:09  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2012/01/17 15:09:14  j_penner
 * using MAX_BASE_2 for the phi coordinate
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2000/09/08  16:07:26  eric
 * Ajout de la base P_COSSIN_I
 *
 * Revision 2.2  2000/03/06  15:59:03  eric
 * *** empty log message ***
 *
 * Revision 2.1  2000/03/06  15:57:29  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/03/06  10:26:54  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Headers Lorene
#include "mtbl_cf.h"
#include "proto.h"

    //-------------------------------------------------------------//
namespace Lorene {
    //	    version for an arbitrary point in (xi,theta',phi')     //
    //-------------------------------------------------------------//

double Mtbl_cf::val_point_asymy(int l, double x, double theta, double phi) 
			 const {

// Routines de sommation
static void (*som_r[MAX_BASE])
	    (double*, const int, const int, const int, const double, double*) ;
static void (*som_tet[MAX_BASE])
	    (double*, const int, const int, const double, double*) ;
static void (*som_phi[MAX_BASE_2])
	    (double*, const int, const double, double*) ;
static int premier_appel = 1 ;

// Initialisations au premier appel
// --------------------------------
    if (premier_appel == 1) {

	premier_appel = 0 ;

	for (int i=0 ; i<MAX_BASE ; i++) {
		if(i%2==0){
	    som_phi[i/2] = som_phi_pas_prevu ;
		}
	    som_tet[i] = som_tet_pas_prevu ;
	    som_r[i]   = som_r_pas_prevu ;
	}

	som_r[R_CHEB >> TRA_R] = som_r_cheb_asymy ;
	som_r[R_CHEBP >> TRA_R] = som_r_chebp ;
	som_r[R_CHEBI >> TRA_R] = som_r_chebi ;
	som_r[R_CHEBU >> TRA_R] = som_r_chebu_asymy ;
	som_r[R_CHEBPIM_P >> TRA_R] = som_r_chebpim_p_asymy ;
	som_r[R_CHEBPIM_I >> TRA_R] = som_r_chebpim_i_asymy ;

	som_tet[T_COS >> TRA_T] = som_tet_cos ;
	som_tet[T_SIN >> TRA_T] = som_tet_sin ;
	som_tet[T_COS_P >> TRA_T] = som_tet_cos_p ;
	som_tet[T_SIN_P >> TRA_T] = som_tet_sin_p ;
	som_tet[T_COSSIN_CP >> TRA_T] = som_tet_cossin_cp_asymy ;
	som_tet[T_COSSIN_CI >> TRA_T] = som_tet_cossin_ci_asymy ;

	som_phi[P_COSSIN >> TRA_P] = som_phi_cossin_asymy ;
	som_phi[P_COSSIN_P >> TRA_P] = som_phi_cossin_p ;
	som_phi[P_COSSIN_I >> TRA_P] = som_phi_cossin_i ;

    }	// fin des operations de premier appel


    assert (etat != ETATNONDEF) ; 

    double resu ;		    // valeur de retour    

// Cas ou tous les coefficients sont nuls :
    if (etat == ETATZERO ) {
	resu = 0 ;
	return resu ; 
    } 
 
// Nombre de points en phi, theta et r :
    int np = mg->get_np(l) ;	 
    int nt = mg->get_nt(l) ;
    int nr = mg->get_nr(l) ;

// Bases de developpement : 
    int base_r = (base.b[l] & MSQ_R) >> TRA_R ;
    int base_t = (base.b[l] & MSQ_T) >> TRA_T ;
    int base_p = (base.b[l] & MSQ_P) >> TRA_P ;
    
//--------------------------------------
// Calcul de la valeur au point demande
//--------------------------------------

// Pointeur sur le tableau contenant les coefficients: 

    assert(etat == ETATQCQ) ; 
    Tbl* tbcf = t[l] ; 
    
    if (tbcf->get_etat() == ETATZERO ) {
	resu = 0 ;
	return resu ; 
    } 


    assert(tbcf->get_etat() == ETATQCQ) ; 

    double* cf = tbcf->t ;

    // Tableaux de travail 
    double* trp = new double [np+2] ;
    double* trtp = new double [(np+2)*nt] ;
 
    if (nr == 1) {
    
// Cas particulier nr = 1 (Fonction purement angulaire) : 
// ----------------------
	
	som_tet[base_t](cf, nt, np, theta, trp) ;   // sommation sur theta
    }
    else {

// Cas general
// -----------

	som_r[base_r](cf, nr, nt, np, x, trtp) ;    // sommation sur r
	som_tet[base_t](trtp, nt, np, theta, trp) ; // sommation sur theta
    }

// Sommation sur phi
// -----------------

    if (np == 1) {
	resu = trp[0] ;		// cas axisymetrique
    }
    else {
	som_phi[base_p](trp, np, phi, &resu) ;	    // sommation sur phi
    }

    // Menage
    delete [] trp ;
    delete [] trtp ;
    
    // Termine
    return resu ;  
 
}
 


    //-------------------------------------------------------------//
    //	    version for an arbitrary point in xi		   //
    //		but collocation point in (theta',phi')		   //
    //-------------------------------------------------------------//

double Mtbl_cf::val_point_jk_asymy(int l, double x, int j0, int k0) const {

// Routines de sommation
static void (*som_r[MAX_BASE])
	    (double*, const int, const int, const int, const double, double*) ;
static int premier_appel = 1 ;

// Initialisations au premier appel
// --------------------------------
    if (premier_appel == 1) {

	premier_appel = 0 ;

	for (int i=0 ; i<MAX_BASE ; i++) {
	    som_r[i]   = som_r_pas_prevu ;
	}

	som_r[R_CHEB >> TRA_R] = som_r_cheb_asymy ;
	som_r[R_CHEBP >> TRA_R] = som_r_chebp ;
	som_r[R_CHEBI >> TRA_R] = som_r_chebi ;
	som_r[R_CHEBU >> TRA_R] = som_r_chebu_asymy ;
	som_r[R_CHEBPIM_P >> TRA_R] = som_r_chebpim_p_asymy ;
	som_r[R_CHEBPIM_I >> TRA_R] = som_r_chebpim_i_asymy ;

    }	// fin des operations de premier appel

    assert (etat != ETATNONDEF) ; 

    double resu ;		    // valeur de retour    

// Cas ou tous les coefficients sont nuls :
    if (etat == ETATZERO ) {
	resu = 0 ;
	return resu ; 
    } 
 
// Nombre de points en phi, theta et r :
    int np = mg->get_np(l) ;	 
    int nt = mg->get_nt(l) ;
    int nr = mg->get_nr(l) ;

// Bases de developpement : 
    int base_r = (base.b[l] & MSQ_R) >> TRA_R ;

//------------------------------------------------------------------------
//  Valeurs des fonctions de base en phi aux points de collocation en phi
//   et des fonctions de base en theta aux points de collocation en theta
//------------------------------------------------------------------------

    Tbl tab_phi = base.phi_functions(l, np) ; 
    Tbl tab_theta = base.theta_functions(l, nt) ; 

    
//--------------------------------------
// Calcul de la valeur au point demande
//--------------------------------------

// Pointeur sur le tableau contenant les coefficients: 

    assert(etat == ETATQCQ) ; 
    Tbl* tbcf = t[l] ; 
    
    if (tbcf->get_etat() == ETATZERO ) {
	resu = 0 ;
	return resu ; 
    } 


    assert(tbcf->get_etat() == ETATQCQ) ; 

    double* cf = tbcf->t ;

    // Tableau de travail 
    double* coef_tp = new double [(np+2)*nt] ;
 

// 1/ Sommation sur r
// ------------------

    som_r[base_r](cf, nr, nt, np, x, coef_tp) ;


// 2/ Sommation sur theta et phi
// -----------------------------
    double* pi = coef_tp ;  // pointeur courant sur les coef en theta et phi

    resu = 0 ; 
 
    if (np > 1) {	// sommation sur phi

    // On saute les coef k=0 (cos(0 phi), k=1 (sin(0 phi) et k=2
    pi += 3*nt ; 

    // Sommation sur le reste des phi (pour k=3,...,np)

	int base_t = base.b[l] & MSQ_T ;

	switch (base_t) {
	    	    	    
	    case T_COSSIN_CP : {

		for (int k=3 ; k<np+1 ; k+=2) { // k+=2 : on saute les cos(m phi)
		    int m_par = (k/2)%2 ;   // 0 pour m pair, 1 pour m impair
		    double somt = 0 ;
		    for (int j=0 ; j<nt ; j++) {
			somt += (*pi) * tab_theta(m_par, j, j0) ;
			pi++ ;  // theta suivant
		    }
		    resu += somt * tab_phi(k, k0) ;

		    // On saute le cos(k*phi) : 
		    pi += nt ;
		}
		break ;
	    }
	 
	    default: {
		cout << "Mtbl_cf::val_point_jk_asymy: unknown theta basis ! " 
		     << endl ;
		abort () ;
	    }
	       
	}   // fin des cas sur base_t 
    
    }	// fin du cas np > 1 


    // Menage
    delete [] coef_tp ;
    
    // Termine
    return resu ;  
 
}
 

}
