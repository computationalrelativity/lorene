/*
 * Member functions of the Mtbl_cf class for computing the value of a field
 *  at an arbitrary point
 *
 * (see file mtbl_cf.h for the documentation).
 */

/*
 *   Copyright (c) 1999-2002 Eric Gourgoulhon
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


char mtbl_cf_val_point_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.12  2009/10/08 16:21:16  j_novak
 * Addition of new bases T_COS and T_SIN.
 *
 * Revision 1.11  2007/12/20 09:11:08  jl_cornou
 * Correction of an error in op_sxpun about Jacobi(0,2) polynomials
 *
 * Revision 1.10  2007/12/11 15:28:16  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.9  2007/10/23 17:15:13  j_novak
 * Added the bases T_COSSIN_C and T_COSSIN_S
 *
 * Revision 1.8  2006/06/06 14:57:01  j_novak
 * Summation functions for angular coefficients at xi=+/-1.
 *
 * Revision 1.7  2006/05/30 08:30:15  n_vasset
 * Implementation of sine-like bases (T_SIN_P, T_SIN_I, T_COSSIN_SI, etc...).
 *
 * Revision 1.6  2005/05/27 14:55:00  j_novak
 * Added new bases T_COSSIN_CI and T_COS_I
 *
 * Revision 1.5  2005/02/16 15:10:39  m_forot
 * Correct the case T_COSSIN_C
 *
 * Revision 1.4  2004/12/17 13:35:03  m_forot
 * Add the case T_LEG
 *
 * Revision 1.3  2002/05/11 12:37:31  e_gourgoulhon
 * Added basis T_COSSIN_SI.
 *
 * Revision 1.2  2002/05/05 16:22:33  e_gourgoulhon
 * Added the case of the theta basis T_COSSIN_SP.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.7  2000/09/08  16:26:43  eric
 * Ajout de la base T_SIN_I.
 *
 * Revision 1.6  2000/09/08  16:07:16  eric
 * Ajout de la base P_COSSIN_I
 *
 * Revision 1.5  2000/09/06  14:00:19  eric
 * Ajout de la base T_COS_I.
 *
 * Revision 1.4  1999/12/29  13:11:42  eric
 * Ajout de la fonction val_point_jk.
 *
 * Revision 1.3  1999/12/07  15:10:45  eric
 * *** empty log message ***
 *
 * Revision 1.2  1999/12/07  14:52:34  eric
 * Changement ordre des arguments (phi,theta,xi) --> (xi,theta,phi)
 *
 * Revision 1.1  1999/12/06  16:47:39  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "mtbl_cf.h"
#include "proto.h"

    //-------------------------------------------------------------//
    //	    version for an arbitrary point in (xi,theta',phi')     //
    //-------------------------------------------------------------//

double Mtbl_cf::val_point(int l, double x, double theta, double phi) const {

// Routines de sommation
static void (*som_r[MAX_BASE])
	    (double*, const int, const int, const int, const double, double*) ;
static void (*som_tet[MAX_BASE])
	    (double*, const int, const int, const double, double*) ;
static void (*som_phi[MAX_BASE])
	    (double*, const int, const double, double*) ;
static int premier_appel = 1 ;

// Initialisations au premier appel
// --------------------------------
    if (premier_appel == 1) {

	premier_appel = 0 ;

	for (int i=0 ; i<MAX_BASE ; i++) {
	    som_phi[i] = som_phi_pas_prevu ;
	    som_tet[i] = som_tet_pas_prevu ;
	    som_r[i]   = som_r_pas_prevu ;
	}

	som_r[R_CHEB >> TRA_R] = som_r_cheb ;
	som_r[R_CHEBP >> TRA_R] = som_r_chebp ;
	som_r[R_CHEBI >> TRA_R] = som_r_chebi ;
	som_r[R_CHEBU >> TRA_R] = som_r_chebu ;
	som_r[R_CHEBPIM_P >> TRA_R] = som_r_chebpim_p ;
	som_r[R_CHEBPIM_I >> TRA_R] = som_r_chebpim_i ;
	som_r[R_CHEBPI_P >> TRA_R] = som_r_chebpi_p ;
	som_r[R_CHEBPI_I >> TRA_R] = som_r_chebpi_i ;
	som_r[R_JACO02 >> TRA_R] = som_r_jaco02 ;

	som_tet[T_COS >> TRA_T] = som_tet_cos ;
	som_tet[T_SIN >> TRA_T] = som_tet_sin ;
	som_tet[T_COS_P >> TRA_T] = som_tet_cos_p ;
	som_tet[T_COS_I >> TRA_T] = som_tet_cos_i ;
	som_tet[T_SIN_P >> TRA_T] = som_tet_sin_p ;
	som_tet[T_SIN_I >> TRA_T] = som_tet_sin_i ;
	som_tet[T_COSSIN_CP >> TRA_T] = som_tet_cossin_cp ;
	som_tet[T_COSSIN_CI >> TRA_T] = som_tet_cossin_ci ;
	som_tet[T_COSSIN_SP >> TRA_T] = som_tet_cossin_sp ;
	som_tet[T_COSSIN_SI >> TRA_T] = som_tet_cossin_si ;
	som_tet[T_COSSIN_C >> TRA_T] = som_tet_cossin_c ;
	som_tet[T_COSSIN_S >> TRA_T] = som_tet_cossin_s ;

	som_phi[P_COSSIN >> TRA_P] = som_phi_cossin ;
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

double Mtbl_cf::val_point_jk(int l, double x, int j0, int k0) const {

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

	som_r[R_CHEB >> TRA_R] = som_r_cheb ;
	som_r[R_CHEBP >> TRA_R] = som_r_chebp ;
	som_r[R_CHEBI >> TRA_R] = som_r_chebi ;
	som_r[R_CHEBU >> TRA_R] = som_r_chebu ;
	som_r[R_CHEBPIM_P >> TRA_R] = som_r_chebpim_p ;
	som_r[R_CHEBPIM_I >> TRA_R] = som_r_chebpim_i ;
	som_r[R_CHEBPI_P >> TRA_R] = som_r_chebpi_p ;
	som_r[R_CHEBPI_I >> TRA_R] = som_r_chebpi_i ;
	som_r[R_JACO02 >> TRA_R] = som_r_jaco02 ;

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

// Sommation sur le premier phi, k=0

    double somt = 0 ;
    for (int j=0 ; j<nt ; j++) {
	somt += (*pi) * tab_theta(0, j, j0) ;
	pi++ ;	// theta suivant
    }
    resu = somt * tab_phi(0, k0)	;
 
    if (np > 1) {	// sommation sur phi

    // On saute le phi suivant (sin(0)), k=1
	pi += nt ;

    // Sommation sur le reste des phi (pour k=2,...,np)

	int base_t = base.b[l] & MSQ_T ;

	switch (base_t) {
	    	    
	    case T_COS :
	    case T_SIN :
	    case T_SIN_P :
	    case T_SIN_I :
	    case T_COS_I : 
	    case T_COS_P : {

		for (int k=2 ; k<np+1 ; k++) {
		    somt = 0 ;
		    for (int j=0 ; j<nt ; j++) {
			somt += (*pi) * tab_theta(0, j, j0) ;
			pi++ ;  // theta suivant
		    }
		    resu += somt * tab_phi(k, k0) ;
		}
		break ;
	    } 
	    
	    case T_COSSIN_C : 
	    case T_COSSIN_S : 
	    case T_COSSIN_SP : 
	    case T_COSSIN_SI : 
	    case T_COSSIN_CI : 
	    case T_COSSIN_CP : {

		for (int k=2 ; k<np+1 ; k++) {
		    int m_par = (k/2)%2 ;   // 0 pour m pair, 1 pour m impair
		    somt = 0 ;
		    for (int j=0 ; j<nt ; j++) {
			somt += (*pi) * tab_theta(m_par, j, j0) ;
			pi++ ;  // theta suivant
		    }
		    resu += somt * tab_phi(k, k0) ;
		}
		break ;
	    }
	 
	    default: {
		cout << "Mtbl_cf::val_point_jk: unknown theta basis ! " << endl ;
		abort () ;
	    }
	       
	}   // fin des cas sur base_t 
    
    }	// fin du cas np > 1 


    // Menage
    delete [] coef_tp ;
    
    // Termine
    return resu ;  
 
}
 

    //-------------------------------------------------------------//
    //	    version for xi = 1	                          	   //
    //		and collocation point in (theta',phi')		   //
    //-------------------------------------------------------------//

double Mtbl_cf::val_out_bound_jk(int l, int j0, int k0) const {

#ifndef NDEBUG
// Bases de developpement : 
    int base_r = (base.b[l] & MSQ_R) >> TRA_R ;
    assert((base_r == R_CHEB) || (base_r == R_CHEBU) || (base_r == R_CHEBP)
	   || (base_r == R_CHEBI) || (base_r == R_CHEBPIM_P) || (base_r == R_CHEBPIM_I)
	   || (base_r == R_CHEBPI_P) || (base_r == R_CHEBPI_I) || (base_r == R_JACO02)) ;
#endif

    int nr = mg->get_nr(l) ;
    double resu = 0 ;
    for (int i=0; i<nr; i++)
	resu += operator()(l, k0, j0, i) ;

    return resu ;
}


    //-------------------------------------------------------------//
    //	    version for xi = -1	                          	   //
    //		and collocation point in (theta',phi')		   //
    //-------------------------------------------------------------//

double Mtbl_cf::val_in_bound_jk(int l, int j0, int k0) const {

#ifndef NDEBUG
// Bases de developpement : 
    int base_r = (base.b[l] & MSQ_R) >> TRA_R ;
    assert((base_r == R_CHEB) || (base_r == R_CHEBU)) ;
#endif

    int nr = mg->get_nr(l) ;
    double resu = 0 ;
    int pari = 1 ;
    for (int i=0; i<nr; i++) {
	resu += pari*operator()(l, k0, j0, i) ;
	pari = - pari ;
    }

    return resu ;
}

