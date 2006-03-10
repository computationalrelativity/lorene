/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
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


char op_d2dtet2_C[] = "$Header$" ;

/*
 * Ensemble des routines de base pour la derivation seconde par rapport a theta
 * (Utilisation interne)
 * 
 *	void _d2sdtet2_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.3  2006/03/10 12:45:38  j_novak
 * Use of C++-style cast.
 *
 * Revision 1.2  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/10/04  11:50:20  eric
 * Ajout d' abort() dans le cas non prevu.
 *
 * Revision 2.3  2000/02/24  16:40:42  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 2.2  1999/11/15  16:37:44  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.1  1999/03/01  15:06:00  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/23  11:23:09  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tbl.h"

// Routine pour les cas non prevus
//--------------------------------
void _d2sdtet2_pas_prevu(Tbl* , int & b) {
    cout << "Unknown theta basis in Mtbl_cf::d2sdt2() !" << endl ;
    cout << " basis: " << hex << b << endl ;
    abort() ; 
}

// cas T_COS_P
//------------
void _d2sdtet2_t_cos_p(Tbl* tb, int &)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cx = 0 ;
    static int nt_pre =0 ;

    // Test sur nt pour initialisation eventuelle
    if (nt > nt_pre) {
	nt_pre = nt ;
	cx = reinterpret_cast<double*>(realloc(cx, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cx[i] = - (2*i) * (2*i) ;
	    }
	}

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;

    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }    
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    // k = 0
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cx[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    for (int k=2 ; k<borne_phi ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cx[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas T_SIN_P
//------------
void _d2sdtet2_t_sin_p(Tbl* tb, int &)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cx = 0 ;
    static int nt_pre =0 ;

    // Test sur nt pour initialisation eventuelle
    if (nt > nt_pre) {
	nt_pre = nt ;
	cx = reinterpret_cast<double*>(realloc(cx, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cx[i] = - (2*i) * (2*i) ;
	}
    }

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    for (int k=0 ; k< borne_phi ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cx[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas T_SIN_I
//------------
void _d2sdtet2_t_sin_i(Tbl* tb, int &)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cx = 0 ;
    static int nt_pre =0 ;

    // Test sur nt pour initialisation eventuelle
    if (nt > nt_pre) {
	nt_pre = nt ;
	cx = reinterpret_cast<double*>(realloc(cx, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cx[i] = - (2*i+1) * (2*i+1) ;
	}
    }

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    for (int k=0 ; k< borne_phi ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cx[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas T_COS_I
//------------
void _d2sdtet2_t_cos_i(Tbl* tb, int &)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cx = 0 ;
    static int nt_pre =0 ;

    // Test sur nt pour initialisation eventuelle
    if (nt > nt_pre) {
	nt_pre = nt ;
	cx = reinterpret_cast<double*>(realloc(cx, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cx[i] = - (2*i+1) * (2*i+1) ;
	}
    }

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    for (int k=0 ; k< borne_phi ; k++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cx[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas T_COSSIN_CP
//----------------
void _d2sdtet2_t_cossin_cp(Tbl* tb, int &)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cxp = 0 ;
    static double* cxi = 0 ;
    static int nt_pre =0 ;

    // Test sur nt pour initialisation eventuelle
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = reinterpret_cast<double*>(realloc(cxp, nt * sizeof(double))) ;
	cxi = reinterpret_cast<double*>(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = - (2*i) * (2*i) ;
	    cxi[i] = - (2*i+1) * (2*i+1) ;
	}
    }

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    // Partie cos(pair)
    int k ;
    for (k=0 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxp[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // Partie sin(impair)
    xci = xi + 2*nr*nt ;
    xco = xo + 2*nr*nt ;
    for (k=2 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxi[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas T_COSSIN_SP
//----------------
void _d2sdtet2_t_cossin_sp(Tbl* tb, int &)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cxp = 0 ;
    static double* cxi = 0 ;
    static int nt_pre =0 ;

    // Test sur nt pour initialisation eventuelle
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = reinterpret_cast<double*>(realloc(cxp, nt * sizeof(double))) ;
	cxi = reinterpret_cast<double*>(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = - (2*i) * (2*i) ;
	    cxi[i] = - (2*i+1) * (2*i+1) ;
	}
    }

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    // Partie sin(pair)
    int k ;
    for (k=0 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxp[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // Partie cos(impair)
    xci = xi + 2*nr*nt ;
    xco = xo + 2*nr*nt ;
    for (k=2 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxi[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas T_COSSIN_SI
//----------------
void _d2sdtet2_t_cossin_si(Tbl* tb, int &)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cxp = 0 ;
    static double* cxi = 0 ;
    static int nt_pre =0 ;

    // Test sur nt pour initialisation eventuelle
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = reinterpret_cast<double*>(realloc(cxp, nt * sizeof(double))) ;
	cxi = reinterpret_cast<double*>(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = - (2*i) * (2*i) ;
	    cxi[i] = - (2*i+1) * (2*i+1) ;
	}
    }

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    // Partie sin(impair)
    int k ;
    for (k=0 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxi[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // Partie cos(pair)
    xci = xi + 2*nr*nt ;
    xco = xo + 2*nr*nt ;
    for (k=2 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxp[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas T_COSSIN_C
//----------------
void _d2sdtet2_t_cossin_c(Tbl* tb, int &)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cxp = 0 ;
    static double* cxi = 0 ;
    static int nt_pre =0 ;

    // Test sur nt pour initialisation eventuelle
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = reinterpret_cast<double*>(realloc(cxp, nt * sizeof(double))) ;
	cxi = reinterpret_cast<double*>(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = - i*i ;
	    cxi[i] = - i*i ;
	}
    }

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    // Partie cos(pair)
    int k ;
    for (k=0 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxp[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // Partie sin(impair)
    xci = xi + 2*nr*nt ;
    xco = xo + 2*nr*nt ;
    for (k=2 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxi[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas T_COSSIN_S
//----------------
void _d2sdtet2_t_cossin_s(Tbl* tb, int &)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // Variables statiques
    static double* cxp = 0 ;
    static double* cxi = 0 ;
    static int nt_pre =0 ;

    // Test sur nt pour initialisation eventuelle
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = reinterpret_cast<double*>(realloc(cxp, nt * sizeof(double))) ;
	cxi = reinterpret_cast<double*>(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = - i*i ;
	    cxi[i] = - i*i ;
	}
    }

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    // Partie sin(pair)
    int k ;
    for (k=0 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxp[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // Partie cos(impair)
    xci = xi + 2*nr*nt ;
    xco = xo + 2*nr*nt ;
    for (k=2 ; k<np+1 ; k += 4) {
     for (int m=0 ; m<2 ; m++) {
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++ ) {
		*xco = cxi[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
     }	  // Fin de la boucle intermediaire
    xci += 2*nr*nt ;
    xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}
