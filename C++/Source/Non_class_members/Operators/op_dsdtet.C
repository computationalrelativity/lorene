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


char op_dsdtet_C[] = "$Header$" ;

/*
 * Ensemble des routines de base pour la derivation par rapport a theta
 * (Utilisation interne)
 * 
 *	void _dsdtet_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.2  2003/12/19 16:21:48  j_novak
 * Shadow hunt
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2000/10/04  11:50:44  eric
 * Ajout d' abort() dans le cas non prevu.
 *
 * Revision 2.5  2000/02/24  16:41:17  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 2.4  1999/11/15  16:38:13  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.3  1999/03/01  15:07:01  eric
 * *** empty log message ***
 *
 * Revision 2.2  1999/02/23  11:04:52  hyc
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/22  17:12:06  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/22  15:51:02  hyc
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
void _dsdtet_pas_prevu(Tbl* , int & b) {
    cout << "Unknown theta basis in Mtbl_cf::dsdt() !" << endl ;
    cout << " basis: " << hex << b << endl ;
    abort() ; 
}

// cas T_COS_P
//------------
void _dsdtet_t_cos_p(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_SIN_P ;
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

    // Test sur np pour initialisation eventuelle
    //#pragma critical (loch_dsdtet_t_cos_p)
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cx = (double*)(realloc(cx, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cx[i] = - (2*i) ;
	    }
	}
    }	    // Fin de region critique

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
    
    // k >=2
    for (int k=2 ; k<np+1 ; k++) {
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
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_SIN_P ;
}

// cas T_SIN_P
//------------
void _dsdtet_t_sin_p(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_COS_P ;
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
    static int nt_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cx = (double*)(realloc(cx, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cx[i] = (2*i) ;
	    }
	}
    }	    // Fin de region critique

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
    
    // k >=2
    for (int k=2 ; k<np+1 ; k++) {
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
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COS_P ;
}

// cas T_SIN_I
//------------
void _dsdtet_t_sin_i(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_COS_I ;
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

    // Test sur np pour initialisation eventuelle
    //#pragma critical (loch_dsdtet_t_cos_p)
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cx = (double*)(realloc(cx, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cx[i] = (2*i+1) ;
	    }
	}
    }	    // Fin de region critique

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
    
    // k >=2
    for (int k=2 ; k<np+1 ; k++) {
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
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COS_I ;
}

// cas T_COS_I
//------------
void _dsdtet_t_cos_i(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_SIN_I ;
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

    // Test sur np pour initialisation eventuelle
    //#pragma critical (loch_dsdtet_t_cos_i)
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cx = (double*)(realloc(cx, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cx[i] = - (2*i+1) ;
	    }
	}
    }	    // Fin de region critique

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
    
    // k >=2
    for (int k=2 ; k<np+1 ; k++) {
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
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_SIN_I ;
}

// cas T_COSSIN_CP
//----------------
void _dsdtet_t_cossin_cp(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_COSSIN_SP ;
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
    static int nt_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    //#pragma critical (loch_dsdtet_t_cossin_cp)
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = (double*)(realloc(cxp, nt * sizeof(double))) ;
	cxi = (double*)(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = - (2*i) ;
	    cxi[i] = (2*i+1) ;
	    }
	}
    }	    // Fin de region critique

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
    double* cx[2] ;	// Tableau des Pointeur de coefficient
    
    // Initialisation des pointeurs de coefficients
    cx[0] = cxp ;	// cos pairs pour m pair
    cx[1] = cxi ;	// sin impair pour m impair

    // k = 0
    // Choix de la parite
    double* cxl = cx[0] ;	// Pointeur de coefficients local
    for (int j=0 ; j<nt ; j++) {
	for (int i=0 ; i<nr ; i++) {
	    *xco = cxl[j] * (*xci) ;
	    xci++ ;
	    xco++ ;
	}	// Fin de la Boucle sur r
    }   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	// Choix de la parite
	int m = (k/2) % 2 ;
	cxl = cx[m] ;	// Pointeur de coefficients local
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++) {
		*xco = cxl[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la Boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_SP ;
}

// cas T_COSSIN_SP
//----------------
void _dsdtet_t_cossin_sp(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_COSSIN_CP ;
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

    // Test sur np pour initialisation eventuelle
    //#pragma critical (loch_dsdtet_t_cossin_sp)
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = (double*)(realloc(cxp, nt * sizeof(double))) ;
	cxi = (double*)(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = (2*i) ;
	    cxi[i] = - (2*i+1) ;
	    }
	}
    }	    // Fin de region critique

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
    double* cx[2] ;	// Tableau des Pointeur de coefficient
    
    // Initialisation des pointeurs de coefficients
    cx[0] = cxp ;	// sin pairs pour m pair
    cx[1] = cxi ;	// cos impair pour m impair

    // k = 0
    // Choix de la parite
    double* cxl = cx[0] ;	// Pointeur de coefficients local
    for (int j=0 ; j<nt ; j++) {
	for (int i=0 ; i<nr ; i++) {
	    *xco = cxl[j] * (*xci) ;
	    xci++ ;
	    xco++ ;
	}	// Fin de la Boucle sur r
    }   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	// Choix de la parite
	int m = (k/2) % 2 ;
	cxl = cx[m] ;	// Pointeur de coefficients local
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++) {
		*xco = cxl[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la Boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_CP ;
}

// cas T_COSSIN_CI
//----------------
void _dsdtet_t_cossin_ci(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_COSSIN_SI ;
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

    // Test sur np pour initialisation eventuelle
    //#pragma critical (loch_dsdtet_t_cossin_ci)
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = (double*)(realloc(cxp, nt * sizeof(double))) ;
	cxi = (double*)(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = (2*i) ;
	    cxi[i] = - (2*i+1) ;
	    }
	}
    }	    // Fin de region critique

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
    double* cx[2] ;	// Tableau des Pointeur de coefficient
    
    // Initialisation des pointeurs de coefficients
    cx[0] = cxi ;	// cos impairs pour m pair
    cx[1] = cxp ;	// sin pair pour m impair

    // k = 0
    // Choix de la parite
    double* cxl = cx[0] ;	// Pointeur de coefficients local
    for (int j=0 ; j<nt ; j++) {
	for (int i=0 ; i<nr ; i++) {
	    *xco = cxl[j] * (*xci) ;
	    xci++ ;
	    xco++ ;
	}	// Fin de la Boucle sur r
    }   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	// Choix de la parite
	int m = (k/2) % 2 ;
	cxl = cx[m] ;	// Pointeur de coefficients local
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++) {
		*xco = cxl[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la Boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_SI ;
}

// cas T_COSSIN_SI
//----------------
void _dsdtet_t_cossin_si(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_COSSIN_CI ;
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

    // Test sur np pour initialisation eventuelle
    //#pragma critical (loch_dsdtet_t_cossin_si)
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = (double*)(realloc(cxp, nt * sizeof(double))) ;
	cxi = (double*)(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = - (2*i) ;
	    cxi[i] = (2*i+1) ;
	    }
	}
    }	    // Fin de region critique

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
    double* cx[2] ;	// Tableau des Pointeur de coefficient
    
    // Initialisation des pointeurs de coefficients
    cx[0] = cxi ;	// sin impair pour m pair
    cx[1] = cxp ;	// cos pairs pour m impair

    // k = 0
    // Choix de la parite
    double* cxl = cx[0] ;	// Pointeur de coefficients local
    for (int j=0 ; j<nt ; j++) {
	for (int i=0 ; i<nr ; i++) {
	    *xco = cxl[j] * (*xci) ;
	    xci++ ;
	    xco++ ;
	}	// Fin de la Boucle sur r
    }   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	// Choix de la parite
	int m = (k/2) % 2 ;
	cxl = cx[m] ;	// Pointeur de coefficients local
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++) {
		*xco = cxl[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la Boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_CI ;
}

// cas T_COSSIN_C
//----------------
void _dsdtet_t_cossin_c(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_COSSIN_S ;
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
    static int nt_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    //#pragma critical (loch_dsdtet_t_cossin_cp)
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = (double*)(realloc(cxp, nt * sizeof(double))) ;
	cxi = (double*)(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = - i ;
	    cxi[i] = i ;
	    }
	}
    }	    // Fin de region critique

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
    double* cx[2] ;	// Tableau des Pointeur de coefficient
    
    // Initialisation des pointeurs de coefficients
    cx[0] = cxp ;	// cos pairs pour m pair
    cx[1] = cxi ;	// sin impair pour m impair

    // k = 0
    // Choix de la parite
    double* cxl = cx[0] ;	// Pointeur de coefficients local
    for (int j=0 ; j<nt ; j++) {
	for (int i=0 ; i<nr ; i++) {
	    *xco = cxl[j] * (*xci) ;
	    xci++ ;
	    xco++ ;
	}	// Fin de la Boucle sur r
    }   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	// Choix de la parite
	int m = (k/2) % 2 ;
	cxl = cx[m] ;	// Pointeur de coefficients local
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++) {
		*xco = cxl[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la Boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_S ;
}

// cas T_COSSIN_S
//----------------
void _dsdtet_t_cossin_s(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	b = base_r | base_p | T_COSSIN_C ;
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

    // Test sur np pour initialisation eventuelle
    //#pragma critical (loch_dsdtet_t_cossin_sp)
    {
    if (nt > nt_pre) {
	nt_pre = nt ;
	cxp = (double*)(realloc(cxp, nt * sizeof(double))) ;
	cxi = (double*)(realloc(cxi, nt * sizeof(double))) ;
	for (int i=0 ; i<nt ; i++) {
	    cxp[i] = i ;
	    cxi[i] = - i ;
	    }
	}
    }	    // Fin de region critique

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
    double* cx[2] ;	// Tableau des Pointeur de coefficient
    
    // Initialisation des pointeurs de coefficients
    cx[0] = cxp ;	// sin pairs pour m pair
    cx[1] = cxi ;	// cos impair pour m impair

    // k = 0
    // Choix de la parite
    double* cxl = cx[0] ;	// Pointeur de coefficients local
    for (int j=0 ; j<nt ; j++) {
	for (int i=0 ; i<nr ; i++) {
	    *xco = cxl[j] * (*xci) ;
	    xci++ ;
	    xco++ ;
	}	// Fin de la Boucle sur r
    }   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	// Choix de la parite
	int m = (k/2) % 2 ;
	cxl = cx[m] ;	// Pointeur de coefficients local
	for (int j=0 ; j<nt ; j++) {
	    for (int i=0 ; i<nr ; i++) {
		*xco = cxl[j] * (*xci) ;
		xci++ ;
		xco++ ;
	    }	// Fin de la Boucle sur r
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_C ;
}
