/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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
 * Ensemble des routines de base de derivation seconde par rapport a r
 * (Utilisation interne)
 * 
 *	void _d2sdx2_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.8  2016/12/05 16:18:07  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2015/03/05 08:49:32  j_novak
 * Implemented operators with Legendre bases.
 *
 * Revision 1.6  2014/10/13 08:53:24  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2013/06/14 15:54:06  j_novak
 * Inclusion of Legendre bases.
 *
 * Revision 1.4  2008/08/27 08:50:10  jl_cornou
 * Added Jacobi(0,2) polynomials
 *
 * Revision 1.3  2007/12/11 15:28:18  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.2  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.7  2000/02/24  16:40:50  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 2.6  1999/11/15  16:37:53  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.5  1999/10/11  15:33:31  phil
 * *** empty log message ***
 *
 * Revision 2.4  1999/10/11  14:25:47  phil
 * realloc -> delete + new
 *
 * Revision 2.3  1999/09/13  11:30:23  phil
 * gestion du cas k==1
 *
 * Revision 2.2  1999/03/01  15:06:31  eric
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/23  10:39:13  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Fichier includes
#include "tbl.h"

namespace Lorene {
void _d2sdx2_1d_r_jaco02(int, double*, double*) ;

// Prototypage
void c_est_pas_fait(char * ) ;

// Routine pour les cas non prevus
//--------------------------------
void _d2sdx2_pas_prevu(Tbl* , int & b) {
    cout << "d2sdx2 pas prevu..." << endl ;
    cout << " base: " << b << endl ;
    abort () ;
}

// cas R_CHEBU - dzpuis = 0
//-------------------------
void _d2sdx2_r_chebu_0(Tbl *tb, int & )
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
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1[i] =  (i+2)*(i+2)*(i+2) ;
	    cx2[i] =  (i+2) ;
	    cx3[i] =  i*i ;
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som1, som2 ;
	   
	    xco[nr-1] = 0 ;
	    som1 = (nr-1)*(nr-1)*(nr-1) * xci[nr-1] ;
	    som2 = (nr-1) * xci[nr-1] ;
	    xco[nr-3] = som1 - (nr-3)*(nr-3)*som2 ;
	    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
		som1 += cx1[i] * xci[i+2] ;
		som2 += cx2[i] * xci[i+2] ;
		xco[i] = som1 - cx3[i] * som2 ;
	    }	// Fin de la premiere boucle sur r
	    xco[nr-2] = 0 ;
	    som1 = (nr-2)*(nr-2)*(nr-2) * xci[nr-2] ;
	    som2 = (nr-2) * xci[nr-2] ;
	    xco[nr-4] = som1 - (nr-4)*(nr-4)*som2 ;
	    for (int i = nr-6 ; i >= 0 ; i -= 2 ) {
		som1 += cx1[i] * xci[i+2] ;
		som2 += cx2[i] * xci[i+2] ;
		xco[i] = som1 - cx3[i] * som2 ;
	    }	// Fin de la deuxieme boucle sur r
	    xco[0] *= .5 ;
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_CHEB
//-----------
void _d2sdx2_r_cheb(Tbl *tb, int & )
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
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;  
    
	for (int i=0 ; i<nr ; i++) {
	    cx1[i] =  (i+2)*(i+2)*(i+2) ;
	    cx2[i] =  (i+2) ;
	    cx3[i] =  i*i ;
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som1, som2 ;
	    
	    xco[nr-1] = 0 ;
	    som1 = (nr-1)*(nr-1)*(nr-1) * xci[nr-1] ;
	    som2 = (nr-1) * xci[nr-1] ;
	    xco[nr-3] = som1 - (nr-3)*(nr-3)*som2 ;
	    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
		som1 += cx1[i] * xci[i+2] ;
		som2 += cx2[i] * xci[i+2] ;
		xco[i] = som1 - cx3[i] * som2 ;
	    }	// Fin de la premiere boucle sur r
	    xco[nr-2] = 0 ;
	    som1 = (nr-2)*(nr-2)*(nr-2) * xci[nr-2] ;
	    som2 = (nr-2) * xci[nr-2] ;
	    xco[nr-4] = som1 - (nr-4)*(nr-4)*som2 ;
	    for (int i = nr-6 ; i >= 0 ; i -= 2 ) {
		som1 += cx1[i] * xci[i+2] ;
		som2 += cx2[i] * xci[i+2] ;
		xco[i] = som1 - cx3[i] * som2 ;
	    }	// Fin de la deuxieme boucle sur r
	    xco[0] *= .5 ;
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_CHEBP
//------------
void _d2sdx2_r_chebp(Tbl *tb, int & )
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
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;  
	for (int i=0 ; i<nr ; i++) {
	    cx1[i] =  8*(i+1)*(i+1)*(i+1) ;
	    cx2[i] =  2*(i+1) ;
	    cx3[i] =  4*i*i ;
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som1, som2 ;
	    
	    xco[nr-1] = 0 ;
	    som1 = 8*(nr-1)*(nr-1)*(nr-1) * xci[nr-1] ;
	    som2 = 2*(nr-1) * xci[nr-1] ;
	    xco[nr-2] = som1 - 4*(nr-2)*(nr-2)*som2 ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		som1 += cx1[i] * xci[i+1] ;
		som2 += cx2[i] * xci[i+1] ;
		xco[i] = som1 - cx3[i] * som2 ;
	    }	// Fin de la boucle sur r
	    xco[0] *= .5 ;
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_CHEBI
//------------
void _d2sdx2_r_chebi(Tbl *tb, int & )
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
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;   
   
	for (int i=0 ; i<nr ; i++) {
	    cx1[i] =  (2*i+3)*(2*i+3)*(2*i+3) ;
	    cx2[i] =  (2*i+3) ;
	    cx3[i] =  (2*i+1)*(2*i+1) ;
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
    
    for (int k=0 ; k<np+1 ; k++)  
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som1, som2 ;
	    
	    xco[nr-1] = 0 ;
	    som1 = (2*nr-1)*(2*nr-1)*(2*nr-1) * xci[nr-1] ;
	    som2 = (2*nr-1) * xci[nr-1] ;
	    xco[nr-2] = som1 - (2*nr-3)*(2*nr-3)*som2 ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		som1 += cx1[i] * xci[i+1] ;
		som2 += cx2[i] * xci[i+1] ;
		xco[i] = som1 - cx3[i] * som2 ;
	    }	// Fin de la boucle sur r
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_CHEBPIM_P
//----------------
void _d2sdx2_r_chebpim_p(Tbl *tb, int & )
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
    static double* cx1p = 0x0 ;
    static double* cx2p = 0x0 ;
    static double* cx3p = 0x0 ;
    static double* cx1i = 0x0 ;
    static double* cx2i = 0x0 ;
    static double* cx3i = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1p != 0x0) delete [] cx1p ;
	if (cx2p != 0x0) delete [] cx2p ;
	if (cx3p != 0x0) delete [] cx3p ;
	if (cx1i != 0x0) delete [] cx1i ;
	if (cx2i != 0x0) delete [] cx2i ;
	if (cx3i != 0x0) delete [] cx3i ;
	cx1p = new double[nr] ;
	cx2p = new double[nr] ;
	cx3p = new double[nr] ;
	cx1i = new double[nr] ;
	cx2i = new double[nr] ;
	cx3i = new double[nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1p[i] =  8*(i+1)*(i+1)*(i+1) ;
	    cx2p[i] =  2*(i+1) ;
	    cx3p[i] =  4*i*i ;

	    cx1i[i] =  (2*i+3)*(2*i+3)*(2*i+3) ;
	    cx2i[i] =  (2*i+3) ;
	    cx3i[i] =  (2*i+1)*(2*i+1) ;
	}
    }
    
    double* cx1t[2] ;
    double* cx2t[2] ;
    double* cx3t[2] ;
    cx1t[0] = cx1p ; cx1t[1] = cx1i ;
    cx2t[0] = cx2p ; cx2t[1] = cx2i ;
    cx3t[0] = cx3p ; cx3t[1] = cx3i ;

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci ;	// Pointeurs
    double* xco ;	//  courants
    
    // Position depart
    xci = xi ;
    xco = xo ;
    
    double *cx1, *cx2, *cx3 ;

    // k = 0
    cx1 = cx1t[0] ;
    cx2 = cx2t[0] ;
    cx3 = cx3t[0] ;
	    for (int j=0 ; j<nt ; j++) {

		double som1 = 0 ;
		double som2 = 0 ;
	    
		xco[nr-1] = 0 ;
		for (int i = nr-2 ; i >= 0 ; i-- ) {
		    som1 += cx1[i] * xci[i+1] ;
		    som2 += cx2[i] * xci[i+1] ;
		    xco[i] = som1 - cx3[i] * som2 ;
		}	// Fin de la boucle sur r
		xco[0] *= .5 ;	// normalisation	    
		xci += nr ;	// au
		xco += nr ;	//  suivant
	    }   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	int m = (k/2) % 2 ;
	cx1 = cx1t[m] ;
	cx2 = cx2t[m] ;
	cx3 = cx3t[m] ;
	    for (int j=0 ; j<nt ; j++) {

		double som1 = 0 ;
		double som2 = 0 ;
	    
		xco[nr-1] = 0 ;
		for (int i = nr-2 ; i >= 0 ; i-- ) {
		    som1 += cx1[i] * xci[i+1] ;
		    som2 += cx2[i] * xci[i+1] ;
		    xco[i] = som1 - cx3[i] * som2 ;
		}	// Fin de la boucle sur r
		if (m == 0) xco[0] *= .5 ;  // normalisation eventuelle
		xci += nr ; // au
		xco += nr ; //	suivant
	    }   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_CHEBPIM_I
//----------------
void _d2sdx2_r_chebpim_i(Tbl *tb, int & )
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
    static double* cx1p = 0x0 ;
    static double* cx2p = 0x0 ;
    static double* cx3p = 0x0 ;
    static double* cx1i = 0x0 ;
    static double* cx2i = 0x0 ;
    static double* cx3i = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1p != 0x0) delete [] cx1p ;
	if (cx2p != 0x0) delete [] cx2p ;
	if (cx3p != 0x0) delete [] cx3p ;
	if (cx1i != 0x0) delete [] cx1i ;
	if (cx2i != 0x0) delete [] cx2i ;
	if (cx3i != 0x0) delete [] cx3i ;
	cx1p = new double[nr] ;
	cx2p = new double[nr] ;
	cx3p = new double[nr] ;
	cx1i = new double[nr] ;
	cx2i = new double[nr] ;
	cx3i = new double[nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1p[i] =  8*(i+1)*(i+1)*(i+1) ;
	    cx2p[i] =  2*(i+1) ;
	    cx3p[i] =  4*i*i ;

	    cx1i[i] =  (2*i+3)*(2*i+3)*(2*i+3) ;
	    cx2i[i] =  (2*i+3) ;
	    cx3i[i] =  (2*i+1)*(2*i+1) ;
	}
    }

    double* cx1t[2] ;
    double* cx2t[2] ;
    double* cx3t[2] ;
    cx1t[1] = cx1p ; cx1t[0] = cx1i ;
    cx2t[1] = cx2p ; cx2t[0] = cx2i ;
    cx3t[1] = cx3p ; cx3t[0] = cx3i ;

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci ;	// Pointeurs
    double* xco ;	//  courants
    
    // Position depart
    xci = xi ;
    xco = xo ;

    double *cx1, *cx2, *cx3 ;

    // k = 0
    cx1 = cx1t[0] ;
    cx2 = cx2t[0] ;
    cx3 = cx3t[0] ;
	    for (int j=0 ; j<nt ; j++) {

		double som1 = 0 ;
		double som2 = 0 ;
	    
		xco[nr-1] = 0 ;
		for (int i = nr-2 ; i >= 0 ; i-- ) {
		    som1 += cx1[i] * xci[i+1] ;
		    som2 += cx2[i] * xci[i+1] ;
		    xco[i] = som1 - cx3[i] * som2 ;
		}	// Fin de la boucle sur r	    
		xci += nr ;
		xco += nr ;
	    }   // Fin de la boucle sur theta

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	int m = (k/2) % 2 ;
	cx1 = cx1t[m] ;
	cx2 = cx2t[m] ;
	cx3 = cx3t[m] ;
	    for (int j=0 ; j<nt ; j++) {

		double som1 = 0 ;
		double som2 = 0 ;
	    
		xco[nr-1] = 0 ;
		for (int i = nr-2 ; i >= 0 ; i-- ) {
		    som1 += cx1[i] * xci[i+1] ;
		    som2 += cx2[i] * xci[i+1] ;
		    xco[i] = som1 - cx3[i] * som2 ;
		}	// Fin de la boucle sur r
		if (m == 1) xco[0] *= .5 ;
		xci += nr ;
		xco += nr ;
	    }   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_CHEBPI_P
//----------------
void _d2sdx2_r_chebpi_p(Tbl *tb, int & )
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
    static double* cx1p = 0x0 ;
    static double* cx2p = 0x0 ;
    static double* cx3p = 0x0 ;
    static double* cx1i = 0x0 ;
    static double* cx2i = 0x0 ;
    static double* cx3i = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1p != 0x0) delete [] cx1p ;
	if (cx2p != 0x0) delete [] cx2p ;
	if (cx3p != 0x0) delete [] cx3p ;
	if (cx1i != 0x0) delete [] cx1i ;
	if (cx2i != 0x0) delete [] cx2i ;
	if (cx3i != 0x0) delete [] cx3i ;
	cx1p = new double[nr] ;
	cx2p = new double[nr] ;
	cx3p = new double[nr] ;
	cx1i = new double[nr] ;
	cx2i = new double[nr] ;
	cx3i = new double[nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1p[i] =  8*(i+1)*(i+1)*(i+1) ;
	    cx2p[i] =  2*(i+1) ;
	    cx3p[i] =  4*i*i ;

	    cx1i[i] =  (2*i+3)*(2*i+3)*(2*i+3) ;
	    cx2i[i] =  (2*i+3) ;
	    cx3i[i] =  (2*i+1)*(2*i+1) ;
	}
    }
    
    double* cx1t[2] ;
    double* cx2t[2] ;
    double* cx3t[2] ;
    cx1t[0] = cx1p ; cx1t[1] = cx1i ;
    cx2t[0] = cx2p ; cx2t[1] = cx2i ;
    cx3t[0] = cx3p ; cx3t[1] = cx3i ;

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci ;	// Pointeurs
    double* xco ;	//  courants
    
    // Position depart
    xci = xi ;
    xco = xo ;
    
    double *cx1, *cx2, *cx3 ;

    // k = 0
    for (int j=0 ; j<nt ; j++) {
      int l = j % 2 ;
      cx1 = cx1t[l] ;
      cx2 = cx2t[l] ;
      cx3 = cx3t[l] ;
      double som1 = 0 ;
      double som2 = 0 ;
      
      xco[nr-1] = 0 ;
      for (int i = nr-2 ; i >= 0 ; i-- ) {
	som1 += cx1[i] * xci[i+1] ;
	som2 += cx2[i] * xci[i+1] ;
	xco[i] = som1 - cx3[i] * som2 ;
      }	// Fin de la boucle sur r
      if (l == 0) xco[0] *= .5 ;	// normalisation	    
      xci += nr ;	// au
      xco += nr ;	//  suivant
    }   // Fin de la boucle sur theta
    
    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
      for (int j=0 ; j<nt ; j++) {
	int l = j % 2 ;
	cx1 = cx1t[l] ;
	cx2 = cx2t[l] ;
	cx3 = cx3t[l] ;
	double som1 = 0 ;
	double som2 = 0 ;
	
	xco[nr-1] = 0 ;
	for (int i = nr-2 ; i >= 0 ; i-- ) {
	  som1 += cx1[i] * xci[i+1] ;
	  som2 += cx2[i] * xci[i+1] ;
	  xco[i] = som1 - cx3[i] * som2 ;
	}	// Fin de la boucle sur r
	if (l == 0) xco[0] *= .5 ;  // normalisation eventuelle
	xci += nr ; // au
	xco += nr ; //	suivant
      }   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_CHEBPI_I
//----------------
void _d2sdx2_r_chebpi_i(Tbl *tb, int & )
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
    static double* cx1p = 0x0 ;
    static double* cx2p = 0x0 ;
    static double* cx3p = 0x0 ;
    static double* cx1i = 0x0 ;
    static double* cx2i = 0x0 ;
    static double* cx3i = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1p != 0x0) delete [] cx1p ;
	if (cx2p != 0x0) delete [] cx2p ;
	if (cx3p != 0x0) delete [] cx3p ;
	if (cx1i != 0x0) delete [] cx1i ;
	if (cx2i != 0x0) delete [] cx2i ;
	if (cx3i != 0x0) delete [] cx3i ;
	cx1p = new double[nr] ;
	cx2p = new double[nr] ;
	cx3p = new double[nr] ;
	cx1i = new double[nr] ;
	cx2i = new double[nr] ;
	cx3i = new double[nr] ;
	for (int i=0 ; i<nr ; i++) {
	    cx1p[i] =  8*(i+1)*(i+1)*(i+1) ;
	    cx2p[i] =  2*(i+1) ;
	    cx3p[i] =  4*i*i ;

	    cx1i[i] =  (2*i+3)*(2*i+3)*(2*i+3) ;
	    cx2i[i] =  (2*i+3) ;
	    cx3i[i] =  (2*i+1)*(2*i+1) ;
	}
    }

    double* cx1t[2] ;
    double* cx2t[2] ;
    double* cx3t[2] ;
    cx1t[1] = cx1p ; cx1t[0] = cx1i ;
    cx2t[1] = cx2p ; cx2t[0] = cx2i ;
    cx3t[1] = cx3p ; cx3t[0] = cx3i ;

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci ;	// Pointeurs
    double* xco ;	//  courants
    
    // Position depart
    xci = xi ;
    xco = xo ;

    double *cx1, *cx2, *cx3 ;

    // k = 0
    for (int j=0 ; j<nt ; j++) {
      int l = j % 2 ;
      cx1 = cx1t[l] ;
      cx2 = cx2t[l] ;
      cx3 = cx3t[l] ;
      double som1 = 0 ;
      double som2 = 0 ;
      
      xco[nr-1] = 0 ;
      for (int i = nr-2 ; i >= 0 ; i-- ) {
	som1 += cx1[i] * xci[i+1] ;
	som2 += cx2[i] * xci[i+1] ;
	xco[i] = som1 - cx3[i] * som2 ;
      }	// Fin de la boucle sur r
      if (l == 1) xco[0] *= .5 ;
      xci += nr ;
      xco += nr ;
    }   // Fin de la boucle sur theta
    
    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
      for (int j=0 ; j<nt ; j++) {
	int l = j % 2 ;
	cx1 = cx1t[l] ;
	cx2 = cx2t[l] ;
	cx3 = cx3t[l] ;
	double som1 = 0 ;
	double som2 = 0 ;
	
	xco[nr-1] = 0 ;
	for (int i = nr-2 ; i >= 0 ; i-- ) {
	  som1 += cx1[i] * xci[i+1] ;
	  som2 += cx2[i] * xci[i+1] ;
	  xco[i] = som1 - cx3[i] * som2 ;
	}	// Fin de la boucle sur r
	if (l == 1) xco[0] *= .5 ;
	xci += nr ;
	xco += nr ;
      }   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}


// cas R_LEG
//-----------
void _d2sdx2_r_leg(Tbl *tb, int & )
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
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;  
    
	for (int i=0 ; i<nr ; i++) {
	  cx1[i] =  (i+2)*(i+3) ;
	  cx2[i] =  i*(i+1) ;
	  cx3[i] =  double(i) + 0.5 ;
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som1, som2 ;
	    
	    xco[nr-1] = 0 ;
	    som1 = (nr-1)*nr * xci[nr-1] ;
	    som2 = xci[nr-1] ;
	    if (nr > 2) 
	      xco[nr-3] = (double(nr) -2.5) * (som1 - (nr-3)*(nr-2)*som2) ;
	    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
		som1 += cx1[i] * xci[i+2] ;
		som2 += xci[i+2] ;
		xco[i] = cx3[i]*(som1 - cx2[i] * som2) ;
	    }	// Fin de la premiere boucle sur r
	    if (nr > 1) xco[nr-2] = 0 ;
	    if (nr > 3) {
	      som1 = (nr-2)*(nr-1)* xci[nr-2] ;
	      som2 = xci[nr-2] ;
	      xco[nr-4] = (double(nr) - 3.5) * (som1 - (nr-4)*(nr-3)*som2) ;
	    }
	    for (int i = nr-6 ; i >= 0 ; i -= 2 ) {
		som1 += cx1[i] * xci[i+2] ;
		som2 += xci[i+2] ;
		xco[i] = cx3[i]*(som1 - cx2[i] * som2) ;
	    }	// Fin de la deuxieme boucle sur r
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_LEGP
//------------
void _d2sdx2_r_legp(Tbl *tb, int & )
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
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;  
	for (int i=0 ; i<nr ; i++) {
	  cx1[i] =  (2*i+2)*(2*i+3) ;
	  cx2[i] =  2*i*(2*i+1) ;
	  cx3[i] =  double(2*i) + 0.5 ;
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som1, som2 ;
	    
	    xco[nr-1] = 0 ;
	    som1 = (2*nr-2)*(2*nr-1)* xci[nr-1] ;
	    som2 = xci[nr-1] ;
	    if (nr > 1)
	      xco[nr-2] = (double(2*nr) - 1.5)*(som1 - 2*(nr-2)*(2*nr-1)*som2) ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		som1 += cx1[i] * xci[i+1] ;
		som2 += xci[i+1] ;
		xco[i] = cx3[i]*(som1 - cx2[i] * som2) ;
	    }	// Fin de la boucle sur r
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_LEGI
//------------
void _d2sdx2_r_legi(Tbl *tb, int & )
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
    static double* cx1 = 0x0 ;
    static double* cx2 = 0x0 ;
    static double* cx3 = 0x0 ;
    static int nr_pre = 0 ;

    // Test sur np pour initialisation eventuelle
    if (nr > nr_pre) {
	nr_pre = nr ;
	if (cx1 != 0x0) delete [] cx1 ;
	if (cx2 != 0x0) delete [] cx2 ;
	if (cx3 != 0x0) delete [] cx3 ;
	cx1 = new double [nr] ;
	cx2 = new double [nr] ;
	cx3 = new double [nr] ;   
   
	for (int i=0 ; i<nr ; i++) {
	  cx1[i] =  (2*i+3)*(2*i+4) ;
	  cx2[i] =  (2*i+1)*(2*i+2) ;
	  cx3[i] =  double(2*i) + 1.5 ;
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
    
    for (int k=0 ; k<np+1 ; k++)  
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som1, som2 ;
	    
	    xco[nr-1] = 0 ;
	    som1 = (2*nr-1)*(2*nr) * xci[nr-1] ;
	    som2 = xci[nr-1] ;
	    if (nr > 1) 
	      xco[nr-2] = (double(nr) - 1.5)*(som1 - (2*nr-3)*(2*nr-2)*som2) ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		som1 += cx1[i] * xci[i+1] ;
		som2 += xci[i+1] ;
		xco[i] = cx3[i]*(som1 - cx2[i] * som2) ;
	    }	// Fin de la boucle sur r
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}



// cas R_JACO02
//-----------
void _d2sdx2_r_jaco02(Tbl *tb, int & )
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double* tb1 = new double[nr] ;
		for (int m =0 ; m<nr ; m++) { 
			tb1[m]=xci[m]; 
		}
	    double* res = new double[nr] ;
	    _d2sdx2_1d_r_jaco02(nr,tb1,res) ;
	    for (int i = 0 ; i<nr ; i++ ) {
		xco[i] = res[i] ;
		}
	    delete [] res ;
	    delete [] tb1 ;		
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}


}
