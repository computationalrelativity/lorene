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


char op_dsdx_C[] = "$Header$" ;

/* 
 * Ensemble des routines de base de derivation par rapport a r
 * (Utilisation interne)
 * 
 *	void _dsdx_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.5  2007/12/20 09:11:08  jl_cornou
 * Correction of an error in op_sxpun about Jacobi(0,2) polynomials
 *
 * Revision 1.4  2007/12/11 15:28:18  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.3  2006/05/17 13:15:18  j_novak
 * Added a test for the pure angular grid case (nr=1), in the shell.
 *
 * Revision 1.2  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2000/09/07  12:49:04  phil
 * *** empty log message ***
 *
 * Revision 2.5  2000/02/24  16:41:25  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 2.4  1999/11/15  16:38:26  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.3  1999/09/13  11:31:49  phil
 * gestion des cas k==1 et k==np+1\
 *
 * Revision 2.2  1999/02/22  17:25:15  hyc
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/22  16:17:36  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

// Fichier includes
#include "tbl.h"
#include <math.h>

// Routine pour les cas non prevus
//--------------------------------
void _dsdx_pas_prevu(Tbl* , int & b) {
    cout << "dsdx pas prevu..." << endl ;
    cout << " base: " << b << endl ;
    abort () ;
}

// cas R_CHEB
//-----------
void _dsdx_r_cheb(Tbl *tb, int & )
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
    if (nr > 2) { // If not an angular grid...
	// On y va...
	double* xi = tb->t ;
	double* xci = xi ;	// Pointeurs
	double* xco = xo ;	//  courants
	
	int borne_phi = np + 1 ; 
	if (np == 1) borne_phi = 1 ; 
	
	for (int k=0 ; k< borne_phi ; k++)
	    // On evite le coefficient de sin(0*phi)
	    if (k==1) {
		xci += nr*nt ;
		xco += nr*nt ;
	    }
	    else {
		for (int j=0 ; j<nt ; j++) {
		    
		    double som ;
		    
		    xco[nr-1] = 0 ;
		    som = 2*(nr-1) * xci[nr-1] ;
		    xco[nr-2] = som ;
		    for (int i = nr-4 ; i >= 0 ; i -= 2 ) {
			som += 2*(i+1) * xci[i+1] ;
			xco[i] = som ;
		    }	// Fin de la premiere boucle sur r
		    som = 2*(nr-2) * xci[nr-2] ;
		    xco[nr-3] = som ;
		    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
			som += 2*(i+1) * xci[i+1] ;
			xco[i] = som ;
		    }	// Fin de la deuxieme boucle sur r
		    xco[0] *= .5 ;
		    
		    xci += nr ;
		    xco += nr ;
		}   // Fin de la boucle sur theta
	    }	// Fin de la boucle sur phi
    }
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

// cas R_CHEBU
//------------
void _dsdx_r_chebu(Tbl *tb, int & )
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
    
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    for (int k=0 ; k< borne_phi ; k++)
    
	// On evite le coefficient de sin(0*phi)
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	}
	
	else { 
	for (int j=0 ; j<nt ; j++) {

	    double som ;
	    
	    xco[nr-1] = 0 ;
	    som = 2*(nr-1) * xci[nr-1] ;
	    xco[nr-2] = som ;
	    for (int i = nr-4 ; i >= 0 ; i -= 2 ) {
		som += 2*(i+1) * xci[i+1] ;
		xco[i] = som ;
	    }	// Fin de la premiere boucle sur r
	    som = 2*(nr-2) * xci[nr-2] ;
	    xco[nr-3] = som ;
	    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
		som += 2*(i+1) * xci[i+1] ;
		xco[i] = som ;
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
void _dsdx_r_chebp(Tbl *tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = b & MSQ_T ;
	int base_p = b & MSQ_P ;
	b = base_p | base_t | R_CHEBI ;
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
    
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    for (int k=0 ; k< borne_phi ; k++)
	
	// On evite le coefficient de sin(0*phi)
	
	if (k==1) {
	    xco += nr*nt ;
	    xci += nr*nt ;
	    }
	    
	     
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som ;
	    
	    xco[nr-1] = 0 ;
	    som = 4*(nr-1) * xci[nr-1] ;
	    xco[nr-2] = som ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		som += 4*(i+1) * xci[i+1] ;
		xco[i] = som ;
	    }	// Fin de la boucle sur r
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // pair -> impair
    int base_t = b & MSQ_T ;
    int base_p = b & MSQ_P ;
    b = base_p | base_t | R_CHEBI ;
}

// cas R_CHEBI
//------------
void _dsdx_r_chebi(Tbl *tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = b & MSQ_T ;
	int base_p = b & MSQ_P ;
	b = base_p | base_t | R_CHEBP ;
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
    
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    for (int k=0 ; k< borne_phi ; k++)
	// On evite le coefficient de sin(0*phi)
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	}
	
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som ;
	    
	    xco[nr-1] = 0 ;
	    som = 2*(2*nr-3) * xci[nr-2] ;
	    xco[nr-2] = som ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		som += 2*(2*i+1) * xci[i] ;
		xco[i] = som ;
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
    // impair -> pair
    int base_t = b & MSQ_T ;
    int base_p = b & MSQ_P ;
    b = base_p | base_t | R_CHEBP ;
}

// cas R_CHEBPIM_P
//----------------
void _dsdx_r_chebpim_p(Tbl *tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = b & MSQ_T ;
	int base_p = b & MSQ_P ;
	b = base_p | base_t | R_CHEBPIM_I ;
	return ;
    }
    
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
    double* xci ;	// Pointeurs
    double* xco ;	//  courants
    
    int auxiliaire ;
    // Partie paire
    xci = xi ;
    xco = xo ;
    for (int k=0 ; k<np+1 ; k += 4) {
	auxiliaire = (k==np) ? 1 : 2 ; // evite de prendre le dernier coef
	for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
	    
	    
		// On evite le sin (0*phi)
	
	if ((k==0) && (kmod == 1)) { 
		xco += nr*nt ;
		xci += nr*nt ;
		}
		
	else 
	    
	    for (int j=0 ; j<nt ; j++) {

		double som ;
	    
		xco[nr-1] = 0 ;
		som = 4*(nr-1) * xci[nr-1] ;
		xco[nr-2] = som ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    som += 4*(i+1) * xci[i+1] ;
		    xco[i] = som ;
		}	// Fin de la boucle sur r
	    
		xci += nr ; // au
		xco += nr ; //	suivant
	    }   // Fin de la boucle sur theta
	}   // Fin de la boucle sur kmod
	xci += 2*nr*nt ;    // saute
	xco += 2*nr*nt ;    //	2 phis
    }	// Fin de la boucle sur phi

    // Partie impaire
    xci = xi + 2*nr*nt ;
    xco = xo + 2*nr*nt ;
    for (int k=2 ; k<np+1 ; k += 4) {
	auxiliaire = (k==np) ? 1 : 2 ; // evite de prendre le dernier coef
	for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
	    for (int j=0 ; j<nt ; j++) {

		double som ;

		xco[nr-1] = 0 ;
		som = 2*(2*nr-3) * xci[nr-2] ;
		xco[nr-2] = som ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    som += 2*(2*i+1) * xci[i] ;
		    xco[i] = som ;
		}	// Fin de la boucle sur r
		xco[0] *= .5 ;
	    
		xci += nr ;
		xco += nr ;
	}   // Fin de la boucle sur theta
	}   // Fin de la boucle sur kmod
	xci += 2*nr*nt ;
	xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // (pair,impair) -> (impair,pair)
    int base_t = b & MSQ_T ;
    int base_p = b & MSQ_P ;
    b = base_p | base_t | R_CHEBPIM_I ;
}

// cas R_CHEBPIM_I
//----------------
void _dsdx_r_chebpim_i(Tbl *tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = b & MSQ_T ;
	int base_p = b & MSQ_P ;
	b = base_p | base_t | R_CHEBPIM_P ;
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort
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
    double* xci ;	// Pointeurs
    double* xco ;	//  courants
    int auxiliaire ;
    
    // Partie impaire
    xci = xi ;
    xco = xo ;
    for (int k=0 ; k<np+1 ; k += 4) {
	auxiliaire = (k==np) ? 1 : 2 ; // evite de prendre le dernier coef
	for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
	      
	    
		// On evite le sin (0*phi)
	
	if ((k==0) && (kmod == 1)) { 
		xco += nr*nt ;
		xci += nr*nt ;
		}
		
	else   
	    
	    for (int j=0 ; j<nt ; j++) {

		double som ;
	
		xco[nr-1] = 0 ;
		som = 2*(2*nr-3) * xci[nr-2] ;
		xco[nr-2] = som ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    som += 2*(2*i+1) * xci[i] ;
		    xco[i] = som ;
		}	// Fin de la boucle sur r
		xco[0] *= .5 ;
	    
		xci += nr ;
		xco += nr ;
	    }   // Fin de la boucle sur theta
	}   // Fin de la boucle sur kmod
	xci += 2*nr*nt ;
	xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi

    // Partie paire
    xci = xi + 2*nr*nt ;
    xco = xo + 2*nr*nt ;
    for (int k=2 ; k<np+1 ; k += 4) {
    	auxiliaire = (k==np) ? 1 : 2 ; // evite de prendre le dernier coef
	for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
	    for (int j=0 ; j<nt ; j++) {

		double som ;
		
		xco[nr-1] = 0 ;
		som = 4*(nr-1) * xci[nr-1] ;
		xco[nr-2] = som ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    som += 4*(i+1) * xci[i+1] ;
		    xco[i] = som ;
		}	// Fin de la boucle sur r
	    
		xci += nr ;
		xco += nr ;
	}   // Fin de la boucle sur theta
	}   // Fin de la boucle sur kmod
	xci += 2*nr*nt ;
	xco += 2*nr*nt ;
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // (impair,pair) -> (pair,impair)
    int base_t = b & MSQ_T ;
    int base_p = b & MSQ_P ;
    b = base_p | base_t | R_CHEBPIM_P ;
}

// cas R_CHEBPI_P
//------------
void _dsdx_r_chebpi_p(Tbl *tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = b & MSQ_T ;
	int base_p = b & MSQ_P ;
	b = base_p | base_t | R_CHEBPI_I ;
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
    
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    for (int k=0 ; k< borne_phi ; k++)
	
	// On evite le coefficient de sin(0*phi)
	
	if (k==1) {
	    xco += nr*nt ;
	    xci += nr*nt ;
	    }
	    
	     
	else {
	for (int j=0 ; j<nt ; j++) {
	    int l = j%2 ;
	    double som ;
	    
	    if(l==0){
	      xco[nr-1] = 0 ;
	      som = 4*(nr-1) * xci[nr-1] ;
	      xco[nr-2] = som ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		som += 4*(i+1) * xci[i+1] ;
		xco[i] = som ;
	      }	// Fin de la boucle sur r
	    } else {
	      xco[nr-1] = 0 ;
	      som = 2*(2*nr-3) * xci[nr-2] ;
	      xco[nr-2] = som ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		som += 2*(2*i+1) * xci[i] ;
		xco[i] = som ;
	      }	// Fin de la boucle sur r
	      xco[0] *= .5 ;
	    }
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // pair -> impair
    int base_t = b & MSQ_T ;
    int base_p = b & MSQ_P ;
    b = base_p | base_t | R_CHEBPI_I ;
}

// cas R_CHEBPI_I
//------------
void _dsdx_r_chebpi_i(Tbl *tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = b & MSQ_T ;
	int base_p = b & MSQ_P ;
	b = base_p | base_t | R_CHEBPI_P ;
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
    
    int borne_phi = np + 1 ; 
    if (np == 1) borne_phi = 1 ; 
    
    for (int k=0 ; k< borne_phi ; k++)
	// On evite le coefficient de sin(0*phi)
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	}
	
	else {
	for (int j=0 ; j<nt ; j++) {
	    int l = j%2 ;
	    double som ;
	    
	    if(l==1){
	      xco[nr-1] = 0 ;
	      som = 4*(nr-1) * xci[nr-1] ;
	      xco[nr-2] = som ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		som += 4*(i+1) * xci[i+1] ;
		xco[i] = som ;
	      }	// Fin de la boucle sur r
	    } else {
	      xco[nr-1] = 0 ;
	      som = 2*(2*nr-3) * xci[nr-2] ;
	      xco[nr-2] = som ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		som += 2*(2*i+1) * xci[i] ;
		xco[i] = som ;
	      }	// Fin de la boucle sur r
	      xco[0] *= .5 ;
	    }
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // impair -> pair
    int base_t = b & MSQ_T ;
    int base_p = b & MSQ_P ;
    b = base_p | base_t | R_CHEBPI_P ;
}


// cas R_JACO02
//-----------
void _dsdx_r_jaco02(Tbl *tb, int & )
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
    if (nr > 2) { // If not an angular grid...
	// On y va...
	double* xi = tb->t ;
	double* xci = xi ;	// Pointeurs
	double* xco = xo ;	//  courants
	
	int borne_phi = np + 1 ; 
	if (np == 1) borne_phi = 1 ; 
	
	for (int k=0 ; k< borne_phi ; k++)
	    // On evite le coefficient de sin(0*phi)
	    if (k==1) {
		xci += nr*nt ;
		xco += nr*nt ;
	    }
	    else {
		for (int j=0 ; j<nt ; j++) {
		    
		    double som ;
		    xco[nr-1] = 0 ;
		 
		    for (int i = 0 ; i < nr-1 ; i++ ) {
		    
		      som = 0 ;
		      for (int m = i+1 ; m < nr ; m++ ) {
			som += (1-pow((-1),(m-i))*(i+1)*(i+2)/double((m+1)*(m+2)))* xci[m] ;
		      } // Fin de la boucle annexe
			xco[i] = (i+1.5)*som ;
		    }// Fin de la boucle sur r
		    		    
		    xci += nr ;
		    xco += nr ;
		}   // Fin de la boucle sur theta
	    }	// Fin de la boucle sur phi
    }
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}
