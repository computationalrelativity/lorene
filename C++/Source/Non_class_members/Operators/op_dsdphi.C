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


char op_dsdphi_C[] = "$Header$" ;

/*
 * Ensemble des routines de base pour la derivation par rapport a phi
 * (Utilisation interne)
 * 
 *	void _dsdphi_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.8  2000/11/14  15:08:16  eric
 * Traitement du cas np=1 dans P_COSSIN_I.
 *
 * Revision 2.7  2000/10/04  15:54:09  eric
 * *** empty log message ***
 *
 * Revision 2.6  2000/10/04  12:50:38  eric
 * Ajout de la base P_COSSIN_I.
 *
 * Revision 2.5  2000/10/04  11:50:32  eric
 * Ajout d' abort() dans le cas non prevu.
 *
 * Revision 2.4  2000/08/18  13:20:01  eric
 * Traitement du cas np = 1.
 *
 * Revision 2.3  2000/02/24  16:41:03  eric
 *  Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 2.2  1999/11/15  16:38:03  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.1  1999/02/23  11:36:34  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/22  17:24:59  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tbl.h"
#include "proto.h"


// Routine pour les cas non prevus
//--------------------------------
void _dsdphi_pas_prevu(Tbl* , int & b) {
    cout << "Unknown phi basis in Mtbl_cf::dsdp() !" << endl ;
    cout << " basis: " << hex << b << endl ;
    abort() ; 
}

			    //------------------//
			    //  cas P_COSSIN    //
			    //------------------//
			    
void _dsdphi_p_cossin(Tbl* tb, int & )
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
    static int np_pre =0 ;

    // Test sur np pour initialisation eventuelle
    #pragma critical (loch_dsdphi_p_cossin)
    {
    if (np > np_pre) {
	np_pre = np ;
	cx = (double*)(realloc(cx, (np+2) * sizeof(double))) ;
	for (int i=0 ; i<np+2 ; i += 2) {
	    cx[i] = - (i/2) ;
	    cx[i+1] = (i/2) ;
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
    int k, j ;
    
    // k = 0: inutile, deviendra k=1
    xci += nr*nt ;  // Positionnement
    xco += nr*nt ;

    // k = 1: 0
    for (j=0 ; j<nr*nt ; j++) {
	xco[j] = 0 ;
    }   // Fin de la boucle sur r * theta
    xci += nr*nt ;  // Positionnement
    xco += nr*nt ;
    
    // 2 =< k <= np-1
    for (k=2 ; k<np ; k++) {
	for (j=0 ; j<nr*nt ; j++) {
	    xco[j] = cx[k] * xci[j] ;
	}   // Fin de la boucle sur r * theta
	xci += nr*nt ;	// Positionnement
	xco += nr*nt ;
    }	// Fin de la boucle sur phi

    // k = np: inutile, deviendra k=np+1
    xci += nr*nt ;  // Positionnement
    xco += nr*nt ;

    // k = np+1
    for (j=0 ; j<nr*nt ; j++) {
	xco[j] = 0 ;
    }   // Fin de la boucle sur r * theta
    
    // Inversion des sinus et cosinus (appel a dswap de BLAS)
    int nbr = nr*nt ;
    int inc = 1 ;
    double* p1 = xo ;
    double* p2 = p1 + nr*nt ;
    for (k=0 ; k<np+1 ; k +=2, p1 += 2*nr*nt, p2 += 2*nr*nt) {
	dswap_(&nbr, p1, &inc, p2, &inc) ;
    }

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

			    //------------------//
			    //  cas P_COSSIN_P  //
			    //------------------//
			    
void _dsdphi_p_cossin_p(Tbl* tb, int & )
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
    
    // Cas particulier de la symetrie axiale : 
    if (np == 1) {
	tb->set_etat_zero() ; 
	return ; 
    }
    
    // Variables statiques
    static double* cx = 0 ;
    static int np_pre =0 ;

    // Test sur np pour initialisation eventuelle
    #pragma critical (loch_dsdphi_p_cossin_p)
    {
    if (np > np_pre) {
	np_pre = np ;
	cx = (double*)(realloc(cx, (np+2) * sizeof(double))) ;
	int i ;
	for (i=0 ; i<np+2 ; i += 2) {
	    cx[i] = - (i/2) ;
	    cx[i+1] = (i/2) ;
	    }
	for (i=0 ; i<np+2 ; i++) {
	    cx[i] = 2 * cx[i] ;
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
    int k, j ;
    
    // k = 0: inutile, deviendra k=1
    xci += nr*nt ;  // Positionnement
    xco += nr*nt ;

    // k = 1: 0
    for (j=0 ; j<nr*nt ; j++) {
	xco[j] = 0 ;
    }   // Fin de la boucle sur r * theta
    xci += nr*nt ;  // Positionnement
    xco += nr*nt ;
    
    // 2 =< k <= np-1
    for (k=2 ; k<np ; k++) {
	for (j=0 ; j<nr*nt ; j++) {
	    xco[j] = cx[k] * xci[j] ;
	}   // Fin de la boucle sur r * theta
	xci += nr*nt ;	// Positionnement
	xco += nr*nt ;
    }	// Fin de la boucle sur phi

    // k = np: inutile, deviendra k=np+1
    xci += nr*nt ;  // Positionnement
    xco += nr*nt ;

    // k = np+1
    for (j=0 ; j<nr*nt ; j++) {
	xco[j] = 0 ;
    }   // Fin de la boucle sur r * theta
    
    // Inversion des sinus et cosinus (appel a dswap de BLAS)
    int nbr = nr*nt ;
    int inc = 1 ;
    double* p1 = xo ;
    double* p2 = p1 + nr*nt ;
    for (k=0 ; k<np+1 ; k +=2, p1 += 2*nr*nt, p2 += 2*nr*nt) {
	dswap_(&nbr, p1, &inc, p2, &inc) ;
    }

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

			    //------------------//
			    //  cas P_COSSIN_I  //
			    //------------------//
			    
void _dsdphi_p_cossin_i(Tbl* tb, int & )
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

    int ntnr = nt*nr ;		    // increment en phi
    
    // Cas particulier de la symetrie axiale : 
    if (np == 1) {
	tb->set_etat_zero() ; 
	return ; 
    }
            
    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // pt. sur le tableau de double entree
    const double* xi = tb->t ;
    
    // k = 0 :  resultat sur cos(phi)
    // ------------------------------
    
    double* xco = xo ;			    // coef de cos(phi)
    const double* xci = xi + 2*ntnr ;	    // coef de sin(phi)
    for (int i=0; i<ntnr; i++) {
	xco[i] = xci[i] ; 
    }

    // k = 1 :  mise a zero
    // --------------------
    
    xco = xo + ntnr ;		
    for (int i=0; i<ntnr; i++) {
	xco[i] = 0 ; 
    }

    // k = 2 :  resultat sur sin(phi)
    // ------------------------------
    
    xco = xo + 2*ntnr ;		    // coef de sin(phi)
    xci = xi ;			    // coef de cos(phi)
    for (int i=0; i<ntnr; i++) {
	xco[i] = - xci[i] ; 
    }

    // k >= 3
    // ------
    
    for (int k=3; k<np; k+=2) {

	// resultat sur cos(k phi) 
	// -----------------------
	
	xco = xo + k*ntnr ;	    // coef de cos(k phi) 
	xci = xi + (k+1)*ntnr ;     // coef de sin(k phi) 

	for (int i=0; i<ntnr; i++) {
	    xco[i] = k * xci[i] ; 
	}
	
	// resultat sur sin(k phi) 
	// -----------------------
	
	xco = xo + (k+1)*ntnr ;	    // coef de sin(k phi) 
	xci = xi + k*ntnr ;	    // coef de cos(k phi) 

	for (int i=0; i<ntnr; i++) {
	    xco[i] = - k * xci[i] ; 
	}
	
    }

    // k = np+1 :  mise a zero
    // -----------------------
    
    xco = xo + (np+1)*ntnr ;		
    for (int i=0; i<ntnr; i++) {
	xco[i] = 0 ; 
    }

    // On remet les choses la ou il faut
    // ---------------------------------
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}
