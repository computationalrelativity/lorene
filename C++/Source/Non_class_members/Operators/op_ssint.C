/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 2004 Michael Forot
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


char op_ssint_C[] = "$Header$" ;

/* 
 * Ensemble des routines de base de division par rapport a sint theta
 * (Utilisation interne)
 * 
 *	void _ssint_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */
 
/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/02/24  16:42:49  eric
 *  Initialisation a zero du tableau xo avant le calcul.
 *
 *
 * $Header$
 *
 */
 
// Fichier includes
#include "tbl.h"

 
		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _ssint_pas_prevu(Tbl * tb, int& base) {
    cout << "ssint pas prevu..." << endl ;
    cout << "Tbl: " << tb << " base: " << base << endl ;
    abort () ;
    exit(-1) ;
}

		    //----------------
		    // cas T_COS_P
		    //----------------
		    
void _ssint_t_cos_p(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
        
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt) ;
	xco += nr * (nt-1) ;
	
	// Initialisation de som 
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;	// mise a zero du dernier coef en theta
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] -= 2*xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	xci -= nr ;
	xci += nr*nt ;
	xco += nr*nt ;

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
    
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt) ;
	xco += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ; // mise a zero dui dernier coefficient en theta.
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] -= 2*xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	xci -= nr ;
	xci += nr*nt ;
	xco += nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
  
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_SIN_I ;
}
			
			//--------------
			// cas T_SIN_P
			//--------------
			
void _ssint_t_sin_p(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
    
    // Dernier j: j = nt-1
    // Positionnement
    xci += nr * (nt-1) ;
    xco += nr * (nt-1) ;
	
    // Initialisation de som
    for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
    xco -= nr ;
    for (int i=0 ; i<nr ; i++) xco[i] = 0 ;
	
    // j suivants
    for (int j=nt-2 ; j >= 1 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	xci -= nr ;
	xci += nr*nt ;
	xco += nr*nt ;

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
    
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt-1) ;
	xco += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
	xco -= nr ;
	for (int i=0 ; i<nr ; i++) xco[i] = 0 ;
	
	// j suivants
	for (int j=nt-2 ; j >= 1 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	xci -= nr ;
	xci += nr*nt ;
	xco += nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COS_I ;
   
}
			//--------------
			// cas T_SIN_I
			//--------------
			
void _ssint_t_sin_i(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
    
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt-1) ;
	xco += nr * (nt-1) ;

	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    xco[i] = 0. ;
	    som[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Normalisation du premier theta
	for (int i=0 ; i<nr ; i++) {
	    xco[i] *= .5 ;
	}
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr*nt ;

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
    
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt-1) ;
	xco += nr * (nt-1) ;

	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    xco[i] = 0. ;
	    som[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Normalisation du premier theta
	for (int i=0 ; i<nr ; i++) {
	    xco[i] *= .5 ;
	}
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr*nt ;
    }	// Fin de la boucle sur phi


    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COS_P ;

}
		    //---------------------
		    // cas T_COS_I
		    //---------------------
void _ssint_t_cos_i(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
    
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt) ;
	xco += nr * (nt-1) ;
	
	// Initialisation de som 
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;	// mise a zero du dernier coef en theta
	}
	
	// j suivants
	for (int j=nt-1 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] -= 2*xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	    xco -= nr ;
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr ;
	xco += nr*nt ;

    // k = 1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
    
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt) ;
	xco += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;	// mise a zero du dernier coef en theta
	}
	
	// j suivants
	for (int j=nt-1 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] -= 2*xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	    xco -= nr ;
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr ;
	xco += nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_SIN_P ;
  
}
			//----------------------
			// cas T_COSSIN_CP
			//----------------------
void _ssint_t_cossin_cp(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
    
    double cx ;

    // k = 0
    int m = 0 ;	    // Parite de l'harmonique en phi
	
	// Dernier j: j = nt-1
	// Positionnement
	    xci += nr * (nt) ;
	    cx = -2 ;
	xco += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += cx * xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	if (m == 0) {
	    xci -= nr ;
	    }
	xci += nr*nt ;
	xco += nr*nt ;

    // k=1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	m = (k/2) % 2 ;	    // Parite de l'harmonique en phi
	
	// Dernier j: j = nt-1
	// Positionnement
	if (m == 0) {
	    xci += nr * (nt) ;
	    cx = -2 ;
	    }
	else {
	    xci += nr * (nt-1) ;
	    cx = 2 ;
	}
	xco += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += cx * xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Normalisation du premier theta dans le cas sin(impair)
	if (m == 1) {
	    for (int i=0 ; i<nr ; i++) {
		xco[i] *= .5 ;
	    }
	}
	// Positionnement phi suivant
	if (m == 0) {
	    xci -= nr ;
	    }
	xci += nr*nt ;
	xco += nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_SI ;
}
			
			//------------------
			// cas T_COSSIN_CI
			//----------------
void _ssint_t_cossin_ci(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
    int m = 0 ;	    // Parite de l'harmonique en phi
	
	
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt) ;
	xco += nr * (nt-1) ;
	
	// Initialisation de som 
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;	// mise a zero du dernier coef en theta
	}
	
	// j suivants
	for (int j=nt-1 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] -= 2*xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	    xco -= nr ;
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr ;
	xco += nr*nt ;

    // k=1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	m = (k/2) % 2 ;	    // Parite de l'harmonique en phi
	
	switch(m) {
	    case 0:	    // m pair: cos(impair)
		// Dernier j: j = nt-1
		// Positionnement
		xci += nr * (nt) ;
		xco += nr * (nt-1) ;
	
		// Initialisation de som 
		for (int i=0 ; i<nr ; i++) {
		    som[i] = 0. ;
		    xco[i] = 0. ;	// mise a zero du dernier coef en theta
		}
	
		// j suivants
		for (int j=nt-1 ; j >= 0 ; j--) {
		    // Positionnement
		    xci -= nr ;
		    // Calcul
		    for (int i=0 ; i<nr ; i++ ) {
			som[i] -= 2*xci[i] ;
			xco[i] = som[i] ;
		    }	// Fin de la boucle sur r
		    xco -= nr ;
		}   // Fin de la boucle sur theta
		// Positionnement phi suivant
		xci += nr*nt ;
		xco += nr ;
		xco += nr*nt ;
	    break ;
	    
	    case 1:	    // m impair: sin(impair)
		// Dernier j: j = nt-1
		// Positionnement
		xci += nr * (nt-1) ;
		xco += nr * (nt-1) ;
	
		// Initialisation de som
		for (int i=0 ; i<nr ; i++) {
		    som[i] = 0. ;
		    xco[i] = 0. ;
		}
		xco -= nr ;
		for (int i=0 ; i<nr ; i++) xco[i] = 0 ;
	
		// j suivants
		for (int j=nt-2 ; j >= 1 ; j--) {
		    // Positionnement
		    xci -= nr ;
		    xco -= nr ;
		    // Calcul
		    for (int i=0 ; i<nr ; i++ ) {
			som[i] += 2*xci[i] ;
			xco[i] = som[i] ;
		    }	// Fin de la boucle sur r
		}   // Fin de la boucle sur theta
		// Positionnement phi suivant
		xci -= nr ;
		xci += nr*nt ;
		xco += nr*nt ;
	    break ;
	}   // Fin du switch
    }   // Fin de la boucle en phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_SP ;
}

			//---------------------	
			// cas T_COSSIN_SI
			//----------------
void _ssint_t_cossin_si(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
    
    double cx ;

    // k = 0
    int m = 0 ;	    // Parite de l'harmonique en phi
	
		
	// Dernier j: j = nt-1
	// Positionnement
	    xci += nr * (nt-1) ;
	    cx = 2 ;
	xco += nr * (nt-1) ;
	

	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += cx * xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Normalisation du premier theta dans le cas sin(impair)
	for (int i=0 ; i<nr ; i++) {
		xco[i] *= .5 ;
	}

	xci += nr*nt ;
	xco += nr*nt ;

    // k=1 
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	m = (k/2) % 2 ;	    // Parite de l'harmonique en phi
	
	// Dernier j: j = nt-1
	// Positionnement
	if (m == 1) {
	    xci += nr * (nt) ;
	    cx = -2 ;
	    }
	else {
	    xci += nr * (nt-1) ;
	    cx = 2 ;
	}
	xco += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += cx * xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Normalisation du premier theta dans le cas sin(impair)
	if (m == 0) {
	    for (int i=0 ; i<nr ; i++) {
		xco[i] *= .5 ;
	    }
	}
	// Positionnement phi suivant
	if (m == 1) {
	    xci -= nr ;
	    }
	xci += nr*nt ;
	xco += nr*nt ;
    }	// Fin de la boucle sur phi

    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_CP ;

    
}
			//---------------------	
			// cas T_COSSIN_SP
			//----------------
void _ssint_t_cossin_sp(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
    
    double cx ;

    // k = 0
    int m = 0 ;	    // Parite de l'harmonique en phi
	
		
	// Dernier j: j = nt-1
	// Positionnement
	    xci += nr * (nt-1) ;
	    cx = 2 ;
	    xco += nr * (nt-1) ;
	

	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
	xco -= nr ;
	for (int i=0 ; i<nr ; i++) xco[i] = 0 ;
	
	// j suivants
	for (int j=nt-2 ; j >= 1 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += cx * xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	xci += nr*nt ;
	xci -= nr ;
	xco += nr*nt ;

    // k=1 
    xci += nr*nt ;
    xco += nr*nt ;
    
    for (int k=2 ; k<np+1 ; k++) {
	m = (k/2) % 2 ;	    // Parite de l'harmonique en phi
	
	switch(m) {
	    case 1:	    // m impair: cos(impair)
		// Dernier j: j = nt-1
		// Positionnement
		xci += nr * (nt) ;
		xco += nr * (nt-1) ;
	
		// Initialisation de som 
		for (int i=0 ; i<nr ; i++) {
		    som[i] = 0. ;
		    xco[i] = 0. ;	// mise a zero du dernier coef en theta
		}
	
		// j suivants
		for (int j=nt-1 ; j >= 0 ; j--) {
		    // Positionnement
		    xci -= nr ;
		    // Calcul
		    for (int i=0 ; i<nr ; i++ ) {
			som[i] -= 2*xci[i] ;
			xco[i] = som[i] ;
		    }	// Fin de la boucle sur r
		    xco -= nr ;
		}   // Fin de la boucle sur theta
		// Positionnement phi suivant
		xci += nr*nt ;
		xco += nr ;
		xco += nr*nt ;
	    break ;
	    
	    case 0:	    // m pair: sin(pair)
		// Dernier j: j = nt-1
		// Positionnement
		xci += nr * (nt-1) ;
		xco += nr * (nt-1) ;
	
		// Initialisation de som
		for (int i=0 ; i<nr ; i++) {
		    som[i] = 0. ;
		    xco[i] = 0. ;
		}
		xco -= nr ;
		for (int i=0 ; i<nr ; i++) xco[i] = 0 ;
	
		// j suivants
		for (int j=nt-2 ; j >= 1 ; j--) {
		    // Positionnement
		    xci -= nr ;
		    xco -= nr ;
		    // Calcul
		    for (int i=0 ; i<nr ; i++ ) {
			som[i] += 2*xci[i] ;
			xco[i] = som[i] ;
		    }	// Fin de la boucle sur r
		}   // Fin de la boucle sur theta
		// Positionnement phi suivant
		xci -= nr ;
		xci += nr*nt ;
		xco += nr*nt ;
	    break ;
	}   // Fin du switch
    }   // Fin de la boucle en phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_CI ;

    
}

			//----------------------
			// cas T_COSSIN_C
			//----------------------
void _ssint_t_cossin_c(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
    
    double cx ;

    // k = 0
    int m = 0 ;	    // Parite de l'harmonique en phi
	
	// Dernier j: j = nt-1
	// Positionnement
	    xci += nr * (nt) ;
	    cx = -2 ;
	xco += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += cx * xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	if (m == 0) {
	    xci -= nr ;
	    }
	xci += nr*nt ;
	xco += nr*nt ;

    // k=1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	m = (k/2) % 2 ;	    // Parite de l'harmonique en phi
	
	// Dernier j: j = nt-1
	// Positionnement
	if (m == 0) {
	    xci += nr * (nt) ;
	    cx = -2 ;
	    }
	else {
	    xci += nr * (nt-1) ;
	    cx = 2 ;
	}
	xco += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += cx * xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Normalisation du premier theta dans le cas sin(impair)
	if (m == 1) {
	    for (int i=0 ; i<nr ; i++) {
		xco[i] *= .5 ;
	    }
	}
	// Positionnement phi suivant
	if (m == 0) {
	    xci -= nr ;
	    }
	xci += nr*nt ;
	xco += nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_S ;
}

			//---------------------	
			// cas T_COSSIN_S
			//----------------
void _ssint_t_cossin_s(Tbl* tb, int & b)
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
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // zone de travail
    double* som = new double [nr] ;
    
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
    
    double cx ;

    // k = 0
    int m = 0 ;	    // Parite de l'harmonique en phi
	
		
	// Dernier j: j = nt-1
	// Positionnement
	    xci += nr * (nt-1) ;
	    cx = 2 ;
	    xco += nr * (nt-1) ;
	

	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    som[i] = 0. ;
	    xco[i] = 0. ;
	}
	xco -= nr ;
	for (int i=0 ; i<nr ; i++) xco[i] = 0 ;
	
	// j suivants
	for (int j=nt-2 ; j >= 1 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += cx * xci[i] ;
		xco[i] = som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	xci += nr*nt ;
	xci -= nr ;
	xco += nr*nt ;

    // k=1 
    xci += nr*nt ;
    xco += nr*nt ;
    
    for (int k=2 ; k<np+1 ; k++) {
	m = (k/2) % 2 ;	    // Parite de l'harmonique en phi
	
	switch(m) {
	    case 1:	    // m impair: cos(impair)
		// Dernier j: j = nt-1
		// Positionnement
		xci += nr * (nt) ;
		xco += nr * (nt-1) ;
	
		// Initialisation de som 
		for (int i=0 ; i<nr ; i++) {
		    som[i] = 0. ;
		    xco[i] = 0. ;	// mise a zero du dernier coef en theta
		}
	
		// j suivants
		for (int j=nt-1 ; j >= 0 ; j--) {
		    // Positionnement
		    xci -= nr ;
		    // Calcul
		    for (int i=0 ; i<nr ; i++ ) {
			som[i] -= 2*xci[i] ;
			xco[i] = som[i] ;
		    }	// Fin de la boucle sur r
		    xco -= nr ;
		}   // Fin de la boucle sur theta
		// Positionnement phi suivant
		xci += nr*nt ;
		xco += nr ;
		xco += nr*nt ;
	    break ;
	    
	    case 0:	    // m pair: sin(pair)
		// Dernier j: j = nt-1
		// Positionnement
		xci += nr * (nt-1) ;
		xco += nr * (nt-1) ;
	
		// Initialisation de som
		for (int i=0 ; i<nr ; i++) {
		    som[i] = 0. ;
		    xco[i] = 0. ;
		}
		xco -= nr ;
		for (int i=0 ; i<nr ; i++) xco[i] = 0 ;
	
		// j suivants
		for (int j=nt-2 ; j >= 1 ; j--) {
		    // Positionnement
		    xci -= nr ;
		    xco -= nr ;
		    // Calcul
		    for (int i=0 ; i<nr ; i++ ) {
			som[i] += 2*xci[i] ;
			xco[i] = som[i] ;
		    }	// Fin de la boucle sur r
		}   // Fin de la boucle sur theta
		// Positionnement phi suivant
		xci -= nr ;
		xci += nr*nt ;
		xco += nr*nt ;
	    break ;
	}   // Fin du switch
    }   // Fin de la boucle en phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] som ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    b = base_r | base_p | T_COSSIN_C ;

    
}
