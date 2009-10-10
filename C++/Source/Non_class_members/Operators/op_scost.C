/*
 *   Copyright (c) 1999-2001 Jerome Novak
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


char op_scost_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2009/10/10 18:28:11  j_novak
 * New bases T_COS and T_SIN.
 *
 * Revision 1.5  2007/12/21 10:43:38  j_novak
 * Corrected some bugs in the case nt=1 (1 point in theta).
 *
 * Revision 1.4  2007/10/05 12:37:20  j_novak
 * Corrected a few errors in the theta-nonsymmetric case (bases T_COSSIN_C and
 * T_COSSIN_S).
 *
 * Revision 1.3  2005/02/16 15:29:40  m_forot
 * Correct T_COSSIN_S and T_COSSIN_C cases
 *
 * Revision 1.2  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2000/02/24  16:42:36  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 1.2  1999/11/22  14:34:16  novak
 * Suppression des variables inutiles
 *
 * Revision 1.1  1999/11/16  13:38:00  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

/* 
 * Ensemble des routines de base de division par rapport a cost theta
 * (Utilisation interne)
 * 
 *	void _scost_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */
 
 
// Fichier includes
#include "tbl.h"

 
		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _scost_pas_prevu(Tbl * tb, int& base) {
    cout << "scost pas prevu..." << endl ;
    cout << "Tbl: " << tb << " base: " << base << endl ;
    abort () ;
    exit(-1) ;
}

		    //--------------
		    // cas T_COS
		    //--------------
		    
void _scost_t_cos(Tbl* tb, int & b)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	switch(base_r){
	    case(R_CHEBPI_P):
		b = R_CHEBPI_I | base_p | T_COS ;
		break ;
	    case(R_CHEBPI_I):
		b = R_CHEBPI_P | base_p | T_COS ;
		break ;  
	    default:
		b = base_r | base_p | T_COS ;
		break;
	}
	return ;
    }
    
    // Protection
    assert(tb->get_etat() == ETATQCQ) ;
    
    // Pour le confort : nbre de points reels.
    int nr = (tb->dim).dim[0] ;
    int nt = (tb->dim).dim[1] ;
    int np = (tb->dim).dim[2] ;
    np = np - 2 ;
    
    // pt. sur le tableau de double resultat
    double* somP = new double [nr] ;
    double* somN = new double [nr] ;
    
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
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;	// mise a zero du dernier coef en theta
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
	      if((j%2) == 1) {
		somN[i] = -somN[i] ;
		somN[i] += 2*xci[i] ;
		xco[i] = somP[i] ;
	      }
	      else {
		somP[i] = -somP[i] ;
		somP[i] += 2*xci[i] ; 
		xco[i] = somN[i] ;
	      }
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// j=0 : le facteur 2...
	for (int i=0 ; i<nr ; i++) xco[i] *= .5 ;
 
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
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ; // mise a zero dui dernier coefficient en theta.
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		if((j%2) == 1) {
		    somN[i] = -somN[i];
		    somN[i] += 2 * xci[i] ;
		    xco[i] = somP[i] ;
		}
		else {
		    somP[i] = -somP[i];
		    somP[i] += 2 * xci[i] ;
		    xco[i] = somN[i];
		}
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// j=0 : facteur 2 ...
	for (int i=0 ; i<nr ; i++) xco[i] *= 0.5 ;

	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr*nt ;
    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
  
    //Menage
    delete [] somP ;
    delete [] somN ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
	switch(base_r){
	    case(R_CHEBPI_P):
		b = R_CHEBPI_I | base_p | T_COS ;
		break ;
	    case(R_CHEBPI_I):
		b = R_CHEBPI_P | base_p | T_COS ;
		break ;  
	    default:
		b = base_r | base_p | T_COS ;
		break;
	}
}
			
			//------------
			// cas T_SIN
			//------------
			
void _scost_t_sin(Tbl* tb, int & b)
{


    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	switch(base_r){
	    case(R_CHEBPI_P):
		b = R_CHEBPI_I | base_p | T_SIN ;
		break ;
	    case(R_CHEBPI_I):
		b = R_CHEBPI_P | base_p | T_SIN ;
		break ;  
	    default:
		b = base_r | base_p | T_SIN ;
		break;
	}
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
    double* somP = new double [nr] ;
    double* somN = new double [nr] ;
    
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
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;
	}
	
    // j suivants
    for (int j=nt-2 ; j >= 0 ; j--) {
      // Positionnement
      xci -= nr ;
      xco -= nr ;
      // Calcul
      for (int i=0 ; i<nr ; i++ ) {
	if((j%2) == 1) {
	  somN[i] = -somN[i] ;
	  somN[i] += 2*xci[i] ;
	  xco[i] = somP[i] ;
	}
	else {
	  somP[i] = -somP[i] ;
	  somP[i] += 2*xci[i] ;		  
	  xco[i] = somN[i] ;
	}
      }	// Fin de la boucle sur r
    }   // Fin de la boucle sur theta
    // j=0 : sin(0*theta)
    for (int i=0 ; i<nr ; i++) xco[i] = 0. ;

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
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		if((j%2) == 1) {
		  somN[i] = -somN[i] ;
		  somN[i] += 2*xci[i] ;
		  xco[i] = somP[i] ;
	      }
	      else {
		  somP[i] = -somP[i] ;
		  somP[i] += 2*xci[i] ;		  
		  xco[i] = somN[i] ;
	      }
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// j=0 : sin(0*theta)
	for (int i=0 ; i<nr ; i++) xco[i] = 0. ;
	
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr*nt ;

    }	// Fin de la boucle sur phi

    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] somP ;
    delete [] somN ;

    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    switch(base_r){
	case(R_CHEBPI_P):
	    b = R_CHEBPI_I | base_p | T_SIN ;
	    break ;
	case(R_CHEBPI_I):
	    b = R_CHEBPI_P | base_p | T_SIN ;
	    break ;  
	default:
	    b = base_r | base_p | T_SIN ;
	    break;
    }
}
		    //----------------
		    // cas T_COS_P
		    //----------------
		    
void _scost_t_cos_p(Tbl* tb, int & b)
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
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = -som[i] ;
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
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = - som[i] ;
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
			// cas T_SIN_P
			//--------------
			
void _scost_t_sin_p(Tbl* tb, int & b)
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
    xci += nr * nt ;
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
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = -som[i] ;
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
	xci += nr * nt ;
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
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = -som[i] ;
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
			// cas T_SIN_I
			//--------------
			
void _scost_t_sin_i(Tbl* tb, int & b)
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
	xci += nr * (nt-1) ;
	xco += nr * (nt-1) ;

	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    xco[i] = 0. ;
	    som[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 1 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = -som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	if (nt > 1) {
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul pour le premier theta
	    for (int i=0 ; i<nr ; i++ ) {
		xco[i] =0. ;
	    }
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
	for (int j=nt-2 ; j >= 1 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = -som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	xci -= nr ;
	xco -= nr ;
	// Calcul pour le premier theta
	for (int i=0 ; i<nr ; i++ ) {
	  xco[i] =0. ;
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
    b = base_r | base_p | T_SIN_P ;

}
		    //---------------------
		    // cas T_COS_I
		    //---------------------
void _scost_t_cos_i(Tbl* tb, int & b)
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
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = - som[i] ;
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
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = -som[i] ;
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
			//----------------------
			// cas T_COSSIN_CP
			//----------------------
void _scost_t_cossin_cp(Tbl* tb, int & b)
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
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = -som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	xci -= nr ;
	xci += nr*nt ;
	xco += nr*nt ;

    // k=1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
      m = (k/2) % 2 ;	    // Parite de l'harmonique en phi
	
      switch(m) {
      case 0:	    // m pair: cos(pair)
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
	    som[i] += 2*xci[i] ;
	    xco[i] = som[i] ;
	    som[i] = -som[i] ;
	  }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Positionnement phi suivant
	xci -= nr ;
	xci += nr*nt ;
	xco += nr*nt ;
	break ;
	      
      case 1:	    // m impair: sin(impair)
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
	for (int j=nt-2 ; j >= 1 ; j--) {
	  // Positionnement
	  xci -= nr ;
	  xco -= nr ;
	  // Calcul
	  for (int i=0 ; i<nr ; i++ ) {
	    som[i] += 2*xci[i] ;
	    xco[i] = som[i] ;
	    som[i] = -som[i] ;
	  }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	xci -= nr ;
	xco -= nr ;
	// Calcul pour le premier theta
	for (int i=0 ; i<nr ; i++ ) {
	  xco[i] =0. ;
	}
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr*nt ;
	break;
      } // Fin du switch
    }	// Fin de la boucle sur phi

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
			
			//------------------
			// cas T_COSSIN_CI
			//----------------
void _scost_t_cossin_ci(Tbl* tb, int & b)
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

    // k = 0
    int m = 0 ;	    // Parite de l'harmonique en phi
	
	
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt-1) ;
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
		som[i] += 2*xci[i] ;
		xco[i] = som[i] ;
		som[i] = - som[i] ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Normalisation du premier theta
	for (int i=0 ; i<nr ; i++) {
	    xco[i] *= .5 ;
	}
	// Positionnement phi suivant
	xci += nr*nt ;
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
	xci += nr * (nt-1) ;
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
	    som[i] += 2*xci[i] ;
	    xco[i] = som[i] ;
	    som[i] = - som[i] ;
	  }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// Normalisation du premier theta
	for (int i=0 ; i<nr ; i++) {
	  xco[i] *= .5 ;
	}
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr*nt ;
	break ;
	
      case 1:
	// m impair: sin(pair)
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * nt ;
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
	    som[i] += 2*xci[i] ;
	    xco[i] = som[i] ;
	    som[i] = -som[i] ;
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
    b = base_r | base_p | T_COSSIN_CP ;
}

			//---------------------	
			// cas T_COSSIN_SI
			//----------------
void _scost_t_cossin_si(Tbl* tb, int & b)
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
  xci += nr * (nt-1) ;
  xco += nr * (nt-1) ;
  
  // Initialisation de som
  for (int i=0 ; i<nr ; i++) {
    xco[i] = 0. ;
    som[i] = 0. ;
  }
  
  // j suivants
  for (int j=nt-2 ; j >= 1 ; j--) {
    // Positionnement
    xci -= nr ;
    xco -= nr ;
    // Calcul
    for (int i=0 ; i<nr ; i++ ) {
      som[i] += 2*xci[i] ;
      xco[i] = som[i] ;
      som[i] = -som[i] ;
    }	// Fin de la boucle sur r
  }   // Fin de la boucle sur theta
  xci -= nr ;
  xco -= nr ;
  // Calcul pour le premier theta
  for (int i=0 ; i<nr ; i++ ) {
    xco[i] =0. ;
  }
  // Positionnement phi suivant
  xci += nr*nt ;
  xco += nr*nt ;
  
  // k=1 
  xci += nr*nt ;
  xco += nr*nt ;
  
  // k >= 2
  for (int k=2 ; k<np+1 ; k++) {
    m = (k/2) % 2 ;	    // Parite de l'harmonique en phi
    
    switch(m) {
    case 0:	    // m pair: sin(impair)
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
      for (int j=nt-2 ; j >= 1 ; j--) {
	// Positionnement
	xci -= nr ;
	xco -= nr ;
	// Calcul
	for (int i=0 ; i<nr ; i++ ) {
	  som[i] += 2*xci[i] ;
	  xco[i] = som[i] ;
	  som[i] = -som[i] ;
	}	// Fin de la boucle sur r
      }   // Fin de la boucle sur theta
      xci -= nr ;
      xco -= nr ;
      // Calcul pour le premier theta
      for (int i=0 ; i<nr ; i++ ) {
	xco[i] =0. ;
      }
      // Positionnement phi suivant
      xci += nr*nt ;
      xco += nr*nt ;
      break ;
      
    case 1: // m impair cos(pair)
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
	  som[i] += 2*xci[i] ;
	  xco[i] = som[i] ;
	  som[i] = - som[i] ;
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
			// cas T_COSSIN_SP
			//----------------
void _scost_t_cossin_sp(Tbl* tb, int & b)
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
  
  // k = 0
  int m = 0 ;	    // Parite de l'harmonique en phi
  
  // Dernier j: j = nt-1
  // Positionnement
  xci += nr * nt ;
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
      som[i] += 2*xci[i] ;
      xco[i] = som[i] ;
      som[i] = -som[i] ;
    }	// Fin de la boucle sur r
  }   // Fin de la boucle sur theta
  // Positionnement phi suivant
  xci -= nr ;
  xci += nr*nt ;
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
      xci += nr * (nt-1) ;
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
	  som[i] += 2*xci[i] ;
	  xco[i] = som[i] ;
	  som[i] = -som[i] ;
	}	// Fin de la boucle sur r
      }   // Fin de la boucle sur theta
      // Normalisation du premier theta
      for (int i=0 ; i<nr ; i++) {
	xco[i] *= .5 ;
      }
      // Positionnement phi suivant
      xci += nr*nt ;
      xco += nr*nt ;
      break ;
      
    case 0:	    // m pair: sin(pair)
      // Dernier j: j = nt-1
      // Positionnement
      xci += nr * nt ;
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
	  som[i] += 2*xci[i] ;
	  xco[i] = som[i] ;
	  som[i] = -som[i] ;
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
  b = base_r | base_p | T_COSSIN_SI ;
    
}

			//----------------------
			// cas T_COSSIN_C
			//----------------------
void _scost_t_cossin_c(Tbl* tb, int & b) {

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	switch(base_r){
	    case(R_CHEBPI_P):
		b = R_CHEBPI_I | base_p | T_COSSIN_C ;
		break ;
	    case(R_CHEBPI_I):
		b = R_CHEBPI_P | base_p | T_COSSIN_C ;
		break ;  
	    default:
		b = base_r | base_p | T_COSSIN_C ;
		break;
	}
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
    double* somP = new double [nr] ;
    double* somN = new double [nr] ;
    
    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) xo[i] = 0 ; 
    
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
	somP[i] = 0. ;
	somN[i] = 0. ;
	xco[i] = 0. ;	// mise a zero du dernier coef en theta
    }
    
    // j suivants
    for (int j=nt-2 ; j >= 0 ; j--) {
	// Positionnement
	xci -= nr ;
	xco -= nr ;
	// Calcul
	for (int i=0 ; i<nr ; i++ ) {
	    if((j%2) == 1) {
		somN[i] = -somN[i] ;
		somN[i] += 2*xci[i] ;
		xco[i] = somP[i] ;
	    }
	    else {
		somP[i] = -somP[i] ;
		somP[i] += 2*xci[i] ; 
		xco[i] = somN[i] ;
	    }
	}	// Fin de la boucle sur r
    }   // Fin de la boucle sur theta
    // j=0 : le facteur 2...
    for (int i=0 ; i<nr ; i++) xco[i] *= .5 ;
 
    // Positionnement phi suivant   
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k=1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	
      	// Dernier j: j = nt-1
	// Positionnement
	xco += nr * (nt-1) ;
	xci += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul	   
	    for (int i=0 ; i<nr ; i++ ) {
		if((j%2) == 1) {
		    somN[i] = -somN[i];
		    somN[i] += 2 * xci[i] ;
		    xco[i] = somP[i] ;
		}
		else {
		    somP[i] = -somP[i];
		    somP[i] += 2 * xci[i] ;
		    xco[i] = somN[i];
		}
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	double fac_m = ( (k/2)%2 == 1 ? 0. : 0.5) ;
	// j=0 : sin(0*theta) ou facteur 2, suivant la parite de m
	for (int i=0 ; i<nr ; i++) xco[i] *= fac_m ;
	
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr*nt ;
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] somP ;
    delete [] somN ;
    
    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    switch(base_r){
	case(R_CHEBPI_P):
	    b = R_CHEBPI_I | base_p | T_COSSIN_C ;
	    break ;
	case(R_CHEBPI_I):
	    b = R_CHEBPI_P | base_p | T_COSSIN_C ;
	    break ;  
	default:
	    b = base_r | base_p | T_COSSIN_C ;
	    break;
    }
}

			//---------------------	
			// cas T_COSSIN_S
			//----------------
void _scost_t_cossin_s(Tbl* tb, int & b) {

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_r = b & MSQ_R ;
	int base_p = b & MSQ_P ;
	switch(base_r){
	    case(R_CHEBPI_P):
		b = R_CHEBPI_I | base_p | T_COSSIN_S ;
		break ;
	    case(R_CHEBPI_I):
		b = R_CHEBPI_P | base_p | T_COSSIN_S ;
		break ;  
	    default:
		b = base_r | base_p | T_COSSIN_S ;
		break;
	}
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
    double* somP = new double [nr] ;
    double* somN = new double [nr] ;

    // pt. sur le tableau de double resultat
    double* xo = new double[(tb->dim).taille] ;
    
    // Initialisation a zero :
    for (int i=0; i<(tb->dim).taille; i++) xo[i] = 0 ; 
    
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
	somP[i] = 0. ;
	somN[i] = 0. ;
	xco[i] = 0. ;	// mise a zero du dernier coef en theta
    }
	
    // j suivants
    for (int j=nt-2 ; j >= 0 ; j--) {
	// Positionnement
	xci -= nr ;
	xco -= nr ;
	// Calcul
	for (int i=0 ; i<nr ; i++ ) {
	    if((j%2) == 1) {
		somN[i] = -somN[i] ;
		somN[i] += 2*xci[i] ;
		xco[i] = somP[i] ;
	    }
	    else {
		somP[i] = -somP[i] ;
		somP[i] += 2*xci[i] ;		  
		xco[i] = somN[i] ;
	    }
	}	// Fin de la boucle sur r
    }   // Fin de la boucle sur theta
	// j=0 : sin(0*theta)
    for (int i=0 ; i<nr ; i++) xco[i] = 0. ;
   
    // Positionnement phi suivant
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k=1
    xci += nr*nt ;
    xco += nr*nt ;
    
    // k >= 2
    for (int k=2 ; k<np+1 ; k++) {
	
      	// Dernier j: j = nt-1
	// Positionnement
	xco += nr * (nt-1) ;
	xci += nr * (nt-1) ;

	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul	   
	    for (int i=0 ; i<nr ; i++ ) {
		if((j%2) == 1) {
		  somN[i] = -somN[i] ;
		  somN[i] += 2*xci[i] ;
		  xco[i] = somP[i] ;
	      }
	      else {
		  somP[i] = -somP[i] ;
		  somP[i] += 2*xci[i] ;		  
		  xco[i] = somN[i] ;
	      }
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	
	double fac_m = ( (k/2)%2 == 0 ? 0. : 0.5) ;
	// j=0 : sin(0*theta) ou facteur 2, suivant la parite de m
	for (int i=0 ; i<nr ; i++) xco[i] *= fac_m ;
	
	// Positionnement phi suivant
	xci += nr*nt ;
	xco += nr*nt ;
  
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    //Menage
    delete [] somP ;
    delete [] somN ;
    
    // base de developpement
    int base_r = b & MSQ_R ;
    int base_p = b & MSQ_P ;
    switch(base_r){
	case(R_CHEBPI_P):
	    b = R_CHEBPI_I | base_p | T_COSSIN_S ;
	    break ;
	case(R_CHEBPI_I):
	    b = R_CHEBPI_P | base_p | T_COSSIN_S ;
	    break ;  
	default:
	    b = base_r | base_p | T_COSSIN_S ;
	    break;
    }
}
