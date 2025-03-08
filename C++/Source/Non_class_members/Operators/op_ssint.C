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
 * Revision 1.10  2024/01/26 17:44:25  g_servignat
 * Updated the Pseudopolytrope_1D class to be consistent with the paper (i.e. with a GPP in the middle)
 *
 * Revision 1.9  2016/12/05 16:18:08  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.8  2014/10/13 08:53:26  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2009/10/10 18:28:11  j_novak
 * New bases T_COS and T_SIN.
 *
 * Revision 1.6  2007/12/21 10:43:38  j_novak
 * Corrected some bugs in the case nt=1 (1 point in theta).
 *
 * Revision 1.5  2007/10/05 12:37:20  j_novak
 * Corrected a few errors in the theta-nonsymmetric case (bases T_COSSIN_C and
 * T_COSSIN_S).
 *
 * Revision 1.4  2005/02/16 15:29:48  m_forot
 * Correct T_COSSIN_S and T_COSSIN_C cases
 *
 * Revision 1.3  2004/12/17 13:20:18  m_forot
 * Modify the case T_COSSIN_C and T_COSSIN_S
 *
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

namespace Lorene {
void _ssint_pas_prevu(Tbl * tb, int& base) {
    cout << "ssint pas prevu..." << endl ;
    cout << "Tbl: " << tb << " base: " << base << endl ;
    abort () ;
    exit(-1) ;
}

		    //--------------
		    // cas T_COS
		    //--------------
		    
void _ssint_t_cos(Tbl* tb, int & b)
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
    
    double cx ;

    // k = 0
        
	// Dernier j: j = nt-1
	// Positionnement
	xci += nr * (nt-1) ;
	cx = -2 ;
	xco += nr * (nt-1) ;
	
	// Initialisation de som 
	for (int i=0 ; i<nr ; i++) {
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;	// mise a zero du dernier coef en theta
	}
	
	// j suivants
	for (int j=nt-2 ; j >0 ; j--) {
	  int l = j % 2 ;
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
	      if(l==1) somN[i] += cx * xci[i] ;
	      else somP[i] += cx * xci[i] ;	      
	      xco[i] = somN[i]*(1-l)+somP[i]*l ; 
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// j = 0 sin(0*theta)
	xci -= nr ;
	xco -= nr ;
	for (int i=0 ; i<nr ; i++) {
	  xco[i] = 0 ;
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
	cx = -2 ;
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
	  int l = j% 2 ;
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
	      if(l==1) somN[i] += -2 * xci[i] ;
	      else somP[i] += -2 * xci[i] ;
	      xco[i] = somN[i]*(1-l)+somP[i]*l ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	for (int i=0 ; i<nr ; i++) {
	  xco[i] = 0. ;
	}
	
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
			
			//------------
			// cas T_SIN
			//------------
			
void _ssint_t_sin(Tbl* tb, int & b)
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
      int l = j % 2 ;
      // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
	      if(l==1) somN[i] += 2 * xci[i] ;
	      else somP[i] += 2 * xci[i] ;	      
	      xco[i] = somN[i]*(1-l)+somP[i]*l ; 
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
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
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	  int l = j % 2 ;
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul
	    for (int i=0 ; i<nr ; i++ ) {
	      if(l==1) somN[i] += 2 * xci[i] ;
	      else somP[i] += 2 * xci[i] ;
	      xco[i] = somN[i]*(1-l)+somP[i]*l ;
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

    // // Regularisation
    // for (int k=0; k<np+1; k++)
    // {
    //   if (k==1) continue;
    //   for (int i=0; i<nr; i++)
    //   {
    //     double val_axis = 0.;
    //     for (int j=0; j<nt-1; j++)
    //       val_axis += (*tb)(k, j, i);
    //     tb->set(k, nt-1, i) = -val_axis;
    //   }
    // }
     
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
    if (nt > 1) {
	xco -= nr ;
	for (int i=0 ; i<nr ; i++) xco[i] = 0 ;
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
	    cx = -2 ;
	    xco += nr * (nt-1) ;

	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >0 ; j--) {
	  int l = j % 2 ;
	  // Positionnement
	  xci -= nr ;
	  xco -= nr ;
	  // Calcul	   
	  for (int i=0 ; i<nr ; i++ ) {
	    if(l==1) somN[i] += cx * xci[i] ;
	    else somP[i] += cx * xci[i] ;	      
	    xco[i] = somN[i]*(1-l)+somP[i]*l ; 
	  }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
	// j = 0 sin(0*theta)
	xci -= nr ;
	xco -= nr ;
	for (int i=0 ; i<nr ; i++) {
	  xco[i] = 0 ;
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
	
	// Dernier j: j = nt-1
	// Positionnement
	double fac_j0 = 0 ;
	if (m == 0) {
	  cx = -2 ;
	    }
	else {
	  cx = 2 ;
	  fac_j0 = 0.5 ;
	}
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
	  int l = j % 2 ;
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul	   
	    for (int i=0 ; i<nr ; i++ ) {
	      if(l==1) somN[i] += cx * xci[i] ;
	      else somP[i] += cx * xci[i] ;
	      xco[i] = somN[i]*(1-l)+somP[i]*l ;
	    }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta

	// Normalisation du premier theta dans le cas sin, mise a zero dans le cas cos
	for (int i=0 ; i<nr ; i++) {
	  xco[i] *= fac_j0 ;
	}
	
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

			//---------------------	
			// cas T_COSSIN_S
			//----------------

void _ssint_t_cossin_s(Tbl* tb, int & b)
{


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
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	  int l = j % 2 ;
	  // Positionnement
	  xci -= nr ;
	  xco -= nr ;
	  // Calcul	   
	  for (int i=0 ; i<nr ; i++ ) {
	    if(l==1) somN[i] += cx * xci[i] ;
	    else somP[i] += cx * xci[i] ;	      
	    xco[i] = somN[i]*(1-l)+somP[i]*l ; 
	  }	// Fin de la boucle sur r
	}   // Fin de la boucle sur theta
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
	
	// Dernier j: j = nt-1
	// Positionnement
	if (m == 0) {
	    xci += nr * (nt-1) ;
	    cx = 2 ;
	    }
	else {
	    xci += nr * (nt-1) ;
	    cx = -2 ;
	}
	xco += nr * (nt-1) ;
	
	// Initialisation de som
	for (int i=0 ; i<nr ; i++) {
	    somP[i] = 0. ;
	    somN[i] = 0. ;
	    xco[i] = 0. ;
	}
	
	// j suivants
	for (int j=nt-2 ; j >= 0 ; j--) {
	  int l = j % 2 ;
	    // Positionnement
	    xci -= nr ;
	    xco -= nr ;
	    // Calcul	   
	    for (int i=0 ; i<nr ; i++ ) {
	      if(l==1) somN[i] += cx * xci[i] ;
	      else somP[i] += cx * xci[i] ;
	      xco[i] = somN[i]*(1-l)+somP[i]*l ;
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
}
