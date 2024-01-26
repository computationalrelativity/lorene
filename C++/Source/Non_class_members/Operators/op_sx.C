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
 * Ensemble des routines de base de division par rapport a x
 * (Utilisation interne)
 * 
 *	void _sx_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * 
 */
 
 /*
 * $Id$
 * $Log$
 * Revision 1.7  2024/01/26 17:44:25  g_servignat
 * Updated the Pseudopolytrope_1D class to be consistent with the paper (i.e. with a GPP in the middle)
 *
 * Revision 1.5  2016/12/05 16:18:08  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2015/03/05 08:49:32  j_novak
 * Implemented operators with Legendre bases.
 *
 * Revision 1.3  2014/10/13 08:53:26  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2004/11/23 15:16:01  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.5  2000/09/07  12:50:48  phil
 * *** empty log message ***
 *
 * Revision 2.4  2000/02/24  16:42:59  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 2.3  1999/11/15  16:39:17  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.2  1999/10/18  13:40:11  eric
 * Suppression de la routine _sx_r_chebu car doublon avec _sxm1_cheb.
 *
 * Revision 2.1  1999/09/13  11:35:42  phil
 * gestion des cas k==1 et k==np
 *
 * Revision 2.0  1999/04/26  14:59:42  phil
 * *** empty log message ***
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
void _sx_pas_prevu(Tbl * tb, int& base) {
    cout << "sx pas prevu..." << endl ;
    cout << "Tbl: " << tb << " base: " << base << endl ;
    abort () ;
    exit(-1) ;
}

			//-------------
			// Identite ---
			//-------------

void _sx_identite(Tbl* , int& ) {
    return ;
}

			//---------------
			// cas R_CHEBP --
			//--------------

void _sx_r_chebp(Tbl* tb, int& base)
    {
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = base & MSQ_T ;
	int base_p = base & MSQ_P ;
	base = base_p | base_t | R_CHEBI ;
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // pt. sur le tableau de double resultat
    double* xo = new double [tb->get_taille()];
    
    // Initialisation a zero :
    for (int i=0; i<tb->get_taille(); i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants

    int borne_phi = np + 1 ; 
    if (np == 1) {
	borne_phi = 1 ; 
    }
    
    for (int k=0 ; k< borne_phi ; k++)
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som ;
	    int sgn = 1 ;
		
	    // Regularisation (cf Jordan Nicoules thesis)
		double val_ori = 0 ;
		for (int i=0 ; i<nr-1 ; i++) {
			val_ori += sgn*xci[i] ;
			sgn = -sgn ;
		}
		sgn = -sgn ;
		xci[nr-1] = sgn*val_ori ;
		sgn = 1 ;

		xco[nr-1] = 0 ;
	    som = 2 * sgn * xci[nr-1] ;
	    xco[nr-2] = som ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		som += 2 * sgn * xci[i+1] ;
		xco[i] = som ;
	    }	// Fin de la premiere boucle sur r
	    for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	    }	// Fin de la deuxieme boucle sur r
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // pair -> impair
    int base_t = base & MSQ_T ;
    int base_p = base & MSQ_P ;
    base = base_p | base_t | R_CHEBI ;

}

			//----------------
			// cas R_CHEBI ---
			//----------------

void _sx_r_chebi(Tbl* tb, int& base)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = base & MSQ_T ;
	int base_p = base & MSQ_P ;
	base = base_p | base_t | R_CHEBP ;
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // pt. sur le tableau de double resultat
    double* xo = new double [tb->get_taille()];
    
    // Initialisation a zero :
    for (int i=0; i<tb->get_taille(); i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    int borne_phi = np + 1 ; 
    if (np == 1) {
	borne_phi = 1 ; 
    }
    
    for (int k=0 ; k< borne_phi ; k++) 
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {
	    double som ;
	    int sgn = 1 ;
	    
	    xco[nr-1] = 0 ;
	    som = 2 * sgn * xci[nr-2] ;
	    xco[nr-2] = som ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		som += 2 * sgn * xci[i] ;
		xco[i] = som ;
	    }	// Fin de la premiere boucle sur r
	    for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
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
    // impair -> pair
    int base_t = base & MSQ_T ;
    int base_p = base & MSQ_P ;    
    base = base_p | base_t | R_CHEBP ;
}

			//--------------------
			// cas R_CHEBPIM_P --
			//------------------

void _sx_r_chebpim_p(Tbl* tb, int& base)
{
  
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = base & MSQ_T ;
	int base_p = base & MSQ_P ;
	base = base_p | base_t | R_CHEBPIM_I ;
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		   
   
    // pt. sur le tableau de double resultat
    double* xo = new double [tb->get_taille()];

    // Initialisation a zero :
    for (int i=0; i<tb->get_taille(); i++) {
	xo[i] = 0 ; 
    }
        
    // On y va...
    double* xi = tb->t ;
    double* xci ;	// Pointeurs
    double* xco ;	//  courants
    
    
    // Partie paire
    xci = xi ;
    xco = xo ;
  
    int auxiliaire ;
    
    for (int k=0 ; k<np+1 ; k += 4)  {
	
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
		int sgn = 1 ;
	    
		xco[nr-1] = 0 ;
		som = 2 * sgn * xci[nr-1] ;
		xco[nr-2] = som ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    sgn = - sgn ;
		    som += 2 * sgn * xci[i+1] ;
		    xco[i] = som ;
		}	// Fin de la premiere boucle sur r
		for (int i=0 ; i<nr ; i+=2) {
		    xco[i] = -xco[i] ;
		}	// Fin de la deuxieme boucle sur r
		
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
		int sgn = 1 ;
	    
		xco[nr-1] = 0 ;
		som = 2 * sgn * xci[nr-2] ;
		xco[nr-2] = som ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    sgn = - sgn ;
		    som += 2 * sgn * xci[i] ;
		    xco[i] = som ;
		}	// Fin de la premiere boucle sur r
		for (int i=0 ; i<nr ; i+=2) {
		    xco[i] = -xco[i] ;
		}	// Fin de la deuxieme boucle sur r
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
    int base_t = base & MSQ_T ;
    int base_p = base & MSQ_P ;
    base = base_p | base_t | R_CHEBPIM_I ;
}

			//-------------------
			// cas R_CHEBPIM_I --
			//-------------------

void _sx_r_chebpim_i(Tbl* tb, int& base)
{

   // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = base & MSQ_T ;
	int base_p = base & MSQ_P ;
	base = base_p | base_t | R_CHEBPIM_P ;
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		   
   
    // pt. sur le tableau de double resultat
    double* xo = new double [tb->get_taille()];
    
    // Initialisation a zero :
    for (int i=0; i<tb->get_taille(); i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci ;	// Pointeurs
    double* xco ;	//  courants
    
    // Partie impaire
    xci = xi ;
    xco = xo ;
    
      
    int auxiliaire ;
    
    for (int k=0 ; k<np+1 ; k += 4)  {
	
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
		int sgn = 1 ;
	    
		xco[nr-1] = 0 ;
		som = 2 * sgn * xci[nr-2] ;
		xco[nr-2] = som ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    sgn = - sgn ;
		    som += 2 * sgn * xci[i] ;
		    xco[i] = som ;
		}	// Fin de la premiere boucle sur r
		for (int i=0 ; i<nr ; i+=2) {
		    xco[i] = -xco[i] ;
		}	// Fin de la deuxieme boucle sur r
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
		int i ;
		int sgn = 1 ;
	    
		xco[nr-1] = 0 ;
		som = 2 * sgn * xci[nr-1] ;
		xco[nr-2] = som ;
		for (i = nr-3 ; i >= 0 ; i-- ) {
		    sgn = - sgn ;
		    som += 2 * sgn * xci[i+1] ;
		    xco[i] = som ;
		}	// Fin de la premiere boucle sur r
		for (i=0 ; i<nr ; i+=2) {
		    xco[i] = -xco[i] ;
		}	// Fin de la deuxieme boucle sur r
	    
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
    int base_t = base & MSQ_T ;
    int base_p = base & MSQ_P ;
    base = base_p | base_t | R_CHEBPIM_P ;
}

			//------------------
			// cas R_CHEBPI_P --
			//------------------

void _sx_r_chebpi_p(Tbl* tb, int& base)
    {
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = base & MSQ_T ;
	int base_p = base & MSQ_P ;
	base = base_p | base_t | R_CHEBPI_I ;
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // pt. sur le tableau de double resultat
    double* xo = new double [tb->get_taille()];
    
    // Initialisation a zero :
    for (int i=0; i<tb->get_taille(); i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants

    int borne_phi = np + 1 ; 
    if (np == 1) {
	borne_phi = 1 ; 
    }
    
    for (int k=0 ; k< borne_phi ; k++)
      if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	}
	else {
	for (int j=0 ; j<nt ; j++) {
	    int l = j%2;
	    if(l==0){
	      double som ;
	      int sgn = 1 ;
	      
	      xco[nr-1] = 0 ;
	      som = 2 * sgn * xci[nr-1] ;
	      xco[nr-2] = som ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		som += 2 * sgn * xci[i+1] ;
		xco[i] = som ;
	      }	// Fin de la premiere boucle sur r
	      for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	      }	// Fin de la deuxieme boucle sur r
	    } else {
	      double som ;
	      int sgn = 1 ;
	      
	      xco[nr-1] = 0 ;
	      som = 2 * sgn * xci[nr-2] ;
	      xco[nr-2] = som ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		som += 2 * sgn * xci[i] ;
		xco[i] = som ;
	      }	// Fin de la premiere boucle sur r
	      for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	      }	// Fin de la deuxieme boucle sur r
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
    int base_t = base & MSQ_T ;
    int base_p = base & MSQ_P ;
    base = base_p | base_t | R_CHEBPI_I ;

}

			//-------------------
			// cas R_CHEBPI_I ---
			//-------------------

void _sx_r_chebpi_i(Tbl* tb, int& base)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = base & MSQ_T ;
	int base_p = base & MSQ_P ;
	base = base_p | base_t | R_CHEBPI_P ;
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // pt. sur le tableau de double resultat
    double* xo = new double [tb->get_taille()];
    
    // Initialisation a zero :
    for (int i=0; i<tb->get_taille(); i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    int borne_phi = np + 1 ; 
    if (np == 1) {
	borne_phi = 1 ; 
    }
    
    for (int k=0 ; k< borne_phi ; k++) 
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {
	    int l = j%2;
	    if(l==1){
	      double som ;
	      int sgn = 1 ;
	      
	      xco[nr-1] = 0 ;
	      som = 2 * sgn * xci[nr-1] ;
	      xco[nr-2] = som ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		som += 2 * sgn * xci[i+1] ;
		xco[i] = som ;
	      }	// Fin de la premiere boucle sur r
	      for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	      }	// Fin de la deuxieme boucle sur r
	    } else {
	      double som ;
	      int sgn = 1 ;
	      
	      xco[nr-1] = 0 ;
	      som = 2 * sgn * xci[nr-2] ;
	      xco[nr-2] = som ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		som += 2 * sgn * xci[i] ;
		xco[i] = som ;
	      }	// Fin de la premiere boucle sur r
	      for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	      }	// Fin de la deuxieme boucle sur r
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
    int base_t = base & MSQ_T ;
    int base_p = base & MSQ_P ;    
    base = base_p | base_t | R_CHEBPI_P ;
}

			//--------------
			// cas R_LEGP --
			//--------------

void _sx_r_legp(Tbl* tb, int& base)

    {
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = base & MSQ_T ;
	int base_p = base & MSQ_P ;
	base = base_p | base_t | R_LEGI ;
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // pt. sur le tableau de double resultat
    double* xo = new double [tb->get_taille()];
    
    // Initialisation a zero :
    for (int i=0; i<tb->get_taille(); i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants

    int borne_phi = np + 1 ; 
    if (np == 1) {
	borne_phi = 1 ; 
    }
    
    for (int k=0 ; k< borne_phi ; k++)
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	}
	else {
	for (int j=0 ; j<nt ; j++) {

	    double som = 0 ;
	    
	    xco[nr-1] = 0 ;
	    for (int i=nr - 2; i>=0; i--) {
	      som += xci[i+1] ;
	      xco[i] = double(4*i+3)/double(2*i+2)*som ;
	      som *= -double(2*i+1)/double(2*i+2) ;
	    }	//Fin de la boucle sur r
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // pair -> impair
    int base_t = base & MSQ_T ;
    int base_p = base & MSQ_P ;
    base = base_p | base_t | R_LEGI ;

}

			//---------------
			// cas R_LEGI ---
			//---------------

void _sx_r_legi(Tbl* tb, int& base)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	int base_t = base & MSQ_T ;
	int base_p = base & MSQ_P ;
	base = base_p | base_t | R_LEGP ;
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		    // Nombre de points physiques
    
    // pt. sur le tableau de double resultat
    double* xo = new double [tb->get_taille()];
    
    // Initialisation a zero :
    for (int i=0; i<tb->get_taille(); i++) {
	xo[i] = 0 ; 
    }
    
    // On y va...
    double* xi = tb->t ;
    double* xci = xi ;	// Pointeurs
    double* xco = xo ;	//  courants
    
    int borne_phi = np + 1 ; 
    if (np == 1) {
	borne_phi = 1 ; 
    }
    
    for (int k=0 ; k< borne_phi ; k++) 
	if (k == 1)  {
		xci += nr*nt ;
		xco += nr*nt ;
		}
	else {
	for (int j=0 ; j<nt ; j++) {
	    double som = 0 ;
	    
	    xco[nr-1] = 0 ;
	    for (int i = nr-2 ; i >= 0 ; i-- ) {
	      som += xci[i] ;
	      xco[i] = double(4*i+1)/double(2*i+1)*som ;
	      som *= -double(2*i)/double(2*i+1) ;
	    }	// Fin de la boucle sur r
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t ;
    tb->t = xo ;
    
    // base de developpement
    // impair -> pair
    int base_t = base & MSQ_T ;
    int base_p = base & MSQ_P ;    
    base = base_p | base_t | R_LEGP ;
}


}
