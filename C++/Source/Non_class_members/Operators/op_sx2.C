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


char op_sx2_C[] = "$Header$" ;

/* 
 * Ensemble des routines de base de division par rapport a x
 * (Utilisation interne)
 * 
 *	void _sx2_XXXX(Tbl * t, int & b)
 *	t	pointeur sur le Tbl d'entree-sortie
 *	b	base spectrale
 * On entend par sx2 l'operateur:
 *
 *  --   {f(x) - f(0) - x f'(0)}/x^2		dans le noyau (cas RARE)
 *  --   Id					dans les coquilles (cas FIN)
 *  --   {f(x) - f(1) - (x-1) f'(1)}/(x-1)^2	dans la zone externe compactifiee 
 *						(cas UNSURR) 
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/11/23 15:16:02  m_forot
 *
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/09/07  12:51:34  phil
 * *** empty log message ***
 *
 * Revision 2.3  2000/02/24  16:43:06  eric
 * Initialisation a zero du tableau xo avant le calcul.
 *
 * Revision 2.2  1999/11/15  16:39:30  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 *
 * $Header$
 *
 */


// Fichier includes
#include "tbl.h" 
#include "proto.h" 



		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _sx2_pas_prevu(Tbl * tb, int& base) {
    cout << "sx pas prevu..." << endl ;
    cout << "Tbl: " << tb << " base: " << base << endl ;
    abort () ;
    exit(-1) ;
}

			//-------------
			// Identite ---
			//-------------

void _sx2_identite(Tbl* , int& ) {
    return ;
}

			//---------------
			// cas R_CHEBP --
			//--------------
void _sx2_r_chebp(Tbl* tb, int&)
    {
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	    }
	    
	else {
	for (int j=0 ; j<nt ; j++) {

	    double somp, somn ;
	    int sgn = 1 ;
	    
	    xco[nr-1] = 0 ;
	    somp = 4 * sgn * (nr-1) * xci[nr-1] ;
	    somn = 2 * sgn * xci[nr-1] ;
	    xco[nr-2] = somp - 2*(nr-2)*somn ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		somp += 4 * (i+1) * sgn * xci[i+1] ;
		somn += 2 * sgn * xci[i+1] ;
		xco[i] = somp - 2*i * somn ;
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
    // inchangee
}

			//----------------
			// cas R_CHEBI ---
			//----------------

void _sx2_r_chebi(Tbl* tb, int&)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	    }
	else {
	for (int j=0 ; j<nt ; j++) {

	    double somp, somn ;
	    int sgn = 1 ;
	    
	    xco[nr-1] = 0 ;
	    somp = 2 * sgn * (2*(nr-1)+1) * xci[nr-1] ;
	    somn = 2 * sgn * xci[nr-1] ;
	    xco[nr-2] = somp - (2*(nr-2)+1)*somn ;
	    for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		somp += 2 * (2*(i+1)+1) * sgn * xci[i+1] ;
		somn += 2 * sgn * xci[i+1] ;
		xco[i] = somp - (2*i+1) * somn ;
	    }	// Fin de la premiere boucle sur r
	    for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	    }	// Fin de la deuxieme boucle sur r
	    
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

			//--------------------
			// cas R_CHEBPIM_P --
			//------------------

void _sx2_r_chebpim_p(Tbl* tb, int&)
{
  
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
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
    int auxiliaire ;
    
    // Partie paire
    xci = xi ;
    xco = xo ;
    for (int k=0 ; k<np+1 ; k += 4) {
	auxiliaire = (k==np) ? 1 : 2 ; // evite de prendre le dernier coef
	for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
	
	    // On evite le coefficient de sin (0*phi)
	if ((k==0) && (kmod==1)) { 
		xco += nr*nt ;
		xci += nr*nt ;
		}
		
	else   
	    for (int j=0 ; j<nt ; j++) {

		double somp, somn ;
		int sgn = 1 ;
	    
		xco[nr-1] = 0 ;
		somp = 4 * sgn * (nr-1) * xci[nr-1] ;
		somn = 2 * sgn * xci[nr-1] ;
		xco[nr-2] = somp - 2*(nr-2)*somn ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    sgn = - sgn ;
		    somp += 4 * (i+1) * sgn * xci[i+1] ;
		    somn += 2 * sgn * xci[i+1] ;
		    xco[i] = somp - 2*i * somn ;
		}	// Fin de la premiere boucle sur r
		for (int i=0 ; i<nr ; i+=2) {
		    xco[i] = -xco[i] ;
		}	// Fin de la deuxieme boucle sur r
		xco[0] *= .5 ;
		
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

		double somp, somn ;
		int sgn = 1 ;
	    
		xco[nr-1] = 0 ;
		somp = 2 * sgn * (2*(nr-1)+1) * xci[nr-1] ;
		somn = 2 * sgn * xci[nr-1] ;
		xco[nr-2] = somp - (2*(nr-2)+1)*somn ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    sgn = - sgn ;
		    somp += 2 * (2*(i+1)+1) * sgn * xci[i+1] ;
		    somn += 2 * sgn * xci[i+1] ;
		    xco[i] = somp - (2*i+1) * somn ;
		}	// Fin de la premiere boucle sur r
		for (int i=0 ; i<nr ; i+=2) {
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
    // inchangee
}


			//-------------------
			// cas R_CHEBPIM_I --
			//-------------------

void _sx2_r_chebpim_i(Tbl* tb, int&)
{

   // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
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

		double somp, somn ;
		int sgn = 1 ;
	    
		xco[nr-1] = 0 ;
		somp = 2 * sgn * (2*(nr-1)+1) * xci[nr-1] ;
		somn = 2 * sgn * xci[nr-1] ;
		xco[nr-2] = somp - (2*(nr-2)+1)*somn ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    sgn = - sgn ;
		    somp += 2 * (2*(i+1)+1) * sgn * xci[i+1] ;
		    somn += 2 * sgn * xci[i+1] ;
		    xco[i] = somp - (2*i+1) * somn ;
		}	// Fin de la premiere boucle sur r
		for (int i=0 ; i<nr ; i+=2) {
		    xco[i] = -xco[i] ;
		}	// Fin de la deuxieme boucle sur r
	    
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

		double somp, somn ;
		int sgn = 1 ;
	    
		xco[nr-1] = 0 ;
		somp = 4 * sgn * (nr-1) * xci[nr-1] ;
		somn = 2 * sgn * xci[nr-1] ;
		xco[nr-2] = somp - 2*(nr-2)*somn ;
		for (int i = nr-3 ; i >= 0 ; i-- ) {
		    sgn = - sgn ;
		    somp += 4 * (i+1) * sgn * xci[i+1] ;
		    somn += 2 * sgn * xci[i+1] ;
		    xco[i] = somp - 2*i * somn ;
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
    // inchangee
}

			//---------------
			// cas R_CHEBU --
			//---------------

void _sx2_r_chebu(Tbl* tb, int&)
{
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
	return ;
    }
    
    // Pour le confort
    int nr = (tb->dim).dim[0] ;	    // Nombre
    int nt = (tb->dim).dim[1] ;	    //	 de points
    int np = (tb->dim).dim[2] ;	    //	    physiques REELS
    np = np - 2 ;		   
      
    int ntnr = nt*nr ; 
    
    for (int k=0 ; k<np+1 ; k++) {
	if (k==1) continue ;	// On ne traite pas le coefficient de sin(0*phi) 
	for (int j=0 ; j<nt ; j++) {
	
	    double* cf = tb->t + k*ntnr + j*nr ;
	    sxm1_1d_cheb(nr, cf) ;	// 1ere division par (x-1)
	    sxm1_1d_cheb(nr, cf) ;	// 2eme division par (x-1) 
	    
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // base de developpement
    // inchangee
}

			//---------------
			// cas R_CHEBPI_P --
			//--------------
void _sx2_r_chebpi_p(Tbl* tb, int&)
    {
    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	    }
	    
	else {
	for (int j=0 ; j<nt ; j++) {
	    int l = j%2 ;
	    if(l==0){
	      double somp, somn ;
	      int sgn = 1 ;
	      
	      xco[nr-1] = 0 ;
	      somp = 4 * sgn * (nr-1) * xci[nr-1] ;
	      somn = 2 * sgn * xci[nr-1] ;
	      xco[nr-2] = somp - 2*(nr-2)*somn ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		somp += 4 * (i+1) * sgn * xci[i+1] ;
		somn += 2 * sgn * xci[i+1] ;
		xco[i] = somp - 2*i * somn ;
	      }	// Fin de la premiere boucle sur r
	      for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	      }	// Fin de la deuxieme boucle sur r
	      xco[0] *= .5 ;
	    } else {
	      double somp, somn ;
	      int sgn = 1 ;
	      
	      xco[nr-1] = 0 ;
	      somp = 2 * sgn * (2*(nr-1)+1) * xci[nr-1] ;
	      somn = 2 * sgn * xci[nr-1] ;
	      xco[nr-2] = somp - (2*(nr-2)+1)*somn ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		somp += 2 * (2*(i+1)+1) * sgn * xci[i+1] ;
		somn += 2 * sgn * xci[i+1] ;
		xco[i] = somp - (2*i+1) * somn ;
	      }	// Fin de la premiere boucle sur r
	      for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	      }	// Fin de la deuxieme boucle sur r
	    }
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

			//----------------
			// cas R_CHEBPI_I ---
			//----------------

void _sx2_r_chebpi_i(Tbl* tb, int&)
{

    // Peut-etre rien a faire ?
    if (tb->get_etat() == ETATZERO) {
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
    
    for (int k=0 ; k<np+1 ; k++) 
	if (k==1) {
	    xci += nr*nt ;
	    xco += nr*nt ;
	    }
	else {
	for (int j=0 ; j<nt ; j++) {
	    int l = j%2 ;
	    if(l==1){
	      double somp, somn ;
	      int sgn = 1 ;
	      
	      xco[nr-1] = 0 ;
	      somp = 4 * sgn * (nr-1) * xci[nr-1] ;
	      somn = 2 * sgn * xci[nr-1] ;
	      xco[nr-2] = somp - 2*(nr-2)*somn ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		somp += 4 * (i+1) * sgn * xci[i+1] ;
		somn += 2 * sgn * xci[i+1] ;
		xco[i] = somp - 2*i * somn ;
	      }	// Fin de la premiere boucle sur r
	      for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	      }	// Fin de la deuxieme boucle sur r
	      xco[0] *= .5 ;
	    } else {
	      double somp, somn ;
	      int sgn = 1 ;
	      
	      xco[nr-1] = 0 ;
	      somp = 2 * sgn * (2*(nr-1)+1) * xci[nr-1] ;
	      somn = 2 * sgn * xci[nr-1] ;
	      xco[nr-2] = somp - (2*(nr-2)+1)*somn ;
	      for (int i = nr-3 ; i >= 0 ; i-- ) {
		sgn = - sgn ;
		somp += 2 * (2*(i+1)+1) * sgn * xci[i+1] ;
		somn += 2 * sgn * xci[i+1] ;
		xco[i] = somp - (2*i+1) * somn ;
	      }	// Fin de la premiere boucle sur r
	      for (int i=0 ; i<nr ; i+=2) {
		xco[i] = -xco[i] ;
	      }	// Fin de la deuxieme boucle sur r
	    }
	    xci += nr ;
	    xco += nr ;
	}   // Fin de la boucle sur theta
    }	// Fin de la boucle sur phi
    
    // On remet les choses la ou il faut
    delete [] tb->t;
    tb->t = xo ;
    
    // base de developpement
    // inchangee
}

