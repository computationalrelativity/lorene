/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char xdsdx_mat_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:16:11  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2007/06/21 20:07:16  k_taniguchi
 * nmax increased to 200
 *
 * Revision 1.2  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/03/16  16:23:17  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

/*
 * Routine renvoyan la matrice de l'operateur x f' = s
 * Pour  l != 1 le resultat est donne en s est donne en Chebyshev et
 * f en Gelerkin (T_i + T_{i+1} pour l pair et (2*i+3)T_i + (2*i+1)T_{i+1} pour 
 * l impair.
 * Pour l=1 pas de probleme de singularite on reste donc en Chebyshev.
 */

//fichiers includes
#include <cstdio>
#include <cstdlib>

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"

		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _xdsdx_mat_pas_prevu(int n, int) {
    cout << "xdsdx_mat pas prevu..." << endl ;
    cout << "n : " << n << endl ;
    abort() ;
    exit(-1) ;
    Matrice res(1, 1) ; // On ne passe jamais ici de toute facon !
    return res;
}


		   //-------------------------
		   //--   CAS R_CHEBP    -----
		   //--------------------------
		    

Matrice _xdsdx_mat_r_chebp (int n, int) {
   const int nmax = 200 ;// Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int nr_dejafait[nmax] ;
    
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if (nr_dejafait[conte] == n)
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_laplacien_nul_mat_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice res(n-1, n-1) ;
    res.set_etat_qcq() ;
    
    double* xdsdx  = new double[n] ;
    
    for (int i=0 ; i<n-1 ; i++) {
	for (int j=0 ; j<n ; j++)
	    xdsdx[j] = 0 ;
	xdsdx[i] = 1 ;
	xdsdx[i+1] = 1 ;
	xdsdx_1d (n, &xdsdx, R_CHEBP) ;
	
	for (int j=0 ; j<n-1 ; j++)
	    res.set(j, i) = xdsdx[j] ; 
    }
    delete [] xdsdx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    }
    
    else
	return *tab[indice] ;
}



		   //------------------------
		   //--   CAS R_CHEBI    ----
		   //------------------------
		    

Matrice _xdsdx_mat_r_chebi (int n, int l) {
   const int nmax = 200 ;// Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int nr_dejafait[nmax] ;
   static int nl_dejafait[nmax] ;
    
   int indice = -1 ;
   // On separe les cas l=1 et l != 1
   int indic_l = (l == 1) ? 1 : 2 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((nr_dejafait[conte] == n) && (nl_dejafait[conte] == indic_l))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_laplacien_nul_mat_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    nr_dejafait[nb_dejafait] = n ;
    nl_dejafait[nb_dejafait] = indic_l ;
    
    // non degenere pour l=1
    int taille = (l==1) ? n : n-1 ;
    Matrice res(taille, taille) ;
    res.set_etat_qcq() ;
  
    double* xdsdx = new double[n] ;
    
    for (int i=0 ; i<taille ; i++) {
	for (int j=0 ; j<n ; j++)
	    xdsdx[j] = 0 ;
	
	// Gelerkin ou Chebyshev ?
	if (taille == n) {
	    xdsdx[i] = 1 ;
	}
	else {
	    xdsdx[i] = 2*i+3 ;
	    xdsdx[i+1] = 2*i+1 ;
	    }
	xdsdx_1d (n, &xdsdx, R_CHEBI) ;  // appel dans le cas impair
	
	for (int j=0 ; j<taille ; j++)
	    res.set(j, i) = xdsdx[j] ;
    }
  
    delete [] xdsdx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    }
    
    else
	return *tab[indice] ;
}


		 //--------------------------
		//- La routine a appeler  ---
	       //----------------------------
Matrice xdsdx_mat(int n, int l, int base_r)
{

		// Routines de derivation
    static Matrice (*xdsdx_mat[MAX_BASE])(int, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    xdsdx_mat[i] = _xdsdx_mat_pas_prevu ;
	}
		// Les routines existantes
	xdsdx_mat[R_CHEBP >> TRA_R] = _xdsdx_mat_r_chebp ;
	xdsdx_mat[R_CHEBI >> TRA_R] = _xdsdx_mat_r_chebi ;
    }
    
    Matrice res(xdsdx_mat[base_r](n, l)) ;
    return res ;
}

