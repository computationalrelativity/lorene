/*
 *  Methods preparing the operators of Ope_pois_vect_r for inversion.
 *
 *    (see file ope_elementary.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char prepa_pvect_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/05/17 15:42:21  j_novak
 * The method 1 of vector Poisson eq. solves directly for F^r.
 * Some bugs were corrected in the operator poisson_vect.
 *
 * Revision 1.1  2004/05/10 15:28:22  j_novak
 * First version of functions for the solution of the r-component of the
 * vector Poisson equation.
 *
 *
 * $Header$
 *
 */

//fichiers includes
#include <stdlib.h>

#include "matrice.h"
#include "type_parite.h"

/*
 * Fonctions supprimant le nombre de colonnes (les premieres)
  et de lignes (les dernieres) a l'operateur renvoye par laplacien_mat, de facon
  a ce qu'il ne soit plus degenere. Ceci doit etre fait apres les combinaisons
  lineaires. La mise a bandes et la decomposition LU sont egalement effectuees ici
  
  
  Entree : lap :  resultat de laplacien_mat
	    l : associe a lap
	    puis : puissance dans la ZEC
	    base_r : base de developpement
	    
  Sortie : renvoie un operateur non degenere ....
 */
 
Matrice _nondeg_pvect_r_pas_prevu(const Matrice &, int , double, int) ;
Matrice _nondeg_pvect_r_cheb (const Matrice&, int, double, int) ;
Matrice _nondeg_pvect_r_chebp (const Matrice&, int, double, int) ;
Matrice _nondeg_pvect_r_chebi (const Matrice&, int, double, int) ;
Matrice _nondeg_pvect_r_chebu (const Matrice&, int, double, int) ; 


		//------------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _nondeg_pvect_r_pas_prevu(const Matrice &lap, int l, double echelle, int puis) {
    cout << "Construction non degeneree pas prevue..." << endl ;
    cout << "l : " << l << endl ;
    cout << "lap : " << lap << endl ;
    cout << "echelle : " << echelle << endl ;
    cout << " puis : " << puis << endl ;
    abort() ;
    exit(-1) ;
    Matrice res(1, 1) ;
    return res;
}



	     	//-------------------
	       //--  R_CHEB   -------
	      //--------------------

Matrice _nondeg_pvect_r_cheb (const Matrice &lap, int l, double echelle, int) {
    
    
    int n = lap.get_dim(0) ;
    
   const int nmax = 200 ; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
   static double vieux_echelle = 0;
   
   // Si on a change l'echelle : on detruit tout :
   if (vieux_echelle != echelle) {
       for (int i=0 ; i<nb_dejafait ; i++) {
	   l_dejafait[i] = -1 ;
	   nr_dejafait[i] = -1 ;
	   delete tab[i] ;	 
       }
	vieux_echelle = echelle ;
	 nb_dejafait = 0 ;
   }
      
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_nondeg_pvect_r_cheb : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
    //assert (l<n) ;
    
    Matrice res(n-2, n-2) ;
    res.set_etat_qcq() ;
    for (int i=0 ; i<n-2 ; i++)
	for (int j=0 ; j<n-2 ; j++)
	    res.set(i, j) = lap(i, j+2) ;
	
    res.set_band(2, 2) ;
    res.set_lu() ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;  
}




	     	//------------------
	       //--  R_CHEBP   ----
	      //------------------
	      
Matrice _nondeg_pvect_r_chebp (const Matrice &lap, int l, double, int) {
    
    int n = lap.get_dim(0) ;
    
       
   const int nmax = 200 ; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
    
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_nondeg_pvect_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    assert (div(l, 2).rem == 1) ;
    
    if (l==1) {
	Matrice res(n-1, n-1) ;
	res.set_etat_qcq() ;
	for (int i=0 ; i<n-1 ; i++)
	    for (int j=0 ; j<n-1 ; j++)
		res.set(i, j) = lap(i, j+1) ;
	res.set_band(3, 0) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
    }
    else {
	Matrice res(n-2, n-2) ;
	res.set_etat_qcq() ;
	for (int i=0 ;i<n-2 ; i++)
	    for (int j=0 ; j<n-2 ; j++)
		res.set(i, j) = lap(i, j+2) ;
	
	res.set_band(2, 1) ;	
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
     }
    }
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}




	     	//-------------------
	       //--  R_CHEBI   -----
	      //-------------------
	      
Matrice _nondeg_pvect_r_chebi (const Matrice &lap, int l, double, int) {
    
    int n = lap.get_dim(0) ;
       
   const int nmax = 200 ; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
    
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_nondeg_pvect_r_chebi : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
    assert (div(l, 2).rem == 0) ;
  //  assert (l<=2*n-1) ;
    
    if (l<=2) { //### to be checked!!!
	Matrice res(n-1, n-1) ;
	res.set_etat_qcq() ;
	for (int i=0 ; i<n-1 ; i++)
	    for (int j=0 ; j<n-1 ; j++)
		res.set(i, j) = lap(i, j+1) ;
	res.set_band(3, 0) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
    }
    else {
	Matrice res(n-2, n-2) ;
	res.set_etat_qcq() ;
	for (int i=0 ;i<n-2 ; i++)
	    for (int j=0 ; j<n-2 ; j++)
		res.set(i, j) = lap(i, j+2) ;
	
	res.set_band(2, 1) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
	} 
    }
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}




	     	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
	      
	      
Matrice _nondeg_pvect_r_chebu (const Matrice &lap, int l, double, int puis) {

  if (puis != 4) {
    cout << "_ope_pvect_r_mat_r_chebu : only the case dzpuis = 4 "
	 << '\n' << "is implemented! \n"
	 << "dzpuis = " << puis << endl ;
    abort() ;
  }
  int n = lap.get_dim(0) ;
    
  const int nmax = 200; // Nombre de Matrices stockees
  static Matrice* tab[nmax] ;  // les matrices calculees
  static int nb_dejafait = 0 ; // nbre de matrices calculees
  static int l_dejafait[nmax] ;
  static int nr_dejafait[nmax] ;
    
  int indice = -1 ;
   
  // On determine si la matrice a deja ete calculee :
  for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
      indice = conte ;
    
  // Calcul a faire : 
  if (indice  == -1) {
    if (nb_dejafait >= nmax) {
      cout << "_nondeg_pvect_r_chebu : trop de matrices" << endl ;
      abort() ;
      exit (-1) ;
    }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice res(n-3, n-3) ;
    res.set_etat_qcq() ;
    for (int i=0 ;i<n-3 ; i++)
      for (int j=0 ; j<n-3 ; j++)
	res.set(i, j) = lap(i, j+3) ;
    
    res.set_band(2, 1) ;
    res.set_lu() ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    
  }
  // Cas ou le calcul a deja ete effectue :
  else
    return *tab[indice] ;
}


	     	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      

Matrice nondeg_pvect_r(const Matrice &lap, int l, double echelle, int puis, int base_r)
{

		// Routines de derivation
    static Matrice (*nondeg_pvect_r[MAX_BASE])(const Matrice&, int, double, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    nondeg_pvect_r[i] = _nondeg_pvect_r_pas_prevu ;
	}
		// Les routines existantes
	nondeg_pvect_r[R_CHEB >> TRA_R] = _nondeg_pvect_r_cheb ;
	nondeg_pvect_r[R_CHEBU >> TRA_R] = _nondeg_pvect_r_chebu ;
	nondeg_pvect_r[R_CHEBP >> TRA_R] = _nondeg_pvect_r_chebp ;
	nondeg_pvect_r[R_CHEBI >> TRA_R] = _nondeg_pvect_r_chebi ;
    }
    
    Matrice res(nondeg_pvect_r[base_r](lap, l, echelle, puis)) ;
    return res ;
}

