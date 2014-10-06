/*
 *   Copyright (c) 2004 Philippe Grandclement
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

char ope_poisson_2d_non_dege_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "proto.h"
#include "ope_elementary.h"


		//------------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _poisson_2d_non_dege_pas_prevu(const Matrice &lap, int, double, int) {
    cout << "Construction non degeneree pas prevue..." << endl ;
    abort() ;
    exit(-1) ;
    return lap ;
}



	     	//-------------------
	       //--  R_CHEB   -------
	      //--------------------

Matrice _poisson_2d_non_dege_r_cheb (const Matrice &lap, int l, double echelle, int) {
    
    
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
	   cout << "_poisson_2d_non_dege_r_cheb : trop de matrices" << endl ;
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
	      
Matrice _poisson_2d_non_dege_r_chebp (const Matrice &lap, int l, double, int) {
    
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
	   cout << "_poisson_2d_non_dege_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    assert (div(l, 2).rem == 0) ;
  //  assert (l<=2*n-2) ;
    
    if (l==0) {
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
	      
Matrice _poisson_2d_non_dege_r_chebi (const Matrice &lap, int l, double, int) {
    
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
	   cout << "_poisson_2d_non_dege_r_chebi : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
    assert (div(l, 2).rem == 1) ;
  //  assert (l<=2*n-1) ;
    
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
	       //--  R_CHEBU   -----
	      //-------------------
Matrice _poisson_2d_non_dege_r_chebu_quatre (const Matrice&, int) ;  
Matrice _poisson_2d_non_dege_r_chebu_trois (const Matrice&, int) ;  
Matrice _poisson_2d_non_dege_r_chebu_deux (const Matrice&, int) ;    
	      
Matrice _poisson_2d_non_dege_r_chebu (const Matrice &lap, int l, double, int puis) {

    switch (puis) {
	case 4 :
	    return _poisson_2d_non_dege_r_chebu_quatre (lap, l) ;
	case 3 :
	    return _poisson_2d_non_dege_r_chebu_trois (lap, l) ;
	case 2 :
	    return _poisson_2d_non_dege_r_chebu_deux (lap, l) ;
	default :
	    abort() ;
	    exit(-1) ;
	    return Matrice(0, 0) ;
    }
}

// Cas dzpuis = 4 ;
Matrice _poisson_2d_non_dege_r_chebu_quatre (const Matrice &lap, int l) {
    
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
	   cout << "_poisson_2d_non_dege_r_chebu : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
  //  assert (l<=n-2) ;
    
    if (l==0) {
	Matrice res(n-2, n-2) ;
	res.set_etat_qcq() ;
	for (int i=0 ; i<n-2 ; i++)
	    for (int j=0 ; j<n-2 ; j++)
		res.set(i, j) = lap(i, j+2) ;
	res.set_band(3, 0) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
    }
    else {
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
    }
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}
// Cas dzpuis = 3 ;
Matrice _poisson_2d_non_dege_r_chebu_trois (const Matrice &lap, int l) {
    
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
	   cout << "_poisson_2d_non_dege_r_chebu_trois : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice res(n-2, n-2) ;
    res.set_etat_qcq() ;
    for (int i=0 ; i<n-2 ; i++)
	for (int j=0 ; j<n-2 ; j++)
	    res.set(i, j) = lap(i, j+2) ;
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

// Cas dzpuis = 2 ;
Matrice _poisson_2d_non_dege_r_chebu_deux (const Matrice &lap, int l) {
    
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
	   cout << "_poisson_2d_non_dege_r_chebu : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
  //  assert (l<=n-2) ;
    
    if (l==0) {
	Matrice res(n-2, n-2) ;
	res.set_etat_qcq() ;
	for (int i=0 ;i<n-2 ; i++)
	    for (int j=0 ; j<n-2 ; j++)
		res.set(i, j) = lap(i, j+2) ;
	res.set_band(3, 2) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
    }
    else {
	Matrice res(n-1, n-1) ;
	res.set_etat_qcq() ;
	for (int i=0 ;i<n-1 ; i++)
	    for (int j=0 ; j<n-1 ; j++)
		res.set(i, j) = lap(i, j+1) ;
	res.set_band(4, 1) ;
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


void Ope_poisson_2d::do_non_dege() const {
  if (ope_cl == 0x0)
    do_ope_cl() ;
  
  if (non_dege != 0x0)
    delete non_dege ;
  
  // Routines de derivation
  static Matrice (*poisson_2d_non_dege[MAX_BASE])(const Matrice&, 
						    int, double,int);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      poisson_2d_non_dege[i] = _poisson_2d_non_dege_pas_prevu ;
    }
    // Les routines existantes
    poisson_2d_non_dege[R_CHEB >> TRA_R] = _poisson_2d_non_dege_r_cheb ;
    poisson_2d_non_dege[R_CHEBP >> TRA_R] = _poisson_2d_non_dege_r_chebp ;
    poisson_2d_non_dege[R_CHEBI >> TRA_R] = _poisson_2d_non_dege_r_chebi ;
    poisson_2d_non_dege[R_CHEBU >> TRA_R] = _poisson_2d_non_dege_r_chebu ;
  }
  non_dege = new Matrice(poisson_2d_non_dege[base_r](*ope_cl, l_quant, 
						       beta/alpha, dzpuis)) ;
}
