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

char ope_poisson_2d_cl_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "proto.h"
#include "ope_elementary.h"

// Version Matrice --> Matrice
namespace Lorene {
Matrice _cl_poisson_2d_pas_prevu (const Matrice & source, int, double,int) {
    cout << "Combinaison lineaire pas prevu..." << endl ;
    abort() ;
    exit(-1) ;
    return source;
}


		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Matrice _cl_poisson_2d_r_cheb (const Matrice &source, int l, double echelle, int) {
   int n = source.get_dim(0) ;assert (n == source.get_dim(1)) ;
                
   const int nmax = 100 ; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
   static double vieux_echelle = 0 ;
   
   // Si on a change l'echelle : on detruit tout :
   if (vieux_echelle != echelle) {
       for (int i=0 ; i<nb_dejafait ; i++) {
	   l_dejafait[i] = -1 ;
	   nr_dejafait[i] = -1 ;
	   delete tab[i] ;
       }
       nb_dejafait = 0 ;
       vieux_echelle = echelle ;
   }
      
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_cl_poisson_r_cheb : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	    barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j))
				/(i+1) ;
	if (i==0) dirac = 0 ;
    }
    
    Matrice res(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = barre(i, j)-barre(i+2, j) ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    }
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;  
}

		//-------------------
	       //--  R_CHEBP   -----
	      //-------------------


Matrice _cl_poisson_2d_r_chebp (const Matrice &source, int l, double, int) {
    
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
   const int nmax = 100 ; // Nombre de Matrices stockees
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
	   cout << "_cl_poisson_2d_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
  
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	    barre.set(i, j) = (1+dirac)*source(i, j)-source(i+2, j) ;
	if (i==0) dirac = 0 ;
    }

    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = barre(i, j)-barre(i+2, j) ;

    Matrice res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = tilde(i, j)-tilde(i+1, j) ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    }
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

		//-------------------
	       //--  R_CHEBI   -----
	      //-------------------


Matrice _cl_poisson_2d_r_chebi (const Matrice &source, int l, double, int) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
       
   const int nmax = 100 ; // Nombre de Matrices stockees
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
	   cout << "_cl_poisson_2d_r_chebi : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
   
    for (int i=0 ; i<n-2 ; i++)
	for (int j=0 ; j<n ; j++)
	    barre.set(i, j) = source(i, j)-source(i+2, j) ;

    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = barre(i, j)-barre(i+2, j) ;    

    Matrice res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = tilde(i, j)-tilde(i+1, j) ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}
		//-------------------
	       //--  R_CHEBU   -----
	      //-------------------

Matrice _cl_poisson_2d_r_chebu_quatre (const Matrice&, int) ;
Matrice _cl_poisson_2d_r_chebu_trois (const Matrice&, int) ;
Matrice _cl_poisson_2d_r_chebu_deux (const Matrice&, int) ;


Matrice _cl_poisson_2d_r_chebu (const Matrice &source, int l, double, int puis) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
    Matrice res(n, n) ;
    res.set_etat_qcq() ;
    
    switch (puis) {
	case 4 :
	    res = _cl_poisson_2d_r_chebu_quatre(source, l) ;
	    break ;
	case 3 :
	    res = _cl_poisson_2d_r_chebu_trois (source, l) ;
	    break ;
	case 2 :
	    res = _cl_poisson_2d_r_chebu_deux(source, l) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
    }
    
    return res ;
}


// Cas dzpuis = 4
Matrice _cl_poisson_2d_r_chebu_quatre (const Matrice &source, int l) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
            
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
	   cout << "_cl_poisson_2d_r_chebu_quatre : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
  
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	     barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j)) ;
	if (i==0) dirac = 0 ;
	}
    
    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = (barre(i, j)-barre(i+2, j)) ;
	    
    Matrice prime(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    prime.set(i, j) = (tilde(i, j)-tilde(i+1, j)) ;
    
    Matrice res(prime) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = (prime(i, j)-prime(i+2, j)) ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

// Cas dzpuis == 3
Matrice _cl_poisson_2d_r_chebu_trois (const Matrice &source, int l) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
            
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
	   cout << "_cl_poisson_2d_r_chebu_trois : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;

    Matrice barre(source) ;
  
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	     barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j)) ;
	if (i==0) dirac = 0 ;
	}
    
    Matrice tilde(barre) ;
    for (int i=0 ; i<n-2 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = (barre(i, j)-barre(i+2, j)) ;

    Matrice res(tilde) ;
    for (int i=0 ; i<n-2 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = (tilde(i, j)+tilde(i+1, j)) ;
    
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}


//Cas dzpuis == 2
Matrice _cl_poisson_2d_r_chebu_deux (const Matrice &source, int l) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
            
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
	   cout << "_cl_poisson_2d_r_chebu_deux : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
  
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	     barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j)) ;
	if (i==0) dirac = 0 ;
	}
   
    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = (barre(i, j)-barre(i+2, j)) ;
	    
    Matrice res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = (tilde(i, j)+tilde(i+1, j)) ;

    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

void Ope_poisson_2d::do_ope_cl() const {
  if (ope_mat == 0x0)
    do_ope_mat() ;
  
  if (ope_cl != 0x0)
    delete ope_cl ;
  
  // Routines de derivation
  static Matrice (*cl_poisson_2d[MAX_BASE])(const Matrice&, int, double, int);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      cl_poisson_2d[i] = _cl_poisson_2d_pas_prevu ;
    }
    // Les routines existantes
    cl_poisson_2d[R_CHEBP >> TRA_R] = _cl_poisson_2d_r_chebp ;
    cl_poisson_2d[R_CHEBI >> TRA_R] = _cl_poisson_2d_r_chebi ;
    cl_poisson_2d[R_CHEB >> TRA_R] = _cl_poisson_2d_r_cheb ;
    cl_poisson_2d[R_CHEBU >> TRA_R] = _cl_poisson_2d_r_chebu ;
  }
  ope_cl = new Matrice(cl_poisson_2d[base_r](*ope_mat, l_quant, beta/alpha, 
					     dzpuis)) ;
}


}
