/*
 *  Methods for linear combinations for Ope_pois_tens_rr
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

char cl_ptens_rr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:34  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:16:13  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/23 16:30:15  j_novak
 * New files and class for the solution of the rr component of the tensor Poisson
 * equation.
 *
 *
 *
 * $Header$
 *
 */

//fichiers includes
#include <cstdlib>

#include "matrice.h"
#include "type_parite.h"

/*   FONCTIONS FAISANT DES COMBINAISONS LINEAIRES DES LIGNES.
 * 
 * Entree :
 *  La Matrice de l'operateur
 *  l : nbre quantique
 *  puis = puissance dans la ZEC
 *  La base de developpement en R
 * 
 * Sortie :
 *  Renvoie la matrice apres Combinaison lineaire
 * 
 */

namespace Lorene {
Matrice _cl_ptens_rr_pas_prevu (const Matrice&, int, double, int) ;
Matrice _cl_ptens_rr_cheb (const Matrice&, int, double, int) ;
Matrice _cl_ptens_rr_chebi (const Matrice&, int, double, int) ;
Matrice _cl_ptens_rr_chebu (const Matrice&, int, double, int) ;
Matrice _cl_ptens_rr_chebp (const Matrice&, int, double, int) ;

// Version Matrice --> Matrice
Matrice _cl_ptens_rr_pas_prevu (const Matrice &source, int l, double echelle, int puis) {
    cout << "Combinaison lineaire pas prevu..." << endl ;
    cout << "Source : " << source << endl ;
    cout << "l : " << l << endl ;
    cout << "dzpuis : " << puis << endl ;
    cout << "Echelle : " << echelle << endl ;
    abort() ;
    exit(-1) ;
    return source;
}


		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Matrice _cl_ptens_rr_cheb (const Matrice &source, int l, double echelle, int) {
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
	   cout << "_cl_ptens_rr_cheb : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=i ; j<(n>(i+7)? i+7 : n) ; j++)
	    barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j))
				/(i+1) ;
	if (i==0) dirac = 0 ;
    }
    
    Matrice res(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=i ; j<(n>(i+5)? i+5 : n) ; j++)
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


Matrice _cl_ptens_rr_chebp (const Matrice &source, int l, double, int) {
    
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
	   cout << "_cl_ptens_rr_chebp : trop de matrices" << endl ;
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


Matrice _cl_ptens_rr_chebi (const Matrice &source, int l, double, int) {
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
	   cout << "_cl_ptens_rr_chebi : trop de matrices" << endl ;
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

Matrice _cl_ptens_rr_chebu (const Matrice &source, int l, double, int puis) {
  int n = source.get_dim(0) ;
  assert (n == source.get_dim(1)) ;
  if (puis != 4) {
    cout << "_ope_ptens_rr_mat_r_chebu : only the case dzpuis = 4 "
	 << '\n' << "is implemented! \n"
	 << "dzpuis = " << puis << endl ;
    abort() ;
  }

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
      cout << "_cl_ptens_rr_chebu_quatre : trop de matrices" << endl ;
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


		//-------------------------
	       //- La routine a appeler ---
	      //---------------------------

Matrice cl_ptens_rr(const Matrice &source, int l, double echelle, 
		   int puis, int base_r) {
    
  // Routines de derivation
  static Matrice (*combinaison[MAX_BASE])(const Matrice &, int, double, int) ;
  static int nap = 0 ;

  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      combinaison[i] = _cl_ptens_rr_pas_prevu ;
    }
    // Les routines existantes
    combinaison[R_CHEB >> TRA_R] = _cl_ptens_rr_cheb ;
    combinaison[R_CHEBU >> TRA_R] = _cl_ptens_rr_chebu ;
    combinaison[R_CHEBP >> TRA_R] = _cl_ptens_rr_chebp ;
    combinaison[R_CHEBI >> TRA_R] = _cl_ptens_rr_chebi ;
  }
  
  Matrice res(combinaison[base_r](source, l, echelle, puis)) ;
  return res ;
}

}
