/*
 *  Methods for computing the homogeneous solutions for Ope_pois_vect_r.
 *
 *    (see file Ope_elementary for documentation).
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

char sh_pvect_r_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:53:34  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:16:13  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2004/12/23 16:00:10  j_novak
 * Modif. comments.
 *
 * Revision 1.2  2004/05/10 15:36:42  j_novak
 * Corrected a missing #include
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
#include <cstdlib>
#include <cmath>

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"

/*
 * 
 * Renvoie une ou 2 solution homogene
 * Si l <> 0 :
 *  Si base_r = R_CHEB deux solutions (x+echelle)^(l-1) dans (0, *) et
 *				1/(x+echelle)^(l+2) dans (1, *)
 *  Si base_r = R_CHEBU 1 solution (x-1)^(l+2) dans (*) 
 *  Si base_r = R_CHEBP ou R_CHEBI x^(l-1) dans (*)
 * 
 * Si l = 0:
 *  Si base_r = R_CHEB deux solutions (x+echelle) dans (0, *) et
 *				1/(x+echelle)^(l+2) dans (1, *)
 *  Si base_r = R_CHEBU 1 solution (x-1)^(l+2) dans (*) 
 *  Si base_r = R_CHEBP ou R_CHEBI x dans (*)
 *
 * Entree : 
 *	n : nbre de points en r
 *	l : nbre quantique associe
 *	echelle : cf ci-dessus, utile que dans le cas R_CHEB
 *	base_r : base de decomposition
 * 
 * Sortie :
 *	Tbl contenant les coefficient de la ou des solution homogenes
 * 
 */

namespace Lorene {
Tbl _sh_pvect_r_pas_prevu (int, int, double) ;
Tbl _sh_pvect_r_cheb (int, int, double) ;
Tbl _sh_pvect_r_chebp (int, int, double) ;
Tbl _sh_pvect_r_chebi (int, int, double) ;
Tbl _sh_pvect_r_chebu (int, int, double) ;

		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _sh_pvect_r_pas_prevu (int n, int l, double echelle) {

    cout << " Solution homogene pas prevue ..... : "<< endl ;
    cout << " N : " << n << endl ;
    cout << " l : " << l << endl ;
    cout << " echelle : " << echelle << endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _sh_pvect_r_cheb (int n, int l, double echelle) {
                
   const int nmax = 200 ; // Nombre de Tbl stockes
   static Tbl* tab[nmax] ;  // les Tbl calcules
   static int nb_dejafait = 0 ; // nbre de Tbl calcules
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
	   cout << "_sh_pvect_r_cheb : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    	
  //  assert (l < n) ;
    
    Tbl res(2, n) ;
    res.set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    //Construction de la premiere solution homogene :
    // cad celle polynomiale.
    
    if (l==0) {
	for (int i=0 ; i<n ; i++)
	    coloc[i] = echelle - cos(M_PI*i/(n-1)) ;
	
	cfrcheb(deg, deg, coloc, deg, coloc) ;
	for (int i=0 ; i<n ;i++)
	    res.set(0, i) = coloc[i] ;
    }
    else if (l==1) {
      res.set(0, 0) = 1 ;
      for (int i=1 ; i<n ; i++)
	res.set(0, i) = 0 ;
    }
    else {

	for (int i=0 ; i<n ; i++)
	    coloc[i] = pow(echelle-cos(M_PI*i/(n-1)), double(l-1)) ;
	
	cfrcheb(deg, deg, coloc, deg, coloc) ;
	for (int i=0 ; i<n ;i++)
	  res.set(0, i) = coloc[i] ;
    }
    
    
    // construction de la seconde solution homogene :
    // cad celle fractionnelle.
    for (int i=0 ; i<n ; i++)
      coloc[i] = 1/pow(echelle-cos(M_PI*i/(n-1)), double(l+2)) ;
	
      cfrcheb(deg, deg, coloc, deg, coloc) ;
      for (int i=0 ; i<n ;i++)
	res.set(1, i) = coloc[i] ;	
    
    delete [] coloc ;
    delete [] deg ;
    tab[nb_dejafait] = new Tbl(res) ;
    nb_dejafait ++ ;
    return res ;
    }
    
    else return *tab[indice] ;
}	
	
		//-------------------
	       //--  R_CHEBP  ------
	      //-------------------

Tbl _sh_pvect_r_chebp (int n, int l, double) {
   
   const int nmax = 200 ; // Nombre de Tbl stockes
   static Tbl* tab[nmax] ;  // les Tbl calcules
   static int nb_dejafait = 0 ; // nbre de Tbl calcules
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
	   cout << "_sh_pvect_r_chebp : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
   	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Tbl res(n) ;
    res.set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    assert (div(l, 2).rem == 1)  ;
   if (l==1) {
      res.set(0) = 1 ;
      for (int i=1 ; i<n ; i++)
	res.set(i) = 0 ;
    }
    else {
      for (int i=0 ; i<n ; i++)
	coloc[i] = pow(sin(M_PI*i/2/(n-1)), double(l-1)) ;
	
	cfrchebp(deg, deg, coloc, deg, coloc) ;
	for (int i=0 ; i<n ;i++)
	    res.set(i) = coloc[i] ;
    }

    delete [] coloc ;
    delete [] deg ;
    tab[nb_dejafait] = new Tbl(res) ;
    nb_dejafait ++ ;
    return res ;
    }
    
    else return *tab[indice] ;
}
	
	
	      	//-------------------
	       //--  R_CHEBI   -----
	      //-------------------
	
Tbl _sh_pvect_r_chebi (int n, int l, double) {
       
   const int nmax = 200 ; // Nombre de Tbl stockes
   static Tbl* tab[nmax] ;  // les Tbl calcules
   static int nb_dejafait = 0 ; // nbre de Tbl calcules
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
	   cout << "_sh_pvect_r_chebi : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
    assert (div(l, 2).rem == 0)  ;
   // assert (l < 2*n) ;
    
    Tbl res(n) ;
    res.set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    if (l==0) {
	res.set(0) = 1 ;
	for (int i=1 ; i<n ; i++)
	    res.set(i) = 0 ;
	    }
    else {
	for (int i=0 ; i<n ; i++)
	    coloc[i] = pow(sin(M_PI*i/2/(n-1)), double(l-1)) ;
	
	cfrchebi(deg, deg, coloc, deg, coloc) ;
	for (int i=0 ; i<n ;i++)
	    res.set(i) = coloc[i] ;
    }
	
    delete [] coloc ;
    delete [] deg ;
    tab[nb_dejafait] = new Tbl(res) ;
    nb_dejafait ++ ;
    return res ;
    }
    
    else return *tab[indice] ;
}
	
	
	
	       	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
	
Tbl _sh_pvect_r_chebu (int n, int l, double) {
           
   const int nmax = 200 ; // Nombre de Tbl stockes
   static Tbl* tab[nmax] ;  // les Tbl calcules
   static int nb_dejafait = 0 ; // nbre de Tbl calcules
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
	   cout << "_sh_pvect_r_chebu : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
	  }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
  //  assert (l < n-1) ;
    Tbl res(n) ;
    res.set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    for (int i=0 ; i<n ; i++)
	coloc[i] = pow(-1-cos(M_PI*i/(n-1)), double(l+2)) ;
	
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
	res.set(i) = coloc[i] ;
	
    delete [] coloc ;
    delete [] deg ;
    tab[nb_dejafait] = new Tbl(res) ;
    nb_dejafait ++ ;
    return res ;
    }
    
    else return *tab[indice] ;
}
	
	
	
	
	      	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
	      
Tbl sh_pvect_r(int n, int l, double echelle, int base_r) {

		// Routines de derivation
    static Tbl (*sh_pvect_r[MAX_BASE])(int, int, double) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    sh_pvect_r[i] = _sh_pvect_r_pas_prevu ;
	}
		// Les routines existantes
	sh_pvect_r[R_CHEB >> TRA_R] = _sh_pvect_r_cheb ;
	sh_pvect_r[R_CHEBU >> TRA_R] = _sh_pvect_r_chebu ;
	sh_pvect_r[R_CHEBP >> TRA_R] = _sh_pvect_r_chebp ;
	sh_pvect_r[R_CHEBI >> TRA_R] = _sh_pvect_r_chebi ;
    }
    
    Tbl res(sh_pvect_r[base_r](n, l, echelle)) ;
    return res ;
}
}
