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

char ope_poisson_2d_solh_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <math.h>
#include <stdlib.h>

#include "proto.h"
#include "ope_elementary.h"


		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solh_poisson_2d_pas_prevu (int, int, double, double, Tbl&) {

    cout << " Solution homogene pas prevue ..... : "<< endl ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solh_poisson_2d_r_cheb (int n, int l, double alpha, double beta, 
			     Tbl& val_lim) {
                

  double echelle = beta / alpha ;

  // CASE l != 0
  if (l > 0) {
    val_lim.set(0,0) = pow(echelle-1, double(l)) ;
    val_lim.set(0,1) = double(l) * pow(echelle-1, double(l-1))/alpha ;
    val_lim.set(0,2) = pow(echelle+1, double(l)) ;
    val_lim.set(0,3) = double(l) * pow(echelle+1, double(l-1))/alpha ;
    
    
    val_lim.set(1,0) = pow(echelle-1, -double(l)) ;
    val_lim.set(1,1) = -double(l) * pow(echelle-1, -double(l+1))/alpha ;
    val_lim.set(1,2) = pow(echelle+1, -double(l)) ;
    val_lim.set(1,3) = -double(l) * pow(echelle+1, -double(l+1))/alpha ;
  }   
  // CASE l =0
  else {
    val_lim.set(0,0) = 1. ;
    val_lim.set(0,1) = 0. ;
    val_lim.set(0,2) = 1. ;
    val_lim.set(0,3) = 0. ;
  
    val_lim.set(1,0) = log(echelle-1) ;
    val_lim.set(1,1) = 1./(echelle-1)/alpha ;
    val_lim.set(1,2) = log(echelle+1) ;
    val_lim.set(1,3) = 1./(echelle+1)/alpha ;
  }

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
      cout << "_solh_poisson_2d_r_cheb : trop de Tbl" << endl ;
      abort() ;
      exit (-1) ;
    }
    
    
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
       
    Tbl res(2, n) ;
    res.set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    //Construction de la premiere solution homogene :
    // cad celle polynomiale.
   
    for (int i=0 ; i<n ; i++)
      if (l!=0)
	coloc[i] = pow(echelle-cos(M_PI*i/(n-1)), double(l)) ;
      else
	coloc[i] = 1 ;
    
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
      res.set(0, i) = coloc[i] ;
    
    // construction de la seconde solution homogene :
    // cad celle fractionnelle (ou log dans le cas l==0)
    for (int i=0 ; i<n ; i++)
      if (l != 0)
	coloc[i] = pow(echelle-cos(M_PI*i/(n-1)), -double(l)) ;
      else
	coloc[i] = log(echelle-cos(M_PI*i/(n-1))) ;

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

Tbl _solh_poisson_2d_r_chebp (int n, int l, double alpha, 
			      double, Tbl& val_lim) {
  
  val_lim.set(0,0) = (l!=0) ? 1 : 0 ;
  val_lim.set(0,1) = (l!=1) ? 0 : 1 ;
  val_lim.set(0,2) = 1. ;
  val_lim.set(0,3) = double(l)/alpha ;

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
	   cout << "_solh_poisson_2d_r_chebp : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
   	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    assert (div(l, 2).rem ==0) ;
    
    Tbl res(n) ;
    res.set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    for (int i=0 ; i<n ; i++)
      if (l!=0)
	coloc[i] = pow(sin(M_PI*i/2/(n-1)), double(l)) ;
      else
	coloc[i] = 1 ;

    cfrchebp(deg, deg, coloc, deg, coloc) ;
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
	       //--  R_CHEBI   -----
	      //-------------------
	
Tbl _solh_poisson_2d_r_chebi (int n, int l, 
			      double alpha, double, Tbl& val_lim) {

  val_lim.set(0,0) = 0 ;
  val_lim.set(0,1) = (l!=1) ? 0 : 1 ;
  val_lim.set(0,2) = 1. ;
  val_lim.set(0,3) = double(l)/alpha ;

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
	   cout << "_solh_r_chebi : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
    assert (div(l, 2).rem == 1)  ;

    Tbl res(n) ;
    res.set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    for (int i=0 ; i<n ; i++)
      coloc[i] = pow(sin(M_PI*i/2/(n-1)), double(l)) ;
    
    cfrchebi(deg, deg, coloc, deg, coloc) ;
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
	       //--  R_CHEBU   -----
	      //-------------------
	
Tbl _solh_poisson_2d_r_chebu (int n, int l, double alpha, 
			      double, Tbl& val_lim) {
    

  if (l==0) {
    cout << "Case l=0 in 2D Poisson not defined in the external compactified domain..." << endl ;
    abort() ;
  }

  val_lim.set(0,0) = pow(-2., double(l)) ;
  val_lim.set(0,1) = -double(l)*alpha*pow(-2, double(l+1.)) ;
  val_lim.set(0,2) = 0 ;
  val_lim.set(0,3) = 0 ;
 
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
	   cout << "_solh_poisson_2d_r_chebu : trop de Tbl" << endl ;
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
	coloc[i] = pow(-1-cos(M_PI*i/(n-1)), double(l)) ;
	
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


Tbl Ope_poisson_2d::get_solh () const {

  // Routines de derivation
  static Tbl (*solh_poisson_2d[MAX_BASE]) (int, int, double, double, Tbl&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solh_poisson_2d[i] = _solh_poisson_2d_pas_prevu ;
    }
    // Les routines existantes
    solh_poisson_2d[R_CHEB >> TRA_R] = _solh_poisson_2d_r_cheb ;
    solh_poisson_2d[R_CHEBP >> TRA_R] = _solh_poisson_2d_r_chebp ;
    solh_poisson_2d[R_CHEBI >> TRA_R] = _solh_poisson_2d_r_chebi ;
    solh_poisson_2d[R_CHEBU >> TRA_R] = _solh_poisson_2d_r_chebu ;
  }
  
  Tbl val_lim (2,4) ;
  val_lim.set_etat_qcq() ;
  Tbl res(solh_poisson_2d[base_r](nr,l_quant, alpha, beta, val_lim)) ;

  s_one_minus  = val_lim(0,0) ;
  ds_one_minus = val_lim(0,1) ;
  s_one_plus   = val_lim(0,2) ;
  ds_one_plus  = val_lim(0,3) ;

  s_two_minus  = val_lim(1,0) ;
  ds_two_minus = val_lim(1,1) ;
  s_two_plus   = val_lim(1,2) ;
  ds_two_plus  = val_lim(1,3) ;

  return res ;
}
