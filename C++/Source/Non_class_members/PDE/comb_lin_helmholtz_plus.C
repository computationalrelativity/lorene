/*
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


char comb_lin_helmholtz_plusC[] = "$Header $" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/12/11 14:48:49  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */

//fichiers includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrice.h"
#include "tbl.h"
#include "type_parite.h"
#include "indent.h"
#include "proto.h"


// Version Matrice --> Matrice
Matrice _cl_helmholtz_plus_pas_prevu (const Matrice& so, double, double, 
				      double) {
  cout << "CL Helmholtz plus not implemented" << endl ;
    abort() ;
    exit(-1) ;
    return so;
}


		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Matrice _cl_helmholtz_plus_r_cheb (const Matrice& source, double alpha, 
				    double beta, double masse) {

 
  int n = source.get_dim(0) ;
  assert (n==source.get_dim(1)) ;
  
  const int nmax = 10 ;// Nombre de Matrices stockees
  static Matrice* tab[nmax] ;  // les matrices calculees
  static int nb_dejafait = 0 ; // nbre de matrices calculees
  static double masse_dejafait[nmax] ;
  
  static double vieux_alpha = 0;
  static double vieux_beta = 0 ;
  static int vieux_n = 0;
  
  // Si on a change alpha ou n, on detruit tout :
  if ((vieux_alpha != alpha) || (vieux_n != n) || (vieux_beta != beta)) {
    for (int i=0 ; i<nb_dejafait ; i++) {
      masse_dejafait[i] = 0 ;
      delete tab[i] ;
    }
    
    nb_dejafait = 0 ;
    vieux_alpha = alpha ;
    vieux_beta = beta ;
    vieux_n = n ;
  }
  
  int indice = -1 ;
  
  // On determine si la matrice a deja ete calculee :
  for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if (masse_dejafait[conte] == masse)
      indice = conte ;
  
  // Calcul a faire : 
  if (indice  == -1) {
    if (nb_dejafait >= nmax) {
      cout << "_cl_helmholtz_plus_r_cheb : trop de matrices" << endl ;
      abort() ;
      exit (-1) ;
    }

    masse_dejafait[nb_dejafait] = masse ;
       
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


                //-------------------------
	       //- La routine a appeler ---
	      //---------------------------

Matrice cl_helmholtz_plus (const Matrice &source, double alpha, double beta, 
			    double masse, int base_r) {
    
		// Routines de derivation
    static Matrice (*cl_helmholtz_plus[MAX_BASE]) (const Matrice &, double, 
						    double, double) ;
    static int nap = 0 ;
    
    // Premier appel
    if (nap==0) {
      nap = 1 ;
      for (int i=0 ; i<MAX_BASE ; i++) {
	cl_helmholtz_plus[i] = _cl_helmholtz_plus_pas_prevu ;
	}
      // Les routines existantes
      cl_helmholtz_plus[R_CHEB >> TRA_R] = _cl_helmholtz_plus_r_cheb ;
    }
    
    Matrice res(cl_helmholtz_plus[base_r](source, alpha, beta, masse)) ;
    return res ;
}


//************************ TBL Versions *************************************




Tbl _cl_helmholtz_plus_pas_prevu (const Tbl &so) {

  cout << "Linear combination for Helmholtz plus not implemented..." << endl ;
  abort() ;
  exit(-1) ;
  return so;
}

               //-------------------
	       //--  R_CHEB  -------
	      //--------------------
Tbl _cl_helmholtz_plus_r_cheb (const Tbl& source) {
  
  int n = source.get_dim(0) ;

  Tbl barre(source) ;
  int dirac = 1 ;
  for (int i=0 ; i<n-2 ; i++) {
    barre.set(i) = ((1+dirac)*source(i)-source(i+2))
      /(i+1) ;
    if (i==0) dirac = 0 ;
  }
  
  Tbl res(barre) ;
  for (int i=0 ; i<n-4 ; i++)
    res.set(i) = barre(i)-barre(i+2) ;

  return res ;
}
		//----------------------------
	       //- Routine a appeler        ---
	      //------------------------------

Tbl cl_helmholtz_plus (const Tbl &source, int base_r) {
    
  // Routines de derivation
  static Tbl (*cl_helmholtz_plus[MAX_BASE])(const Tbl &) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      cl_helmholtz_plus[i] = _cl_helmholtz_plus_pas_prevu ;
    }
    // Les routines existantes
    cl_helmholtz_plus[R_CHEB >> TRA_R] = _cl_helmholtz_plus_r_cheb ;
  }
    
    Tbl res(cl_helmholtz_plus[base_r](source)) ;
    return res ;
}
