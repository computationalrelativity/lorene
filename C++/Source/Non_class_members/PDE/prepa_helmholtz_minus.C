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


char prepa_helmholtz_minus_C[] = "$Header$" ;

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




		//------------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _prepa_helmholtz_minus_nondege_pas_prevu(const Matrice &so, 
						double, double, double) {
  
    cout << "Unknown case for prepa_helmholtz_minus_nondege" << endl ;
    abort() ;
    exit(-1) ;
    return so;
}



	     	//-------------------
	       //--  R_CHEB   -------
	      //--------------------

Matrice _prepa_helmholtz_minus_nondege_r_cheb (const Matrice &lap, 
					       double alpha, double beta, 
					       double masse) {
    
  int n = lap.get_dim(0) ;

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
      cout << "_prepa_helmholtz_minus_nondege_r_cheb trop de matrices" << endl ;
      abort() ;
      exit (-1) ;
    }
    
    masse_dejafait[nb_dejafait] = masse ;
    
    
    int non_dege = 2 ;
    
    Matrice res(n-non_dege, n-non_dege) ;
    res.set_etat_qcq() ;
    for (int i=0 ; i<n-non_dege ; i++)
      for (int j=0 ; j<n-non_dege ; j++)
	res.set(i, j) = lap(i, j+non_dege) ;
    
    res.set_band (4,4) ;
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
	       //--  R_CHEBU  -------
	      //--------------------
Matrice _prepa_helmholtz_minus_nondege_r_chebu (const Matrice &lap, 
						double alpha, double, 
						double masse) {
    
  int n = lap.get_dim(0) ;

  const int nmax = 10 ;// Nombre de Matrices stockees
  static Matrice* tab[nmax] ;  // les matrices calculees
  static int nb_dejafait = 0 ; // nbre de matrices calculees
  static double masse_dejafait[nmax] ;
  
  static double vieux_alpha = 0;
  static int vieux_n = 0;
  
  // Si on a change alpha ou n, on detruit tout :
  if ((vieux_alpha != alpha) || (vieux_n != n)) {
    for (int i=0 ; i<nb_dejafait ; i++) {
      masse_dejafait[i] = 0 ;
      delete tab[i] ;
    }
    
    nb_dejafait = 0 ;
    vieux_alpha = alpha ;
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
      cout << "_prepa_helmholtz_minus_nondege_r_chebu : trop de matrices" << endl ;
      abort() ;
      exit (-1) ;
    }
    
    masse_dejafait[nb_dejafait] = masse ;
    
    
    int non_dege = 1 ;
    
    Matrice res(n-non_dege, n-non_dege) ;
    res.set_etat_qcq() ;
    for (int i=0 ; i<n-non_dege ; i++)
      for (int j=0 ; j<n-non_dege ; j++)
	res.set(i, j) = lap(i, j+non_dege) ;
    
    res.set_band (5,3) ;
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
	      
Matrice prepa_helmholtz_minus_nondege(const Matrice &ope, double alpha, 
				      double beta, double masse, int base_r) {

  // Routines de derivation
  static Matrice (*prepa_helmholtz_minus_nondege[MAX_BASE])
    (const Matrice&, double, double, double) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      prepa_helmholtz_minus_nondege[i] = 
	_prepa_helmholtz_minus_nondege_pas_prevu ;
    }
    // Les routines existantes
    prepa_helmholtz_minus_nondege[R_CHEB >> TRA_R] = 
      _prepa_helmholtz_minus_nondege_r_cheb ;
    prepa_helmholtz_minus_nondege[R_CHEBU >> TRA_R] = 
      _prepa_helmholtz_minus_nondege_r_chebu ;
  }
  
  Matrice res(prepa_helmholtz_minus_nondege[base_r](ope, alpha, beta, masse)) ;
  return res ;
}

