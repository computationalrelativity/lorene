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


char lap_cpt_mat_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  2000/03/16  16:23:08  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */




//fichiers includes
#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"

/*
 * Routine renvoyant la matrice de l'operateur (1-x^2)*Laplacien f = s
 * Pour  l != 1 le resultat est donne en s est donne en Chebyshev et
 * f en Gelerkin (T_i + T_{i+1} pour l pair et (2*i+3)T_i + (2*i+1)T_{i+1} pour 
 * l impair.
 * Pour l=1 pas de probleme de singularite on reste donc en Chebyshev.
 */
  

		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _lap_cpt_mat_pas_prevu(int n, int l) {
    cout << "laplacien * (1-r^2/R_0^2) pas prevu..." << endl ;
    cout << "n : " << n << endl ;
    cout << "l : " << l << endl ;
    abort() ;
    exit(-1) ;
    Matrice res(1, 1) ; // On ne passe jamais ici de toute facon !
    return res;
}


		   //-------------------------
		   //--   CAS R_CHEBP    -----
		   //--------------------------
		    

Matrice _lap_cpt_mat_r_chebp (int n, int l) {
 
   const int nmax = 100 ;// Nombre de Matrices stockees
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
	   cout << "_laplacien_nul_mat_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       

    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice res(n-1, n-1) ;
    res.set_etat_qcq() ;
    
    
    double* xdsdx = new double[n] ;
    double* x2d2sdx2 = new double[n] ;
    double* d2sdx2 = new double[n] ;
    double* sxdsdx = new double[n] ;
    double* sx2 = new double [n] ;
    
    for (int i=0 ; i< n-1 ; i++) {
    for (int j=0 ; j<n ; j++)
	    xdsdx[j] = 0 ;
    xdsdx[i] = 1 ;
    xdsdx[i+1] = 1 ;
    xdsdx_1d (n, &xdsdx, R_CHEBP) ;
	
	
    for (int j=0 ; j<n ; j++)
	x2d2sdx2[j] = 0 ;
    x2d2sdx2[i] = 1 ;
    x2d2sdx2[i+1] = 1 ;
	
    d2sdx2_1d(n, &x2d2sdx2, R_CHEBP) ;
    for (int j=0 ; j<n ; j++)
	d2sdx2[j] = x2d2sdx2[j] ;
    multx2_1d(n, &x2d2sdx2, R_CHEBP) ;
	

    for (int j=0 ; j<n ; j++)
	sxdsdx[j] = 0 ;
    sxdsdx[i] = 1 ;
    sxdsdx[i+1] = 1 ;
    sxdsdx_1d(n, &sxdsdx, R_CHEBP) ;
	
	
    for (int j=0 ; j<n ; j++)
	sx2[j] = 0 ;
    sx2[i] = 1 ;
    sx2[i+1] = 1 ;
    sx2_1d(n, &sx2, R_CHEBP) ;
	
    for (int j=0 ; j<n-1 ; j++)
	res.set(j, i) = (d2sdx2[j] + 2*sxdsdx[j] - l*(l+1)*sx2[j])
		      - (x2d2sdx2[j]+2*xdsdx[j]) ;
    res.set(i, i) += l*(l+1) ;
    if (i < n-2)
	res.set(i+1, i) += l*(l+1) ;
    }
    delete [] d2sdx2 ;
    delete [] x2d2sdx2 ;
    delete [] sxdsdx ;
    delete [] xdsdx ;
    delete [] sx2 ;
    
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
		    

Matrice _lap_cpt_mat_r_chebi (int n, int l) {
  
   const int nmax = 100 ;// Nombre de Matrices stockees
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
	   cout << "_laplacien_nul_mat_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       

    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    // Non degenere si l = 1
    int taille = (l == 1) ? n : n-1 ;
    Matrice res(taille, taille) ;
    res.set_etat_qcq() ;
    
    
    double* xdsdx = new double[n] ;
    double* x2d2sdx2 = new double[n] ;
    double* d2sdx2 = new double[n] ;
    double* sxdsdx = new double[n] ;
    double* sx2 = new double [n] ;
    
    int f_un, f_deux ;
    
    for (int i=0 ; i<taille ; i++) {
    
    // Gelerkin ou Chebyshev ????
    if (taille == n) {
	f_un = 1 ;
	f_deux = 0 ;
	}
    else {
	f_un = 2*i+3 ;
	f_deux = 2*i+1 ;
	}
	
    for (int j=0 ; j<n ; j++)
	    xdsdx[j] = 0 ;
    xdsdx[i] = f_un ;
    if (i+1 < n)
	xdsdx[i+1] = f_deux ;
    xdsdx_1d (n, &xdsdx, R_CHEBI) ;
	
	
    for (int j=0 ; j<n ; j++)
	x2d2sdx2[j] = 0 ;
    x2d2sdx2[i] = f_un ;
    if (i+1 < n)
	x2d2sdx2[i+1] = f_deux ;
	
    d2sdx2_1d(n, &x2d2sdx2, R_CHEBI) ;
    for (int j=0 ; j<n ; j++)
	d2sdx2[j] = x2d2sdx2[j] ;
    multx2_1d(n, &x2d2sdx2, R_CHEBI) ;
	

    for (int j=0 ; j<n ; j++)
	sxdsdx[j] = 0 ;
    sxdsdx[i] = f_un ;
    if (i+1 < n)
	sxdsdx[i+1] = f_deux ;
    sxdsdx_1d(n, &sxdsdx, R_CHEBI) ;
	
	
    for (int j=0 ; j<n ; j++)
	sx2[j] = 0 ;
    sx2[i] = f_un ;
    if (i+1 < n)
	sx2[i+1] = f_deux ;
    sx2_1d(n, &sx2, R_CHEBI) ;
	
    for (int j=0 ; j<taille ; j++)
	res.set(j, i) = (d2sdx2[j] + 2*sxdsdx[j] - l*(l+1)*sx2[j])
		      - (x2d2sdx2[j]+2*xdsdx[j]) ;
    res.set(i, i) += l*(l+1)*f_un ;
    if (i < taille-1)
	res.set(i+1, i) += l*(l+1)*f_deux ;
    }
    delete [] d2sdx2 ;
    delete [] x2d2sdx2 ;
    delete [] sxdsdx ;
    delete [] xdsdx ;
    delete [] sx2 ;
    
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
Matrice lap_cpt_mat(int n, int l, int base_r)
{

		// Routines de derivation
    static Matrice (*lap_cpt_mat[MAX_BASE])(int, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    lap_cpt_mat[i] = _lap_cpt_mat_pas_prevu ;
	}
		// Les routines existantes
	lap_cpt_mat[R_CHEBP >> TRA_R] = _lap_cpt_mat_r_chebp ;
	lap_cpt_mat[R_CHEBI >> TRA_R] = _lap_cpt_mat_r_chebi ;
    }
    
    Matrice res(lap_cpt_mat[base_r](n, l)) ;
    return res ;
}

