/*
 *   Copyright (c) 2003 Philippe Grandclement
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

char ope_poisson_2d_mat_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <math.h>
#include <stdlib.h>

#include "proto.h"
#include "ope_elementary.h"

		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _poisson_2d_mat_pas_prevu(int, int, double, double, int) {
    cout << "laplacien pas prevu..." << endl ;
    abort() ;
    exit(-1) ;
    Matrice res(1, 1) ;
    return res;
}


		   //-------------------------
		   //--   CAS R_CHEBP    -----
		   //--------------------------
		    

Matrice _poisson_2d_mat_r_chebp (int n, int l, double, double, int) {
   
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
	   cout << "_poisson_2d_mat_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       

    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

   double* vect  = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBP) ;
	
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j] ; 
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxdsdx_1d (n, &vect, R_CHEBP) ;
	    for (int j=0 ; j<n ; j++)
		xd.set(j, i) = vect[j] ;
	
	}
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEBP) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
   
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+xd-l*l*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    }
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}



		   //------------------------
		   //--   CAS R_CHEBI    ----
		   //------------------------
		    

Matrice _poisson_2d_mat_r_chebi (int n, int l, double, double, int) {
   
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
	   cout << "_poisson_2d_mat_r_chebi : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBI) ;  // appel dans le cas impair
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxdsdx_1d (n, &vect, R_CHEBI) ;
	    for (int j=0 ; j<n ; j++)
		xd.set(j, i) = vect[j] ;
	}
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEBI) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+xd-l*l*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}




		   //-------------------------
		   //--   CAS R_CHEBU    -----
		   //-------------------------

Matrice _poisson_2d_mat_r_chebu_deux(int,int) ;
Matrice _poisson_2d_mat_r_chebu_trois(int,int) ;
Matrice _poisson_2d_mat_r_chebu_quatre(int,int) ;

Matrice _poisson_2d_mat_r_chebu( int n, int l, double, double, int puis) {
    Matrice res(n, n) ;
    res.set_etat_qcq() ;
    switch (puis) {
	case 4 :
	    res = _poisson_2d_mat_r_chebu_quatre (n, l) ;
	    break ;
	case 3 :
	    res = _poisson_2d_mat_r_chebu_trois (n, l) ;
	    break ;
	case 2 :
	    res = _poisson_2d_mat_r_chebu_deux (n, l) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
    }
    return res ;
}
    
    // Cas ou dzpuis = 4
Matrice _poisson_2d_mat_r_chebu_quatre (int n, int l) {
        
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_poisson_2d_mat_r_chebu_quatre : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice dx (n,n) ;
    dx.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      sxdsdx_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
      for (int j=0 ; j<n ; j++)
	dx.set(j, i) = vect[j] ;
    }

    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEBU) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+dx-l*l*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

// Cas ou dzpuis =3 
Matrice _poisson_2d_mat_r_chebu_trois (int n, int l) {
        
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_poisson_2d_mat_r_chebu_trois : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice dx(n, n) ;
    dx.set_etat_qcq() ;  
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    double* auxi = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	
      for (int j=0 ; j<n ; j++)
	  vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
      mult_xm1_1d_cheb (n, vect, auxi) ;
      for (int j=0 ; j<n ; j++)
	dd.set(j, i) = auxi[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      dsdx_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
      for (int j=0 ; j<n ; j++)
	dx.set(j, i) = vect[j] ;
    }

 
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxm1_1d_cheb (n, vect) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    delete [] auxi ;
    
    Matrice res(n, n) ;
    res = dd+dx-l*l*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}


    //Cas ou dzpuis = 2
Matrice _poisson_2d_mat_r_chebu_deux (int n, int l) {
        
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_poisson_2d_mat_r_chebu_deux : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;

    Matrice dx (n,n) ;
    dx.set_etat_qcq() ;
    Matrice res(n, n) ;
    res.set_etat_qcq() ;


    double* vect = new double[n] ;
    
    double* x2vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
	mult2_xm1_1d_cheb (n, vect, x2vect) ; // multiplication par (x-1)^2
	for (int j=0 ; j<n ; j++)
	    res.set(j, i) = x2vect[j] ;
    }
     
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      dsdx_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
      mult_xm1_1d_cheb (n, vect, x2vect) ; // multiplication par (x-1)
      for (int j=0 ; j<n ; j++)
	dx.set(j, i) = x2vect[j] ;
    }

    delete [] vect ;
    delete [] x2vect ;

    res = res + dx ;

    for (int i=0 ; i<n ; i++)
	res.set(i, i) -= l*l ;
    
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

		   //-------------------------
		   //--   CAS R_CHEB    -----
		   //-----------------------
		    

Matrice _poisson_2d_mat_r_cheb (int n, int l, double alf, double bet, int) {
            
  double echelle = bet / alf ;

   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_poisson_2d_mat_r_cheb : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j]*echelle*echelle ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
	multx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	    dd.set(j, i) += 2*echelle*vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
	multx_1d (n, &vect, R_CHEB) ;
	multx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	    dd.set(j, i) += vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxdsdx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<n ; j++)
	    xd.set(j, i) = vect[j]*echelle ;
	}
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxdsdx_1d (n, &vect, R_CHEB) ;
	multx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	    xd.set(j, i) += vect[j] ;
	}
	   
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+xd-l*l*xx ;   
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;  
}


void Ope_poisson_2d::do_ope_mat() const {
  if (ope_mat != 0x0) 
    delete ope_mat ;

  // Routines de derivation
  static Matrice (*poisson_2d_mat[MAX_BASE])(int, int, double, double, int);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      poisson_2d_mat[i] = _poisson_2d_mat_pas_prevu ;
    }
    // Les routines existantes
    poisson_2d_mat[R_CHEB >> TRA_R] = _poisson_2d_mat_r_cheb ;
    poisson_2d_mat[R_CHEBP >> TRA_R] = _poisson_2d_mat_r_chebp ;
    poisson_2d_mat[R_CHEBI >> TRA_R] = _poisson_2d_mat_r_chebi ;
    poisson_2d_mat[R_CHEBU >> TRA_R] = _poisson_2d_mat_r_chebu ;
  }
  ope_mat = new Matrice(poisson_2d_mat[base_r](nr, l_quant, alpha, beta, dzpuis)) ;
}
