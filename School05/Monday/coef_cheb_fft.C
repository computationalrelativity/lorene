/*
 * Chebyshev transform via a FFT 
 *
 */

/*
 *   Copyright (c) 2005 Eric Gourgoulhon & Jerome Novak
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

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:09:47  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/11/14 14:12:10  e_gourgoulhon
 * Added include <assert.h>
 *
 * Revision 1.1  2005/11/14 01:56:58  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */


#include <iostream>

using namespace std ;

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <fftw3.h>

fftw_plan prepare_fft(int n, double*& pg) {

    static const int nmax = 50 ; //Maximal number of FFT sizes 
    static int nworked = 0 ;
    static double* tab_tab[nmax] ;
    static fftw_plan plan_fft[nmax] ;
    static int nb_fft[nmax] ;

  int index = -1 ;
  for (int i=0; ((i<nworked) && (index<0)); i++) 
    if (nb_fft[i] == n) index = i ; //Has the plan already been estimated?

  if (index <0) { //New plan needed
    index = nworked ;
    if (index >= nmax) {
      cout << "prepare_fft: " << endl ;
      cout << "too many plans!" << endl ;
      abort() ;
    }
        
    double* tab = new double[n] ;
    tab_tab[index] = tab ; 

    plan_fft[index] = 
      fftw_plan_r2r_1d(n, tab, tab, FFTW_R2HC, FFTW_MEASURE) ;
      
    nb_fft[index] = n ;
    nworked++ ;
  }
  assert((index>=0)&&(index<nmax)) ;
  pg = tab_tab[index] ;
  return plan_fft[index] ;
}


double* cheb_ini(int n) {

    const int nmax = 30 ; /* Nombre maximun de dimensions differentes */
    static	double*	table_sin[nmax] ;	/* Tableau des pointeurs sur les tableaux */
    static	int	nwork = 0 ;		/* Nombre de tableaux deja initialises */
    static	int	tbn[nmax] ;		/* Tableau des points deja initialises */

    // Ce nombre de points a-t-il deja ete utilise ?
    int indice = -1 ;
    for (int i=0 ; i < nwork ; i++ ) {
	    if ( tbn[i] == n ) indice = i ;
	}

    // Initialisation
    if (indice == -1) {		    /* Il faut une nouvelle initialisation */
	    if ( nwork >= nmax ) {
	        cout << "cheb_ini : nwork >= nmax !" << endl ; 
	        abort() ; 
	    }
	    indice = nwork ; nwork++ ; tbn[indice] = n ;

	    int nm1s2 = (n-1) / 2 ;  		
	    table_sin[indice] = new double[nm1s2] ; 

	    double xx = M_PI / double(n-1);
	    for (int i = 0; i < nm1s2 ; i++ ) {
	        table_sin[indice][i] = sin( xx * i );
	    }
	}

    return table_sin[indice] ;
}


void coef_cheb_fft(int nr, const double* ff, double* cf) {

    // Nombre de points pour la FFT:
    int nm1 = nr - 1;
    int nm1s2 = nm1 / 2;

    // Recherche des tables pour la FFT:
    double* pg = 0x0 ;
    fftw_plan p = prepare_fft(nm1, pg) ;
    double* g = pg ;

    // Recherche de la table des sin(psi) :
    double* sinp = cheb_ini(nr);		

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a x par  x = - cos(psi)   et F(psi) = f(x(psi)).  
 */
 
    // Valeur en psi=0 de la partie antisymetrique de F, F_ :
    double fmoins0 = 0.5 * ( ff[0] - ff[nm1] );

// Fonction G(psi) = F+(psi) + F_(psi) sin(psi) 
//---------------------------------------------
    	    for (int i = 1; i < nm1s2 ; i++ ) {
// ... indice du pt symetrique de psi par rapport a pi/2:
    int isym = nm1 - i ; 
// ... F+(psi)	
    double fp = 0.5 * ( ff[i] + ff[isym] ) ;	
// ... F_(psi) sin(psi)
    double fms = 0.5 * ( ff[i] - ff[isym] ) * sinp[i] ; 
    g[i] = fp + fms ;
    g[isym] = fp - fms ;
    	    }
//... cas particuliers:
    g[0] = 0.5 * ( ff[0] + ff[nm1] );
    g[nm1s2] = ff[nm1s2];

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------

    fftw_execute(p) ;

// Coefficients pairs du developmt. de Tchebyshev de f
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en cosinus du developpement
//  de G en series de Fourier (le facteur 2/nm1 vient de la normalisation
//  de fftw) :

	    double fac = 2./double(nm1) ;
	    cf[0] = g[0] / double(nm1) ;
	    for (int i=2; i<nm1; i += 2) cf[i] = fac*g[i/2] ;
	    cf[nm1] = g[nm1s2] / double(nm1) ;

// Coefficients impairs du developmt. de Tchebyshev de f
//------------------------------------------------------
// 1. Coef. c'_k (recurrence amorcee a partir de zero):
//    NB: Le 4/nm1 en facteur de g(i) est du a la normalisation de fftw
//  (si fftw donnait reellement les coef. en sinus, il faudrait le
//   remplacer par un +2.) 
	    fac *= -2. ;
    	    cf[1] = 0 ;
    	    double som = 0;
    	    for (int i = 3; i < nr; i += 2 ) {
	      cf[i] = cf[i-2] + fac * g[nm1-i/2] ;
	    	som += cf[i] ;
    	    }

// 2. Calcul de c_1 :
	    double c1 = - ( fmoins0 + som ) / nm1s2 ;

// 3. Coef. c_k avec k impair:	
    	    cf[1] = c1 ;
    	    for (int i = 3; i < nr; i += 2 ) cf[i] += c1 ;

}
