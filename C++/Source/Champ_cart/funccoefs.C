/*
 * Methods of class FuncSpec to compute coefficients.
 *
 *  (see file funcspec.h for documentation)
 *
 */

/*
 *   Copyright (c) 2019 Jerome Novak
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


#include<cmath>
#include<exception>
#include "funcspec.h"

namespace Lorene {

  // Local prototypes
  //------------------
  void check_fft_size(int) ;
  void cheb_transform(double* t, int np, int ntrans, int nskip, int next) ;
  void cheb_inv_transform(double* t, int np, int ntrans, int nskip, int next) ;
  double* cheb_ini(const int n );
  fftw_plan prepare_fft(int n, TabSpec*& pg) ;
  fftw_plan back_fft(int n, TabSpec*& pg) ;

  namespace {
    const int nmax = 500 ; //Maximal number of FFT sizes 
    int nworked = 0 ;
    int nwork_back = 0 ;
    int nwork_sin = 0 ;
    TabSpec* tab_tab[nmax] ;
    fftw_plan plan_fft[nmax] ;
    int nb_fft[nmax] ;
    TabSpec* tab_back[nmax] ;
    fftw_plan plan_fft_back[nmax] ;
    int nb_fft_back[nmax] ;
    double* table_sin[nmax] ;
    int nb_sin[nmax] ;
  }

  void FuncSpec::compute_coefs() const {

    if (!coefs_up_to_date) {
      if (!values_up_to_date)
	throw(std::runtime_error("Values not up to date in FuncSpec::compute_coefs")) ;  
      coefs = values ; // Initialisation with the values

#ifndef NDEBUG
      try {
	check_fft_size(nx) ;
      }
      catch(std::invalid_argument& e) {
	cerr << "Error in FFT transform: " << e.what() << " in x-direction" << endl ;
	throw ;
      }
#endif
      cheb_transform(coefs.tableau, nx, ny*nz, ny*nz, 1) ;

      if (ny > 1) {
#ifndef NDEBUG
	try {
	  check_fft_size(ny) ;
	}
	catch(std::invalid_argument& e) {
	  cerr << "Error in FFT transform: " << e.what() << " in y-direction"
	       << endl ;
	  throw ;
	}
#endif
	cheb_transform(coefs.tableau, ny, nx*nz, nz, ny*nz) ;
      }

      if (nz > 1) {
#ifndef NDEBUG
	try {
	  check_fft_size(nz) ;
	}
	catch(std::invalid_argument& e) {
	  cerr << "Error in FFT transform: " << e.what() << " in z-direction"
	       << endl ;
	  throw ;
	}
#endif
	cheb_transform(coefs.tableau, nz, nx*ny, 1, nz) ;
      }
      coefs_up_to_date = true ;
      del_deriv() ;
    }
  
    return ;
  
  }

  void cheb_transform(double* t, int np, int ntrans, int nskip, int next) {

    int nm1 = np - 1 ;
    int nm1s2 = nm1 / 2 ;

    TabSpec* pg = 0x0 ;
    fftw_plan p = prepare_fft(nm1, pg) ;
    TabSpec& g = *pg ;

    double* sinp = cheb_ini(np) ;

    int ntot = ntrans*np ; // total number of points/coefs
    for (int n=0; n<ntrans; n++) {
      // Starting index
      int ind = (n*next) % ntot + (n*next/ntot) ;
      double* ff0 = t + ind ;

      // Valeur en psi=0 de la partie antisymetrique de F, F_ :
      double fmoins0 = 0.5 * ( ff0[0] - ff0[nm1*nskip] );

      // Fonction G(psi) = F+(psi) + F_(psi) sin(psi) 
      //---------------------------------------------
      for (int i = 1; i < nm1s2 ; i++ ) {
	// ... indice du pt symetrique de psi par rapport a pi/2:
	int isym = nm1 - i ; 
	// ... F+(psi)	
	double fp = 0.5 * ( ff0[i*nskip] + ff0[isym*nskip] ) ;	
	// ... F_(psi) sin(psi)
	double fms = 0.5 * ( ff0[i*nskip] - ff0[isym*nskip] ) * sinp[i] ; 
	g.set(i) = fp + fms ;
	g.set(isym) = fp - fms ;
      }
      //... cas particuliers:
      g.set(0) = 0.5 * ( ff0[0] + ff0[nm1*nskip] );
      g.set(nm1s2) = ff0[nm1s2*nskip];

      // Developpement de G en series de Fourier par une FFT
      //----------------------------------------------------

      fftw_execute(p) ;

      // Coefficients pairs du developmt. de Tchebyshev de f
      //----------------------------------------------------
      //  Ces coefficients sont egaux aux coefficients en cosinus du developpement
      //  de G en series de Fourier (le facteur 2/nm1 vient de la normalisation
      //  de fftw) :

      double fac = 2./double(nm1) ;
      ff0[0] = g(0) / double(nm1) ;
      for (int i=2; i<nm1; i += 2) ff0[i*nskip] = fac*g(i/2) ;
      ff0[nm1*nskip] = g(nm1s2) / double(nm1) ;

      // Coefficients impairs du developmt. de Tchebyshev de f
      //------------------------------------------------------
      // 1. Coef. c'_k (recurrence amorcee a partir de zero):
      //    NB: Le 4/nm1 en facteur de g(i) est du a la normalisation de fftw
      //  (si fftw donnait reellement les coef. en sinus, il faudrait le
      //   remplacer par un +2.) 
      fac *= -2. ;
      ff0[nskip] = 0 ;
      double som = 0;
      for (int i = 3; i < np; i += 2 ) {
	ff0[i*nskip] = ff0[(i-2)*nskip] + fac * g(nm1-i/2) ;
	som += ff0[i*nskip] ;
      }
    
      // 2. Calcul de c_1 :
      double c1 = - ( fmoins0 + som ) / nm1s2 ;
    
      // 3. Coef. c_k avec k impair:	
      ff0[nskip] = c1 ;
      for (int i = 3; i < np; i += 2 ) ff0[i*nskip] += c1 ;

    }

  }

  fftw_plan prepare_fft(int n, TabSpec*& pg) {
    int index = -1 ;
    for (int i=0; ((i<nworked) && (index<0)); i++) 
      if (nb_fft[i] == n) index = i ; //Has the plan already been estimated?

    if (index <0) { //New plan needed
      index = nworked ;
      if (index >= nmax) {
	throw overflow_error("prepare_fft:: too many plans!") ;
      }
      tab_tab[index] = new TabSpec(n) ;
      TabSpec& tab = (*tab_tab[index]) ;
      plan_fft[index] = 
	fftw_plan_r2r_1d(n, tab.tableau, tab.tableau, FFTW_R2HC, FFTW_ESTIMATE) ;
      nb_fft[index] = n ;
      nworked++ ;
    }
    pg = tab_tab[index] ;
    return plan_fft[index] ;
  }

  double* cheb_ini(const int n ) {

    // Ce nombre de points a-t-il deja ete utilise ?
    int indice = -1 ;
    int i ;
    for ( i=0 ; i < nwork_sin ; i++ ) {
      if ( nb_sin[i] == n ) indice = i ;
    }
  
    // Initialisation
    if (indice == -1) {		    /* Il faut une nouvelle initialisation */
      if ( nwork_sin >= nmax ) {
	throw(std::range_error("cheb_ini : nwork_sin >= nmax !")) ; 
      }
      indice = nwork_sin ; nwork_sin++ ; nb_sin[indice] = n ;
    
      int nm1s2 = (n-1) / 2 ;  		
      table_sin[indice] = new double[nm1s2] ; 
    
      double xx = M_PI / double(n-1);
      for ( i = 0; i < nm1s2 ; i++ ) {
	table_sin[indice][i] = std::sin( xx * i );
      }
    }
    // Valeurs de retour
    return table_sin[indice] ;

  }


  void FuncSpec::compute_values() const {

    if (!values_up_to_date) {
      if (!coefs_up_to_date)
	throw(std::runtime_error("Coefficients not up to date in FuncSpec::compute_values")) ;
      values = coefs ; // initialization with the coefs

#ifndef NDEBUG
      try {
	check_fft_size(nx) ;
      }
      catch(std::invalid_argument& e) {
	cerr << "Error in FFT inverse transform: " << e.what() << " in x-direction" << endl ;
	throw ;
      }
#endif
      cheb_inv_transform(values.tableau, nx, ny*nz, ny*nz, 1) ;

      if (ny > 1) {
#ifndef NDEBUG
	try {
	  check_fft_size(ny) ;
	}
	catch(std::invalid_argument& e) {
	  cerr << "Error in inverse FFT transform: " << e.what() << " in y-direction"
	       << endl ;
	  throw ;
	}
#endif
	cheb_inv_transform(values.tableau, ny, nx*nz, nz, ny*nz) ;
      }

      if (nz > 1) {
#ifndef NDEBUG
	try {
	  check_fft_size(nz) ;
	}
	catch(std::invalid_argument& e) {
	  cerr << "Error in inverse FFT transform: " << e.what() << " in z-direction"
	       << endl ;
	  throw ;
	}
#endif
	cheb_inv_transform(values.tableau, nz, nx*ny, 1, nz) ;
      }
      values_up_to_date = true ;
      del_deriv() ;
    }
  
    return ;
  
  }

  void cheb_inv_transform(double* t, int np, int ntrans, int nskip, int next) {

    int nm1 = np - 1 ;
    int nm1s2 = nm1 / 2 ;

    TabSpec* pg = 0x0 ;
    fftw_plan p = back_fft(nm1, pg) ;
    TabSpec& g = *pg ; 

    double* sinp = cheb_ini(np) ;

    int ntot = ntrans*np ; // total number of points/coefs
    for (int n=0; n<ntrans; n++) {
      // Starting index
      int ind = (n*next) % ntot + (n*next/ntot) ;
      double* cf0 = t + ind ;

      // Calcul des coefficients de Fourier de la fonction 
      //   G(psi) = F+(psi) + F_(psi) sin(psi)
      // en fonction des coefficients de Tchebyshev de f:

      // Coefficients impairs de G
      //--------------------------
 
      double c1 = cf0[nskip] ;

      double som = 0;
      cf0[nskip] = 0 ;
      for (int i = 3; i < np; i += 2 ) {
	cf0[i*nskip] = cf0[i*nskip] - c1 ;
	som += cf0[i*nskip] ;
      }	
      // Valeur en psi=0 de la partie antisymetrique de F, F_ :
      double fmoins0 = - nm1s2 * c1 - som ;

      // Coef. impairs de G
      // NB: le facteur -0.25 est du a la normalisation de fftw; si fftw
      //     donnait exactement les coef. des sinus, ce facteur serait +0.5.
      for (int i = 3; i < np; i += 2 ) {
	g.set(nm1-i/2) = -0.25 * ( cf0[i*nskip] - cf0[(i-2)*nskip] ) ;
      }

      // Coefficients pairs de G
      //------------------------
      //  Ces coefficients sont egaux aux coefficients pairs du developpement de
      //   f en polynomes de Tchebyshev.
      // NB: le facteur 0.5 est du a la normalisation de fftw; si fftw
      //     donnait exactement les coef. des cosinus, ce facteur serait 1.

      g.set(0) = cf0[0] ;
      for (int i=1; i<nm1s2; i ++ ) g.set(i) = 0.5 * cf0[2*i*nskip] ;	
      g.set(nm1s2) = cf0[nm1*nskip] ;

      // Transformation de Fourier inverse de G 
      //---------------------------------------

      // FFT inverse
      fftw_execute(p) ;

      // Valeurs de f deduites de celles de G
      //-------------------------------------

      for (int i = 1; i < nm1s2 ; i++ ) {
	// ... indice du pt symetrique de psi par rapport a pi/2:
	int isym = nm1 - i ; 
	
	double fp = .5 * ( g(i) + g(isym) ) ;
	double fm = .5 * ( g(i) - g(isym) ) / sinp[i] ;

	cf0[i*nskip] = fp + fm ;
	cf0[isym*nskip] = fp - fm ;
      }
	
      //... cas particuliers:
      cf0[0] = g(0) + fmoins0 ;
      cf0[nm1*nskip] = g(0) - fmoins0 ;
      cf0[nm1s2*nskip] = g(nm1s2) ;

    }

  }

  fftw_plan back_fft(int n, TabSpec*& pg) {
    int index = -1 ;
    for (int i=0; ((i<nwork_back) && (index<0)); i++) 
      if (nb_fft_back[i] == n) index = i ; //Has the plan already been estimated?

    if (index <0) { //New plan needed
      index = nwork_back ;
      if (index >= nmax) {
	throw overflow_error("back_fft : too many plans!") ;
      }
      tab_back[index] = new TabSpec(n) ;
      TabSpec& tab = (*tab_back[index]) ;
      plan_fft_back[index] = 
	fftw_plan_r2r_1d(n, tab.tableau, tab.tableau, FFTW_HC2R, FFTW_ESTIMATE) ;
      nb_fft_back[index] = n ;
      nwork_back++ ;
    }
    pg = tab_back[index] ;
    return plan_fft_back[index] ;
  }

  void check_fft_size(int n) {
    if (n < 5) 
      throw(std::invalid_argument("Too small number of points for FFT transform")) ;
    // Division by 2   
    int reste = n % 2 ; 
    if (reste != 1) 
      throw(std::invalid_argument("Should have odd number of points for FFT transform")) ;
  }
  
} // Namespace Lorene
