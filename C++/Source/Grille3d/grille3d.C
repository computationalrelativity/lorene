/*
 *  Methods of class Grille3d and derived classes
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char grille3d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2008/01/08 13:53:29  j_novak
 * Special treatment of the case nt=1.
 *
 * Revision 1.3  2007/12/11 15:28:13  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.2  2002/10/16 14:36:36  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.10  2001/05/26  14:50:46  eric
 * *** empty log message ***
 *
 * Revision 2.9  2001/05/26  13:25:59  eric
 * Ajout du membre g_twice (grille double pour le desaliasing)
 * Modif de la declaration de g_angu (pointeur mutable)
 *   g_twice et g_angu ne sont calcules que si necessaire (cad si
 *   on appelle la fonction get_twice() ou get_angu()).
 *
 * Revision 2.8  2000/03/22  13:38:51  eric
 * Remplacement des iendl par endl dans <<
 *
 * Revision 2.7  1999/10/12  15:04:29  eric
 * *** empty log message ***
 *
 * Revision 2.6  1999/10/12  15:03:30  eric
 * *** empty log message ***
 *
 * Revision 2.5  1999/09/30  14:58:16  eric
 * Operator!= declare const
 *
 * Revision 2.4  1999/09/30  14:12:04  eric
 * sauve declaree const.
 *
 * Revision 2.3  1999/09/30  12:52:52  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.2  1999/03/01  14:35:21  eric
 * Modif affichage (operator<<)
 *
 *
 * $Header$
 *
 */


// Fichiers include
// ----------------
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "grilles.h"
#include "type_parite.h"
#include "proto.h"

		    	//-------------//
		    	// Mono-grille //
		    	//-------------//

			// Classe de base

// Constructeur
//-------------
Grille3d::Grille3d(int nrs, int nts, int nps) 
    : nr(nrs), nt(nts), np(nps)
{
    x = new double[nr] ; 
    tet = new double[nt] ;
    phi = new double[np] ;
}
    
// Destructeur
//------------
Grille3d::~Grille3d() {
    delete [] x ; 
    delete [] tet ; 
    delete [] phi ; 
}


			// Cas rare + sans symetrie

// Constructeur
//-------------
Grille3d_r::Grille3d_r(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    // Cette grille n'a pas de sens si np = 1
    assert(nps != 1) ;	// ??

      type_r = RARE ; 		// echantillonnage radial rarefie en 0     
      type_t = NONSYM ;  	// echantillonnage en theta sur [0, pi] 
      type_p = NONSYM ;  	// echantillonnage en phi sur [0, 2 pi[

    // Partie radiale
    double xx = 0 ;
    if (nr>1) xx = M_PI/double(2*(nr-1)) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = sin(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(nt-1) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = 2.*M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_r::~Grille3d_r() { }	// ne fait rien (c'est le destructeur
				// de la classe de base, Grille3d, qui 
				// fait tout. 


// cas fin + sans symetrie

// Constructeur
//-------------
Grille3d_f::Grille3d_f(int nrs, int nts, int nps)
    : Grille3d(nrs, nts, nps)
{
    // Cette grille n'a pas de sens si np = 1
    assert(nps != 1) ;	// ??

    double xx ;

    type_r = FIN ;	// echantillonnage radial fin sur [r_min,r_max]   
    type_t = NONSYM ;	// echantillonnage en theta sur [0, pi] 
    type_p = NONSYM ;	// echantillonnage en phi sur [0, 2 pi[

    // Partie radiale
    xx = 0 ;
    if (nr>1)  xx = M_PI/double(nr-1) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = -cos(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(nt-1) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = 2.*M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_f::~Grille3d_f() { }	// ne fait rien (c'est le destructeur
				// de la classe de base, Grille3d, qui 
				// fait tout. 

// cas echantillonnage en 1/r + sans symetrie

// Constructeur
//-------------
Grille3d_i::Grille3d_i(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{
    // Cette grille n'a pas de sens si np = 1
    assert(nps != 1) ;	// ??

    double xx ;

    type_r = UNSURR ;	// echantillonnage radial en u=1/r (fin sur [u_min,u_max]) 
    type_t = NONSYM ;	// echantillonnage en theta sur [0,pi] 
    type_p = NONSYM ;	// echantillonnage en phi sur [0,pi[

    // Partie radiale
    xx = 0 ;
    if (nr>1)  xx = M_PI/double(nr-1) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = -cos(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(nt-1) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = 2.*M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_i::~Grille3d_i() { }	// ne fait rien (c'est le destructeur
				// de la classe de base, Grille3d, qui 
				// fait tout. 

// Cas rare + symetrie equatoriale 

// Constructeur
//-------------
Grille3d_req::Grille3d_req(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;

    type_r = RARE ;	// echantillonnage radial rarefie en 0     
    type_t = SYM ;	// echantillonnage en theta sur [0, pi/2] 
    type_p = NONSYM ;	// echantillonnage en phi sur [0, 2 pi[
     
    // Partie radiale
    xx = 0 ;
    if (nr>1) xx = M_PI/double(2*(nr-1)) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = sin(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(2*(nt-1)) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = 2.*M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_req::~Grille3d_req() { }	// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 


// cas fin + symetrie equatoriale 

// Constructeur
//-------------
Grille3d_feq::Grille3d_feq(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;

    type_r = FIN ;	// echantillonnage radial fin sur [r_min,r_max]   
    type_t = SYM ;	// echantillonnage en theta sur [0, pi/2] 
    type_p = NONSYM ;	// echantillonnage en phi sur [0, 2 pi[

    // Partie radiale
    xx = 0 ;
    if (nr>1)  xx = M_PI/double(nr-1) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = -cos(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(2*(nt-1)) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = 2.*M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_feq::~Grille3d_feq() { }	// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 

// cas echantillonnage en 1/r + symetrie equatoriale 

// Constructeur
//-------------
Grille3d_ieq::Grille3d_ieq(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;

    type_r = UNSURR ;	// echantillonnage radial en u=1/r (fin sur [u_min,u_max]) 
    type_t = SYM ;	// echantillonnage en theta sur [0,pi/2] 
    type_p = NONSYM ;	// echantillonnage en phi sur [0,pi[

    // Partie radiale
    xx = 0 ;
    if (nr>1)  xx = M_PI/double(nr-1) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = -cos(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(2*(nt-1)) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = 2.*M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_ieq::~Grille3d_ieq() { }	// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 

// Cas rare + 2 phi

// Constructeur
//-------------
Grille3d_r2p::Grille3d_r2p(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;

    type_r = RARE ;	// echantillonnage radial rarefie en 0     
    type_t = NONSYM ;	// echantillonnage en theta sur [0,pi] 
    type_p = SYM ;	// echantillonnage en phi sur [0,pi[
     
    // Partie radiale
    xx = 0 ;
    if (nr>1)     xx = M_PI/double(2*(nr-1)) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = sin(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(nt-1) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_r2p::~Grille3d_r2p() { }	// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 

// cas fin + 2 phi

// Constructeur
//-------------
Grille3d_f2p::Grille3d_f2p(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;

    type_r = FIN ;	// echantillonnage radial fin sur [r_min,r_max]   
    type_t = NONSYM ;	// echantillonnage en theta sur [0,pi] 
    type_p = SYM ;	// echantillonnage en phi sur [0,pi[

    // Partie radiale
    xx = 0 ;
    if (nr>1)     xx = M_PI/double(nr-1) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = -cos(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(nt-1) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_f2p::~Grille3d_f2p() { }	// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 

// cas echantillonnage en 1/r + 2 phi

// Constructeur
//-------------
Grille3d_i2p::Grille3d_i2p(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;

    type_r = UNSURR ;	// echantillonnage radial en u=1/r (fin sur [u_min,u_max]) 
    type_t = NONSYM ;	// echantillonnage en theta sur [0,pi] 
    type_p = SYM ;	// echantillonnage en phi sur [0,pi[

    // Partie radiale
    xx = 0 ;
    if (nr>1) xx = M_PI/double(nr-1) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = -cos(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(nt-1) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_i2p::~Grille3d_i2p() { }	// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 

// Cas rare supersymmetrique

// Constructeur
//-------------
Grille3d_rs::Grille3d_rs(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;

    type_r = RARE ;	// echantillonnage radial rarefie en 0     
    type_t = SYM ;	// echantillonnage en theta sur [0,pi/2] 
    type_p = SYM ;	// echantillonnage en phi sur [0,pi[
     
    // Partie radiale
    xx = 0 ;
    if (nr>1)     xx = M_PI/double(2*(nr-1)) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = sin(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(2*(nt-1)) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_rs::~Grille3d_rs() { }		// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 

// cas fin supersymetrique

// Constructeur
//-------------
Grille3d_fs::Grille3d_fs(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;

    type_r = FIN ;	// echantillonnage radial fin sur [r_min,r_max]   
    type_t = SYM ;	// echantillonnage en theta sur [0,pi/2] 
    type_p = SYM ;	// echantillonnage en phi sur [0,pi[

    // Partie radiale
    xx = 0 ;
    if (nr>1)     xx = M_PI/double(nr-1) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = -cos(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(2*(nt-1)) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_fs::~Grille3d_fs() { }		// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 

// cas echantillonnage en 1/r supersymetrique

// Constructeur
//-------------
Grille3d_is::Grille3d_is(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;

    type_r = UNSURR ;	// echantillonnage radial en u=1/r (fin sur [u_min,u_max]) 
    type_t = SYM ;	// echantillonnage en theta sur [0,pi/2] 
    type_p = SYM ;	// echantillonnage en phi sur [0,pi[

    // Partie radiale
    xx = 0 ;
    if (nr>1)     xx = M_PI/double(nr-1) ;
    for (int i=0 ; i<nr ; i++) {
	x[i] = -cos(xx*i) ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(2*(nt-1)) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_is::~Grille3d_is() { }		// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 


		


                     // Cas fin de jacobi + sans symetrie

// Constructeur
//-------------
Grille3d_fj::Grille3d_fj(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    // Cette grille n'a pas de sens si np = 1
    assert(nps != 1) ;	// ??

    double xx ;
    double* yy ;

      type_r = FINJAC ; 	// echantillonnage radial fin avec Jacobi     
      type_t = NONSYM ;  	// echantillonnage en theta sur [0, pi] 
      type_p = NONSYM ;  	// echantillonnage en phi sur [0, 2 pi[

    // Partie radiale
      yy = 0;
      if (nr>1) { yy = pointsgausslobatto(nr-1); }
      for (int i=0 ; i<nr ; i++) {
    	x[i] = yy[i] ;
    }

    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(nt-1) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = 2.*M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_fj::~Grille3d_fj() { }	// ne fait rien (c'est le destructeur
				// de la classe de base, Grille3d, qui 
				// fait tout. 



// cas fin de jacobi + symetrie equatoriale 

// Constructeur
//-------------
Grille3d_fjeq::Grille3d_fjeq(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;
    double* yy ;
    

    type_r = FINJAC ;	// echantillonnage radial fin avec Jacobi   
    type_t = SYM ;	// echantillonnage en theta sur [0, pi/2] 
    type_p = NONSYM ;	// echantillonnage en phi sur [0, 2 pi[

    // Partie radiale
    yy = 0 ;
    if (nr>1) yy = pointsgausslobatto(nr-1);
      for (int i=0 ; i<nr ; i++) {
    	x[i] = yy[i] ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(2*(nt-1)) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = 2.*M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_fjeq::~Grille3d_fjeq() { }	// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 



// cas fin de jacobi + 2 phi

// Constructeur
//-------------
Grille3d_fj2p::Grille3d_fj2p(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double* yy ;
    double xx ;

    type_r = FINJAC ;	// echantillonnage radial fin avec Jacobi   
    type_t = NONSYM ;	// echantillonnage en theta sur [0,pi] 
    type_p = SYM ;	// echantillonnage en phi sur [0,pi[

    // Partie radiale
    yy = 0 ;
    if (nr>1) { yy = pointsgausslobatto(nr-1); }
      for (int i=0 ; i<nr ; i++) {
    	x[i] = yy[i] ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(nt-1) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_fj2p::~Grille3d_fj2p() { }	// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 



// cas fin de jacobi supersymetrique

// Constructeur
//-------------
Grille3d_fjs::Grille3d_fjs(int nrs, int nts, int nps) 
    : Grille3d(nrs, nts, nps)
{

    double xx ;
    double* yy ;

    type_r = FINJAC ;	// echantillonnage radial fin avec Jacobi   
    type_t = SYM ;	// echantillonnage en theta sur [0,pi/2] 
    type_p = SYM ;	// echantillonnage en phi sur [0,pi[

    // Partie radiale
    yy = 0 ;
    if (nr>1) { yy = pointsgausslobatto(nr-1); }
      for (int i=0 ; i<nr ; i++) {
    	x[i] = yy[i] ;
    }
    // Partie en theta
    if (nt > 1) 
	xx = M_PI/double(2*(nt-1)) ;
    for (int i=0 ; i<nt ; i++) {
	tet[i] = xx*i ;
    }
    // Partie longitudinale
    xx = M_PI/double(np) ;
    for (int i=0 ; i<np ; i++) {
	phi[i] = xx*i ;
    }
}

// Destructeur 
//-------------
Grille3d_fjs::~Grille3d_fjs() { }	// ne fait rien (c'est le destructeur
					// de la classe de base, Grille3d, qui 
					// fait tout. 
