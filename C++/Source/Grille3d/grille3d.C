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
 * Revision 1.10  2014/10/13 08:52:59  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2013/06/07 14:44:33  j_novak
 * Coefficient computation for even Legendre basis.
 *
 * Revision 1.8  2013/06/06 15:31:32  j_novak
 * Functions to compute Legendre coefficients (not fully tested yet).
 *
 * Revision 1.7  2013/06/05 15:00:26  j_novak
 * Suppression of all classes derived from Grille3d. Now Grille3d is no
 * longer an abstract class. r-samplings are only one of RARE, FIN or
 * UNSURR (FINJAC has been removed). Instead, Mg3d possesses a new member
 * colloc_r[nzone] defining the type of collocation points (spectral
 * bases) in each domain.
 *
 * Revision 1.6  2008/08/27 08:47:38  jl_cornou
 * Added R_JACO02 case
 *
 * Revision 1.5  2008/01/09 14:04:03  j_novak
 * Initialization of xx
 *
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
#include <cmath>

#include "tbl.h"
#include "grilles.h"
#include "proto.h"


		    	//-------------//
		    	// Mono-grille //
		    	//-------------//

// Constructeur
//-------------
namespace Lorene {
Grille3d::Grille3d(int nrs, int nts, int nps, int typer, int typet,
		   int typep, int baser) 
  : nr(nrs), nt(nts), np(nps), type_r(typer), type_t(typet),
    type_p(typep), base_r(baser)
{

  //Radial part
  assert(nr > 0) ;
  x = new double[nr] ;
  x[0] = 0. ;
  if (nr > 1) compute_radial_grid() ;

  //Theta part
  assert(nt > 0) ;
  tet = new double[nt] ;
  double fac_tet = M_PI ;
  if (type_t == SYM) fac_tet *= 0.5 ;
  if (nt == 1) 
    fac_tet = 0 ;
  else 
    fac_tet /= double(nt-1) ;
  for (int i=0; i<nt; i++)
    tet[i] = double(i)*fac_tet ;
  if ( (type_t != SYM) && (type_t != NONSYM) ) {
    cout << "Grille3d: unknown type in theta!" << endl ;
    abort() ;
  }

  //Phi part
  assert(np > 0) ;
  phi = new double[np] ;
  double fac_phi = M_PI / double(np) ;
  if (type_p == NONSYM) fac_phi *= 2. ;
  for (int i=0; i<np; i++)
    phi[i] = double(i)*fac_phi ;
  if ( (type_p != SYM) && (type_p != NONSYM) ) {
    cout << "Grille3d: unknown type in phi!" << endl ;
    abort() ;
  }
  
}
    
// Destructeur
//------------
Grille3d::~Grille3d() {
    delete [] x ; 
    delete [] tet ; 
    delete [] phi ; 
}

void Grille3d::compute_radial_grid() {

  assert(nr > 1) ;
  double xx = 0 ;

  switch (base_r) {
  case BASE_CHEB :
    switch (type_r) {
    case RARE:
      xx = M_PI/double(2*(nr-1)) ;
      for (int i=0; i<nr; i++)
	x[i] = sin(xx*double(i)) ;
      break ;
    case FIN: case UNSURR :
      xx = M_PI/double(nr-1) ;
      for (int i=0 ; i<nr ; i++) 
	x[i] = -cos(xx*double(i)) ;
      break ;
    default:
      cout << "Grille3d::compute_radial_grid : " << endl ;
      cout << "Unknown type of sampling for the Chebyshev basis!" << endl ;
      abort() ;
    }
    break ;
  case BASE_LEG :
    switch (type_r) {
    case FIN:
      legendre_collocation_points(nr, x) ;
      break ;
    case RARE: {
      Tbl full_x(2*nr-1) ;
      full_x.set_etat_qcq() ;
      legendre_collocation_points(2*nr - 1, full_x.t) ;
      for (int i=0; i<nr; i++)
	x[i] = full_x(i+nr-1) ;
      break ;
    }
    default:
      cout << "Grille3d::compute_radial_grid : " << endl ;
      cout << "Unknown type of sampling for the Legendre basis!" << endl ;
      abort() ;
    }
    break ;
  case BASE_JAC02 : {
    double* yy = pointsgausslobatto(nr-1); 
    for (int i=0 ; i<nr ; i++) {
      x[i] = yy[i] ;
    }
    delete [] yy ;
    break ;
  }
  default : 
    cout << "Grille3d::compute_radial_grid : " << endl ;
    cout << "Unknown type of basis!" << endl ;
    abort() ;
    break ;
  }
}


}
