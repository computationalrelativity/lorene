/*
 *   Copyright (c) 2007 Jean-Louis Cornou
 *   Copyright (c) 2013 Jerome Novak
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
 * Ensemble des routines de manipulation de base spectrales dans 
 * le cas scalaire et la base Jacobi(0,2).
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.3  2016/12/05 16:18:02  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:53:12  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2013/06/05 15:10:43  j_novak
 * Suppression of FINJAC sampling in r. This Jacobi(0,2) base is now
 * available by setting colloc_r to BASE_JAC02 in the Mg3d constructor.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>

// Lorene
#include "headcpp.h"
#include "type_parite.h"


		    //------------------------------//
		    // Le plus simple: cas une zone //
		    //------------------------------//

// Cree la base standart pour une zone
namespace Lorene {
int jac02_base_scal_1z(int type_r, int type_t, int type_p) {
    
  // Base d'echantillonnage en (r,theta,phi) a determiner :
  int base_l  = 0 ;
  
  // proccess phi
  switch ( type_p ) {
  case NONSYM : 	
    // Cas sans symetrie sur phi : phi dans [0, 2 pi[
    base_l = P_COSSIN ;	    // developpement en cos,sin(m*phi)
    // Base en theta:
    switch ( type_t ) {
    case NONSYM :  	
      // pas de symetrie en theta : theta dans [0,pi]    
      base_l = base_l | T_COSSIN_C ; // developpement en 
      //	cos(l*theta) pour m pair
      //	sin(l*theta) pour m impair
      
      // Base en r :
      switch ( type_r ) {
      case FIN : 			 
	base_l = base_l | R_JACO02  ; // developpement en J_k(x) 
	break ;

      default : 
	cout << 
	  "jac02_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p << " " << type_t << " " <<  type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = NONSYM
      
    case SYM :  	// en theta
      // symetrie theta -> pi - theta :  theta dans [0, pi/2]	    
      base_l = base_l | T_COSSIN_CP ; // developpement en 
      // cos(2*l*theta) pour m pair
      // sin((2*l+1)*theta) pour m impair

      // Base en r :
      switch ( type_r ) {
      case FIN :   			 
	base_l = base_l | R_JACO02  ; // developpement en J_k(x) 
	break ;

      default : 
	cout << 
	  "jac02_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p<< " " << type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
      
    default : 
      cout << 
	"jac02_base_scal : le cas type_p, type_t = " 
	   << type_p<< " " <<type_t << endl ;
      cout << " n'est pas prevu ! " << endl ;
      abort () ;
    }	// fin des cas sur type_t 
    break ;	// fin du cas sans symetrie pour phi 
    
    
  case SYM : 	// en phi
    // Cas symetrie phi -> phi + pi :  phi in [0, pi]
    base_l = P_COSSIN_P ;	    // developpement en cos,sin(2*m*phi)
    // Base en theta:
    switch ( type_t ) {
    case NONSYM :  	
      // pas de symetrie en theta : theta dans [0,pi]    
      base_l = base_l | T_COS ;   // developpement en cos(l*theta) seulement
      //  (puisque m est toujours pair) 
      // Base en r :
      switch ( type_r ) {
      case FIN :  
	base_l = base_l | R_JACO02  ; // developpement en J_k(x) 
	break ;
	
      default : 
	cout << 
	  "jac02_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p<< " " <<type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = NONSYM
      
    case SYM :  // symetrie theta -> pi - theta :  theta dans [0, pi/2]
      base_l = base_l | T_COS_P ;	// developpement en cos(2*l*theta)
      //  (puisque m est toujours pair)
      // Base en r :
      switch ( type_r ) {
      case FIN : 
	base_l = base_l | R_JACO02  ; // developpement en J_k(x) 
	break ;
	
      default : 
	cout << 
	  "jac02_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p<< " " <<type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
		    
    default : 
      cout << 
	"jac02_base_scal : le cas type_p, type_t = " 
	   << type_p<< " " <<type_t << endl ;
      cout << " n'est pas prevu ! " << endl ;
      abort () ;
    }	// fin des cas sur type_t 
    break ;	// fin du cas symetrie phi -> phi + pi
    
  default : 
    cout << 
      "jac02_base_scal : le cas type_p = " << type_p << endl ;
    cout << " n'est pas prevu ! " << endl ;
    abort () ;
  }	// Fin des cas en phi
  
  // On range le resultat
  return base_l ;
}
}
