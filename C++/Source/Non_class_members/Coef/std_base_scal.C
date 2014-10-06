/*
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


char std_base_scal_C[] = "$Header$" ;

/*
 * Ensemble des routines de manipulation de base spectrales dans 
 * le cas scalaire.
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/06 15:16:04  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2013/06/05 14:54:46  j_novak
 * Removed the FINJAC sampling (now BASE_JAC02 in Mg3d).
 *
 * Revision 1.5  2007/12/11 15:28:17  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.4  2005/10/25 08:56:37  p_grandclement
 * addition of std_spectral_base in the case of odd functions near the origin
 *
 * Revision 1.3  2004/11/23 15:13:50  m_forot
 * Added the bases for the cases without any equatorial symmetry
 * (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.2  2002/10/16 14:36:57  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  1999/10/20  15:31:52  eric
 * La routine Valeur::std_base_scal() se trouve desormais dans le
 * fichier valeur.C.
 *
 * Revision 2.1  1999/03/01  15:00:43  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/22  15:30:33  hyc
 * *** empty log message ***
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
int std_base_scal_1z(int type_r, int type_t, int type_p) {
    
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
	// echantillonnage fin
	base_l = base_l | R_CHEB  ; // developpement en T_k(x) 
	break ;

      case RARE : 		 
	// echantillonnage rarefie
	base_l = base_l | R_CHEBPI_P ;  // developpement en 
	//  T_{2k}(x) pour l pair
	//  T_{2k+1}(x) pour l impair
	break ;
	
      case UNSURR : 		    
	// echantillonnage fin (1/r)
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "std_base_scal : le cas type_p, type_t, type_r = " 
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
	// echantillonnage fin
	base_l = base_l | R_CHEB  ; // developpement en T_k(x) 
	break ;

      case RARE : 		 
	// echantillonnage rarefie
	base_l = base_l | R_CHEBPIM_P ;  // developpement en 
	//  T_{2k}(x) pour m pair
	//  T_{2k+1}(x) pour m impair
	break ;
	
      case UNSURR : 		    
	// echantillonnage fin (1/r)
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "std_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p<< " " << type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
      
    default : 
      cout << 
	"std_base_scal : le cas type_p, type_t = " 
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
      case FIN : 	    // echantillonnage fin 
	base_l = base_l | R_CHEB  ; // developpement en T_k(x) 
	break ;

      case RARE :	    // echantillonnage rarefie		 
	base_l = base_l | R_CHEBPI_P ;  // developpement en 
	//  T_{2k}(x) pour l pair
	//  T_{2k+1}(x) pour l impair
	break ;
	
      case UNSURR :   // echantillonnage fin (1/r)	    
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "std_base_scal : le cas type_p, type_t, type_r = " 
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
      case FIN :	// echantillonnage fin	 
	base_l = base_l | R_CHEB  ; // developpement en T_k(x) 
	break ;
	
      case RARE :	    // echantillonnage rarefie	 
	base_l = base_l | R_CHEBP ;  // developpement en T_{2k}(x) 
	break ;
	
      case UNSURR :   // echantillonnage fin (1/r) 
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "std_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p<< " " <<type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
		    
    default : 
      cout << 
	"std_base_scal : le cas type_p, type_t = " 
	   << type_p<< " " <<type_t << endl ;
      cout << " n'est pas prevu ! " << endl ;
      abort () ;
    }	// fin des cas sur type_t 
    break ;	// fin du cas symetrie phi -> phi + pi
    
  default : 
    cout << 
      "std_base_scal : le cas type_p = " << type_p << endl ;
    cout << " n'est pas prevu ! " << endl ;
    abort () ;
  }	// Fin des cas en phi
  
  // On range le resultat
  return base_l ;
}

  		     //----------------------------------------//
		    // Le plus simple: cas une zone cas impair //
		    //----------------------------------------//

// Cree la base standart pour une zone
int std_base_scal_odd_1z(int type_r, int type_t, int type_p) {
    
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
	// echantillonnage fin
	base_l = base_l | R_CHEB  ; // developpement en T_k(x) 
	break ;	

      case RARE : 		 
	// echantillonnage rarefie
	base_l = base_l | R_CHEBPI_I ;  // developpement en 
	//  T_{2k}(x) pour l impair
	//  T_{2k+1}(x) pour l pair
	break ;
	
      case UNSURR : 		    
	// echantillonnage fin (1/r)
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "std_base_scal : le cas type_p, type_t, type_r = " 
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
	// echantillonnage fin
	base_l = base_l | R_CHEB  ; // developpement en T_k(x) 
	break ;
	
      case RARE : 		 
	// echantillonnage rarefie
	base_l = base_l | R_CHEBPIM_I ;  // developpement en 
	//  T_{2k}(x) pour m impair
	//  T_{2k+1}(x) pour m pair
	break ;
	
      case UNSURR : 		    
	// echantillonnage fin (1/r)
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "std_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p<< " " << type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
      
    default : 
      cout << 
	"std_base_scal : le cas type_p, type_t = " 
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
      case FIN : 	    // echantillonnage fin 
	base_l = base_l | R_CHEB  ; // developpement en T_k(x) 
	break ;

      case RARE :	    // echantillonnage rarefie		 
	base_l = base_l | R_CHEBPI_I ;  // developpement en 
	//  T_{2k}(x) pour l impair
	//  T_{2k+1}(x) pour l pair
	break ;
	
      case UNSURR :   // echantillonnage fin (1/r)	    
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "std_base_scal : le cas type_p, type_t, type_r = " 
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
      case FIN :	// echantillonnage fin	 
	base_l = base_l | R_CHEB  ; // developpement en T_k(x) 
	break ;
	
      case RARE :	    // echantillonnage rarefie	 
	base_l = base_l | R_CHEBI ;  // developpement en T_{2k+1}(x) 
	break ;
	
      case UNSURR :   // echantillonnage fin (1/r) 
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "std_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p<< " " <<type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
		    
    default : 
      cout << 
	"std_base_scal : le cas type_p, type_t = " 
	   << type_p<< " " <<type_t << endl ;
      cout << " n'est pas prevu ! " << endl ;
      abort () ;
    }	// fin des cas sur type_t 
    break ;	// fin du cas symetrie phi -> phi + pi
    
  default : 
    cout << 
      "std_base_scal : le cas type_p = " << type_p << endl ;
    cout << " n'est pas prevu ! " << endl ;
    abort () ;
  }	// Fin des cas en phi
  
  // On range le resultat
  return base_l ;
}

