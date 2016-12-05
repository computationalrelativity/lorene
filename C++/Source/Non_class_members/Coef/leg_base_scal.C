/*
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
 * le cas scalaire, avec base radiale Legendre .
 * 
 */

/*
 * $Id$
 * $Log $
 *
 */
//C headers
#include<cstdlib>

// Lorene
#include "headcpp.h"
#include "type_parite.h"


		    //------------------------------//
		    // Le plus simple: cas une zone //
		    //------------------------------//

// Cree la base standart pour une zone
namespace Lorene {
int leg_base_scal_1z(int type_r, int type_t, int type_p) {
    
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
	base_l = base_l | R_LEG  ; // developpement en P_k(x) 
	break ;

      // case RARE : 		 
      // 	// echantillonnage rarefie
      // 	base_l = base_l | R_LEGPI_P ;  // developpement en 
      // 	//  P_{2k}(x) pour l pair
      // 	//  P_{2k+1}(x) pour l impair
      // 	break ;
	
      case UNSURR : 		    
	// echantillonnage fin (1/r)
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "leg_base_scal : le cas type_p, type_t, type_r = " 
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
	base_l = base_l | R_LEG  ; // developpement en T_k(x) 
	break ;

      // case RARE : 		 
      // 	// echantillonnage rarefie
      // 	base_l = base_l | R_LEGPIM_P ;  // developpement en 
      // 	//  T_{2k}(x) pour m pair
      // 	//  T_{2k+1}(x) pour m impair
      // 	break ;
	
      case UNSURR : 		    
	// echantillonnage fin (1/r)
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "leg_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p<< " " << type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
      
    default : 
      cout << 
	"leg_base_scal : le cas type_p, type_t = " 
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
	base_l = base_l | R_LEG  ; // developpement en P_k(x) 
	break ;

      // case RARE :	    // echantillonnage rarefie		 
      // 	base_l = base_l | R_LEGPI_P ;  // developpement en 
      // 	//  P_{2k}(x) pour l pair
      // 	//  P_{2k+1}(x) pour l impair
      // 	break ;
	
      case UNSURR :   // echantillonnage fin (1/r)	    
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "leg_base_scal : le cas type_p, type_t, type_r = " 
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
	base_l = base_l | R_LEG  ; // developpement en P_k(x) 
	break ;
	
      case RARE :	    // echantillonnage rarefie	 
	base_l = base_l | R_LEGP ;  // developpement en P_{2k}(x) 
	break ;
	
      case UNSURR :   // echantillonnage fin (1/r) 
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "leg_base_scal : le cas type_p, type_t, type_r = " 
	     << type_p<< " " <<type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
		    
    default : 
      cout << 
	"leg_base_scal : le cas type_p, type_t = " 
	   << type_p<< " " <<type_t << endl ;
      cout << " n'est pas prevu ! " << endl ;
      abort () ;
    }	// fin des cas sur type_t 
    break ;	// fin du cas symetrie phi -> phi + pi
    
  default : 
    cout << 
      "leg_base_scal : le cas type_p = " << type_p << endl ;
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
int leg_base_scal_odd_1z(int type_r, int type_t, int type_p) {
    
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
	base_l = base_l | R_LEG  ; // developpement en P_k(x) 
	break ;	

      // case RARE : 		 
      // 	// echantillonnage rarefie
      // 	base_l = base_l | R_LEGPI_I ;  // developpement en 
      // 	//  P_{2k}(x) pour l impair
      // 	//  P_{2k+1}(x) pour l pair
      // 	break ;
	
      case UNSURR : 		    
	// echantillonnage fin (1/r)
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "leg_base_scal_odd : le cas type_p, type_t, type_r = " 
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
	base_l = base_l | R_LEG  ; // developpement en P_k(x) 
	break ;
	
      // case RARE : 		 
      // 	// echantillonnage rarefie
      // 	base_l = base_l | R_LEGPIM_I ;  // developpement en 
      // 	//  P_{2k}(x) pour m impair
      // 	//  P_{2k+1}(x) pour m pair
      // 	break ;
	
      case UNSURR : 		    
	// echantillonnage fin (1/r)
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "leg_base_scal_odd : le cas type_p, type_t, type_r = " 
	     << type_p<< " " << type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
      
    default : 
      cout << 
	"leg_base_scal_odd : le cas type_p, type_t = " 
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
	base_l = base_l | R_LEG  ; // developpement en P_k(x) 
	break ;

      // case RARE :	    // echantillonnage rarefie		 
      // 	base_l = base_l | R_LEGPI_I ;  // developpement en 
      // 	//  P_{2k}(x) pour l impair
      // 	//  P_{2k+1}(x) pour l pair
      // 	break ;
	
      case UNSURR :   // echantillonnage fin (1/r)	    
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "leg_base_scal_odd : le cas type_p, type_t, type_r = " 
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
	base_l = base_l | R_LEG  ; // developpement en T_k(x) 
	break ;
	
      case RARE :	    // echantillonnage rarefie	 
	base_l = base_l | R_LEGI ;  // developpement en T_{2k+1}(x) 
	break ;
	
      case UNSURR :   // echantillonnage fin (1/r) 
	base_l = base_l | R_CHEBU  ; // developpement en T_k(x) 
	break ;
	
      default : 
	cout << 
	  "leg_base_scal_odd : le cas type_p, type_t, type_r = " 
	     << type_p<< " " <<type_t<< " " <<type_r << endl ;
	cout << " n'est pas prevu ! " << endl ;
	abort () ;
      }
      break ;	// fin du cas type_t = SYM
		    
    default : 
      cout << 
	"leg_base_scal_odd : le cas type_p, type_t = " 
	   << type_p<< " " <<type_t << endl ;
      cout << " n'est pas prevu ! " << endl ;
      abort () ;
    }	// fin des cas sur type_t 
    break ;	// fin du cas symetrie phi -> phi + pi
    
  default : 
    cout << 
      "leg_base_scal_odd : le cas type_p = " << type_p << endl ;
    cout << " n'est pas prevu ! " << endl ;
    abort () ;
  }	// Fin des cas en phi
  
  // On range le resultat
  return base_l ;
}

}
