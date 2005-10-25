/*
 * Methods of class Mg3d to get the standard spectral bases for scalar and
 *  vector fields.
 */

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


char mg3d_std_base_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2005/10/25 08:56:37  p_grandclement
 * addition of std_spectral_base in the case of odd functions near the origin
 *
 * Revision 1.5  2005/02/16 15:09:16  m_forot
 * Add R_CHEBPI_I and R_CHEBPI_P cases
 *
 * Revision 1.4  2004/11/04 15:21:42  e_gourgoulhon
 * The case without any symmetry in theta is now treated.
 *
 * Revision 1.3  2003/12/19 16:21:45  j_novak
 * Shadow hunt
 *
 * Revision 1.2  2002/08/13 08:02:45  j_novak
 * Handling of spherical vector/tensor components added in the classes
 * Mg3d and Tenseur. Minor corrections for the class Metconf.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/09/27  15:07:40  eric
 * Correction dans le cas type_p = SYM.
 *
 * Revision 1.1  1999/10/12  14:54:43  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// headers C++

// headers C
#include <assert.h>

// headers Lorene
#include "grilles.h"
#include "base_val.h"
#include "type_parite.h"

int std_base_scal_1z(int type_r, int type_t, int type_p) ; 
int std_base_scal_odd_1z(int type_r, int type_t, int type_p) ; 

		    //-----------------------------//
		    //	Bases for a scalar field   //
		    //-----------------------------//
		    

Base_val Mg3d::std_base_scal() const {
          
    Base_val base(nzone) ;  
     
    for (int l=0; l<nzone; l++) {
	base.b[l] = std_base_scal_1z(type_r[l], type_t, type_p) ;
    }
    
    return base ; 
     
}    

Base_val Mg3d::std_base_scal_odd() const {
          
    Base_val base(nzone) ;  
     
    for (int l=0; l<nzone; l++) {
	base.b[l] = std_base_scal_odd_1z(type_r[l], type_t, type_p) ;
    }
    
    return base ; 
     
}    

		    //---------------------------------------//
		    //	Bases for the Cartesian components   //
		    //    of a vector field		     //
		    //---------------------------------------//
		    

/*
 * Calcul les bases spectrales associees aux composantes cartesiennes d'un vecteur
 * antisymetrique en z (pour la composante z)
 * 
 * (*THIS) est la grille du calcul
 * SORTIE : un tableau sur les 3 compsantes (x=1, y=2, z=3) contenant les bases
 * de decomposition
 * 
 */

Base_val** Mg3d::std_base_vect_cart() const {
     
     // nbre de zones :
     int nz = get_nzone() ;
     
     // Tableau contenant le resultat...
     Base_val** bases = new (Base_val* [3]) ;
     for (int i=0 ; i<3 ; i++)
	bases[i] = new Base_val(nz) ;
     	        
    // Boucle sur les differentes zones :
    for (int l=0; l<nzone; l++) {

    // Type d'echantillonnage de la zone l :
	int type_r0 = get_type_r(l) ;

    // Bases de developpement en (r,theta,phi) a determiner pour les composantes
    // (1,2,3) du vecteur : 
	
	int base1,  base2, base3 ;	
	switch ( type_p ) {
	
	    case NONSYM : 	
//---------------------------------------------------------
// Cas sans symetrie sur phi : phi dans [0, 2 pi[
//---------------------------------------------------------

// Base en phi:
//-------------
	    base1 = P_COSSIN ;	    
	    base2 = P_COSSIN ;	    
	    base3 = P_COSSIN ;	    
    
    
// Base en theta:
//---------------
	    switch ( type_t ) {
		case SYM :  	
// symetrie theta -> pi - theta :  theta dans [0, pi/2]	    
//------------------------------------------------------
		base1 = base1 | T_COSSIN_CP ; 
		base2 = base2 | T_COSSIN_CP ; 
		base3 = base3 | T_COSSIN_CI ; 

// Base en r :
//------------
		switch ( type_r0 ) {
		    
		    case FIN : 			 
// echantillonnage fin
		    
		    base1 = base1 | R_CHEB  ;  
		    base2 = base2 | R_CHEB  ;  
		    base3 = base3 | R_CHEB  ;  
		    break ;

		    case RARE : 		 
// echantillonnage rarefie

		    base1 = base1 | R_CHEBPIM_P ;  
		    base2 = base2 | R_CHEBPIM_P ;  
		    base3 = base3 | R_CHEBPIM_I ;  
		    
		    break ;

		    case UNSURR : 		    
// echantillonnage fin (1/r)

		    base1 = base1 | R_CHEBU  ;  
		    base2 = base2 | R_CHEBU  ;  
		    base3 = base3 | R_CHEBU  ;  
		    break ;


		    default : 
			cout << 
	    "Mg3d::std_base_vect_cart : le cas type_p, type_t, type_r = " 
			  << type_p<< " " << type_t<< " " <<type_r0 << endl ;
			cout << 
			  " dans la zone l = " << l << " n'est pas prevu ! " 
			  << endl ;
			  abort () ;
		}

		break ;	// fin du cas type_t = SYM
		    

		case NONSYM :  	
// pas de symetrie en theta  :  theta dans [0, pi]	    
//------------------------------------------------
		base1 = base1 | T_COSSIN_C ; 
		base2 = base2 | T_COSSIN_C ; 
		base3 = base3 | T_COSSIN_C ; 

// Base en r :
//------------
		switch ( type_r0 ) {
		    
		    case FIN : 			 
// echantillonnage fin
		    
		    base1 = base1 | R_CHEB  ;  
		    base2 = base2 | R_CHEB  ;  
		    base3 = base3 | R_CHEB  ;  
		    break ;

		    case RARE : 		 
// echantillonnage rarefie

		    base1 = base1 | R_CHEBPI_P ;  
		    base2 = base2 | R_CHEBPI_P ;  
		    base3 = base3 | R_CHEBPI_P ;  
		    
		    break ;

		    case UNSURR : 		    
// echantillonnage fin (1/r)

		    base1 = base1 | R_CHEBU  ;  
		    base2 = base2 | R_CHEBU  ;  
		    base3 = base3 | R_CHEBU  ;  
		    break ;


		    default : 
			cout << 
	    "Mg3d::std_base_vect_cart : le cas type_p, type_t, type_r = " 
			  << type_p<< " " << type_t<< " " <<type_r0 << endl ;
			cout << 
			  " dans la zone l = " << l << " n'est pas prevu ! " 
			  << endl ;
			  abort () ;
		}

		break ;	// fin du cas type_t = NONSYM
		    


		    
		default : 
		    cout << 
		"Mg3d::std_base_vect_cart : le cas type_p, type_t = " 
			  << type_p<< " " <<type_t << endl ;
		    cout << 
			  " dans la zone l = " << l << " n'est pas prevu ! " 
			  << endl ;
		    abort () ;

	    }	// fin des cas sur type_t 



	    break ;	// fin du cas sans symetrie pour phi 


	    case SYM : 	
//---------------------------------------------------------
// Cas symetrie phi -> phi + pi :  phi in [0, pi]
//---------------------------------------------------------

// Base en phi:
//-------------
	    base1 = P_COSSIN_I ;	   
	    base2 = P_COSSIN_I ;	   
	    base3 = P_COSSIN_P ;	   
    
    
// Base en theta:
//---------------
	    switch ( type_t ) {

		case SYM :  	
// symetrie theta -> pi - theta :  theta dans [0, pi/2]	    
//------------------------------------------------------
		base1 = base1 | T_SIN_I ;	    
		base2 = base2 | T_SIN_I ;	    
		base3 = base3 | T_COS_I;	    

// Base en r :
//------------
		switch ( type_r0 ) {
		    
		    case FIN : 			 
// echantillonnage fin
		    
		    base1 = base1 | R_CHEB  ;  
		    base2 = base2 | R_CHEB  ;  
		    base3 = base3 | R_CHEB  ;  
		    break ;

		    case RARE : 		 
// echantillonnage rarefie

		    base1 = base1 | R_CHEBI ;  
		    base2 = base2 | R_CHEBI ;  
		    base3 = base3 | R_CHEBI ;  
		    break ;

		    case UNSURR : 		    
// echantillonnage fin (1/r)

		    base1 = base1 | R_CHEBU  ;  
		    base2 = base2 | R_CHEBU  ;  
		    base3 = base3 | R_CHEBU  ;  
		    break ;


		    default : 
			cout << 
	    "Mg3d::std_base_vect_cart : le cas type_p, type_t, type_r = " 
			  << type_p<< " " <<type_t<< " " <<type_r0 << endl ;
			cout << 
			  " dans la zone l = " << l << " n'est pas prevu ! " 
			  << endl ;
			  abort () ;
		}

		break ;	// fin du cas type_t = SYM
		    
		    
		default : 
		    cout << 
	    "Mg3d::std_base_vect_cart : le cas type_p, type_t = " 
			  << type_p<< " " <<type_t << endl ;
		    cout << 
			  " dans la zone l = " << l << " n'est pas prevu ! " 
			  << endl ;
		    abort () ;

	    }	// fin des cas sur type_t 



	break ;	// fin du cas symetrie phi -> phi + pi


	default : 
	    cout << 
	    "Mg3d::std_base_vect_cart : le cas type_p = " << type_p << endl ;
	    cout << " dans la zone l = " << l << " n'est pas prevu ! " 
		 << endl ;
	    abort () ;
	    
	    
	}	// Fin des cas en phi
	
	bases[0]->b[l] = base1 ;
	bases[1]->b[l] = base2 ;
	bases[2]->b[l] = base3 ;
    }	//fin de la boucle sur les zones.
  
    return bases ;
}
		    //---------------------------------------//
		    //	Bases for the spherical components   //
		    //    of a vector field		     //
		    //---------------------------------------//
		    

/*
 * Calcul les bases spectrales associees aux composantes spheriques d'un 
 * vecteur antisymetrique en z (pour la composante z)
 * 
 * (*THIS) est la grille du calcul
 * SORTIE : un tableau sur les 3 compsantes (r=1, theta=2, phi=3) contenant 
 * les bases de decomposition
 * 
 */

Base_val** Mg3d::std_base_vect_spher() const {
     
  // nbre de zones :
  int nz = get_nzone() ;
  
  // Tableau contenant le resultat...
  Base_val** bases = new (Base_val* [3]) ;
  for (int i=0 ; i<3 ; i++)
    bases[i] = new Base_val(nz) ;
  
  // Boucle sur les differentes zones :
  for (int l=0; l<nzone; l++) {
    
    // Type d'echantillonnage de la zone l :
    int type_r0 = get_type_r(l) ;
    
    // Bases de developpement en (r,theta,phi) a determiner pour les 
    // composantes (1,2,3) du vecteur : 
	
    int base1,  base2, base3 ;	
    switch ( type_p ) {
	
    case NONSYM : 	
//---------------------------------------------------------
// Cas sans symetrie sur phi : phi dans [0, 2 pi[
//---------------------------------------------------------

      // Base en phi:
      //-------------
      base1 = P_COSSIN ;	    
      base2 = P_COSSIN ;	    
      base3 = P_COSSIN ;	    
    
    
      // Base en theta:
      //---------------
      switch ( type_t ) {
      case SYM :  	
// symetrie theta -> pi - theta :  theta dans [0, pi/2]	    
//------------------------------------------------------
	base1 = base1 | T_COSSIN_CP ; 
	base2 = base2 | T_COSSIN_SP ; 
	base3 = base3 | T_COSSIN_SI ; 

	// Base en r :
	//------------
	switch ( type_r0 ) {
		    
	case FIN : 			 
// echantillonnage fin
		    
	  base1 = base1 | R_CHEB  ;  
	  base2 = base2 | R_CHEB  ;  
	  base3 = base3 | R_CHEB  ;  
	  break ;

	case RARE : 		 
// echantillonnage rarefie

	  base1 = base1 | R_CHEBPIM_I ;  
	  base2 = base2 | R_CHEBPIM_I ;  
	  base3 = base3 | R_CHEBPIM_I ;  
		    
	  break ;

	case UNSURR : 		    
// echantillonnage fin (1/r)

	  base1 = base1 | R_CHEBU  ;  
	  base2 = base2 | R_CHEBU  ;  
	  base3 = base3 | R_CHEBU  ;  
	  break ;


	default : 
	  cout << 
	    "Mg3d::std_base_vect_sphere : le cas type_p, type_t, type_r = " 
	       << type_p<< " " << type_t<< " " <<type_r0 << endl ;
	  cout << " dans la zone l = " << l << " n'est pas prevu ! " 
	       << endl ;
	  abort () ;
	}

	break ;	// fin du cas type_t = SYM
		    
      case NONSYM :  	
// pas de symetrie en theta  :  theta dans [0, pi]	    
//------------------------------------------------

	base1 = base1 | T_COSSIN_C ; 
	base2 = base2 | T_COSSIN_S ; 
	base3 = base3 | T_COSSIN_S ; 

	// Base en r :
	//------------
	switch ( type_r0 ) {
		    
	case FIN : 			 
// echantillonnage fin
		    
	  base1 = base1 | R_CHEB  ;  
	  base2 = base2 | R_CHEB  ;  
	  base3 = base3 | R_CHEB  ;  
	  break ;

	case RARE : 		 
// echantillonnage rarefie

	  base1 = base1 | R_CHEBPI_I ;  
	  base2 = base2 | R_CHEBPI_I ;  
	  base3 = base3 | R_CHEBPI_P ;  
		    
	  break ;

	case UNSURR : 		    
// echantillonnage fin (1/r)

	  base1 = base1 | R_CHEBU  ;  
	  base2 = base2 | R_CHEBU  ;  
	  base3 = base3 | R_CHEBU  ;  
	  break ;


	default : 
	  cout << 
	    "Mg3d::std_base_vect_sphere : le cas type_p, type_t, type_r = " 
	       << type_p<< " " << type_t<< " " <<type_r0 << endl ;
	  cout << " dans la zone l = " << l << " n'est pas prevu ! " 
	       << endl ;
	  abort () ;
	}

	break ;	// fin du cas type_t = SYM
		    
		    
      default : 
	cout << "Mg3d::std_base_vect_spher : le cas type_p, type_t = " 
	     << type_p<< " " <<type_t << endl ;
	cout << " dans la zone l = " << l << " n'est pas prevu ! " 
	     << endl ;
	abort () ;

      }	// fin des cas sur type_t 

      break ;	// fin du cas sans symetrie pour phi 


    case SYM : 	
//---------------------------------------------------------
// Cas symetrie phi -> phi + pi :  phi in [0, pi]
//---------------------------------------------------------

      // Base en phi:
      //-------------
      base1 = P_COSSIN_P ;	   
      base2 = P_COSSIN_P ;	   
      base3 = P_COSSIN_P ;	   
    
    
      // Base en theta:
      //---------------
      switch ( type_t ) {

      case SYM :  	
// symetrie theta -> pi - theta :  theta dans [0, pi/2]	    
//------------------------------------------------------
	base1 = base1 | T_COS_P ;	    
	base2 = base2 | T_SIN_P ;	    
	base3 = base3 | T_SIN_I;	    

	// Base en r :
	//------------
	switch ( type_r0 ) {
		    
	case FIN : 			 
// echantillonnage fin
		    
	  base1 = base1 | R_CHEB  ;  
	  base2 = base2 | R_CHEB  ;  
	  base3 = base3 | R_CHEB  ;  
	  break ;

	case RARE : 		 
// echantillonnage rarefie

	  base1 = base1 | R_CHEBI ;  
	  base2 = base2 | R_CHEBI ;  
	  base3 = base3 | R_CHEBI ;  
	  break ;

	case UNSURR : 		    
// echantillonnage fin (1/r)

	  base1 = base1 | R_CHEBU  ;  
	  base2 = base2 | R_CHEBU  ;  
	  base3 = base3 | R_CHEBU  ;  
	  break ;


	default : 
	  cout << 
	    "Mg3d::std_base_vect_spher : le cas type_p, type_t, type_r = " 
	       << type_p<< " " <<type_t<< " " <<type_r0 << endl ;
	  cout << " dans la zone l = " << l << " n'est pas prevu ! " 
	       << endl ;
	  abort () ;
	}

	break ;	// fin du cas type_t = SYM
		    
      default : 
	cout << "Mg3d::std_base_vect_spher : le cas type_p, type_t = " 
	     << type_p<< " " <<type_t << endl ;
	cout << " dans la zone l = " << l << " n'est pas prevu ! " 
	     << endl ;
	abort () ;
	
      }	// fin des cas sur type_t 

      break ;	// fin du cas symetrie phi -> phi + pi

    default : 
      cout << 
	"Mg3d::std_base_vect_spher : le cas type_p = " << type_p << endl ;
      cout << " dans la zone l = " << l << " n'est pas prevu ! " 
	   << endl ;
      abort () ;
	    
	    
    }	// Fin des cas en phi
	
    bases[0]->b[l] = base1 ;
    bases[1]->b[l] = base2 ;
    bases[2]->b[l] = base3 ;
  }	//fin de la boucle sur les zones.
  
  return bases ;
}

