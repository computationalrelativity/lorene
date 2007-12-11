/*
 *  Methods of class Base_val
 *
 *   (see file base_val.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char base_val_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2007/12/11 15:28:09  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.8  2004/12/17 13:35:01  m_forot
 * Add the case T_LEG
 *
 * Revision 1.7  2004/11/04 15:19:02  e_gourgoulhon
 * operator<< : added names R_CHEBPI_P, R_CHEBPI_I, T_COSSIN_C, T_COSSIN_S.
 *
 * Revision 1.6  2004/08/24 09:14:41  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.5  2003/09/16 08:54:09  j_novak
 * Addition of the T_LEG_II base (odd in theta, only for odd m) and the
 * transformation functions to and from the T_SIN_P base.
 *
 * Revision 1.4  2002/11/13 15:05:59  j_novak
 * Affichage de la base T_COS
 *
 * Revision 1.3  2002/10/16 14:36:30  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2001/12/04 21:27:52  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.8  2000/09/28  10:20:19  eric
 * Affichage: nouvelles bases T_LEG_IP et T_LEG_PI.
 *
 * Revision 2.7  2000/09/08  10:06:45  eric
 * Ajout des methodes set_base_r, etc...
 *
 * Revision 2.6  1999/12/28  12:57:26  eric
 * Reorganisation des headers.
 *
 * Revision 2.5  1999/10/12  10:02:51  eric
 * Implementation de sauve().
 * Modif de << : affichage du nom des bases et non plus de leur numero.
 *
 * Revision 2.4  1999/10/01  15:56:21  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.3  1999/09/13  14:57:06  phil
 * ajout de l'operateur ==
 *
 * Revision 2.2  1999/03/02  15:22:23  eric
 * Affichage des bases.
 *
 * Revision 2.1  1999/03/01  14:54:01  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/22  15:19:20  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdio.h>
#include <assert.h>

// Headers Lorene
#include "headcpp.h"
#include "indent.h"
#include "type_parite.h"
#include "base_val.h"
#include "utilitaires.h"


			//---------------//
			// Constructeurs //
			//---------------//

// Constructeur
Base_val::Base_val(int n) : nzone(n) {
    b = new int[nzone] ;
    for (int i=0 ; i<nzone ; i++) {    // Boucle sur les zones
	b[i] = NONDEF ;
    }
}

// Copie
Base_val::Base_val(const Base_val & bi) : nzone(bi.nzone) {
    b = new int[nzone] ;
    for (int i=0 ; i<nzone ; i++) {    // Boucle sur les zones
	b[i] = bi.b[i] ;
    }
}
	
// From file
Base_val::Base_val(FILE* fd) {
  fread_be(&nzone, sizeof(int), 1, fd) ;		// nzone
  b = new int[nzone] ; 
  fread_be(b, sizeof(int), nzone, fd) ;		// b[]
}

			//--------------//
			// Destructeurs //
			//--------------//

// Destructeur
Base_val::~Base_val() {
    delete [] b ;
}

			//-------------//
			// Affectation //
			//-------------//

void Base_val::set_base_r(int l, int base_r) {
    
    assert( (l>=0) && (l<nzone) ) ; 

    int base_t = b[l] & MSQ_T ;
    int base_p = b[l] & MSQ_P ;
    b[l] = base_p | base_t | base_r ;
    
}

void Base_val::set_base_t(int base_t) {
    
    for (int l=0; l<nzone; l++) {
	int base_r = b[l] & MSQ_R ;
	int base_p = b[l] & MSQ_P ;
	b[l] = base_p | base_t | base_r ;
    }
}

void Base_val::set_base_p(int base_p) {
    
    for (int l=0; l<nzone; l++) {
	int base_r = b[l] & MSQ_R ;
	int base_t = b[l] & MSQ_T ;
	b[l] = base_p | base_t | base_r ;
    }
}


// From Base_val
void Base_val::operator=(const Base_val & bi) {
    // Protection
    assert(nzone == bi.nzone) ;
    for (int i=0 ; i<nzone ; i++) {    // Boucle sur les zones
	b[i] = bi.b[i] ;
    }
}
    
			//------------//
			// Sauvegarde //
			//------------//

// Save in a file
void Base_val::sauve(FILE* fd) const {
    fwrite_be(&nzone, sizeof(int), 1, fd) ;	    // nzone
    fwrite_be(b, sizeof(int), nzone, fd) ;	            // b[]
}
    
			//------------//
			// Impression //
			//------------//

// Operateurs <<
ostream& operator<<(ostream& o, const Base_val & bi) {
  
  static bool premier_appel = true ; 
  static char* nom_r[MAX_BASE] ;
  static char* nom_t[MAX_BASE] ;
  static char* nom_p[MAX_BASE] ;

  if (premier_appel) {   // First call initializations

    premier_appel = false ;

    for (int i=0; i<MAX_BASE; i++) {
      nom_r[i] = "UNDEFINED" ;    
      nom_t[i] = "UNDEFINED" ;    
      nom_p[i] = "UNDEFINED" ;    
    }

    nom_r[R_CHEB >> TRA_R] =      "R_CHEB     " ; 
    nom_r[R_CHEBU >> TRA_R] =     "R_CHEBU    " ; 
    nom_r[R_CHEBP >> TRA_R] =     "R_CHEBP    " ; 
    nom_r[R_CHEBI >> TRA_R] =     "R_CHEBI    " ; 
    nom_r[R_CHEBPIM_P >> TRA_R] = "R_CHEBPIM_P" ; 
    nom_r[R_CHEBPIM_I >> TRA_R] = "R_CHEBPIM_I" ; 
    nom_r[R_CHEBPI_P >> TRA_R] =  "R_CHEBPI_P " ; 
    nom_r[R_CHEBPI_I >> TRA_R] =  "R_CHEBPI_I " ; 
    nom_r[R_JACO02 >> TRA_R] =    "R_JACO02   " ;
  
    nom_t[T_COS >> TRA_T] =       "T_COS      " ; 
    nom_t[T_COS_P >> TRA_T] =     "T_COS_P    " ; 
    nom_t[T_COS_I >> TRA_T] =     "T_COS_I    " ; 
    nom_t[T_SIN_P >> TRA_T] =     "T_SIN_P    " ; 
    nom_t[T_SIN_I >> TRA_T] =     "T_SIN_I    " ; 
    nom_t[T_COSSIN_CP >> TRA_T] = "T_COSSIN_CP" ;
    nom_t[T_COSSIN_SI >> TRA_T] = "T_COSSIN_SI" ;
    nom_t[T_COSSIN_SP >> TRA_T] = "T_COSSIN_SP" ;
    nom_t[T_COSSIN_CI >> TRA_T] = "T_COSSIN_CI" ;
    nom_t[T_COSSIN_C >> TRA_T] =  "T_COSSIN_C " ;
    nom_t[T_COSSIN_S >> TRA_T] =  "T_COSSIN_S " ;
    nom_t[T_LEG >> TRA_T] =       "T_LEG      " ;
    nom_t[T_LEG_P >> TRA_T] =     "T_LEG_P    " ;
    nom_t[T_LEG_PP >> TRA_T] =    "T_LEG_PP   " ;
    nom_t[T_LEG_I >> TRA_T] =     "T_LEG_I    " ;
    nom_t[T_LEG_IP >> TRA_T] =    "T_LEG_IP   " ;
    nom_t[T_LEG_PI >> TRA_T] =    "T_LEG_PI   " ;
    nom_t[T_LEG_II >> TRA_T] =    "T_LEG_II   " ;
    nom_t[T_CL_COS_P >> TRA_T] =  "T_CL_COS_P " ;
    nom_t[T_CL_SIN_P >> TRA_T] =  "T_CL_SIN_P " ;
    nom_t[T_CL_COS_I >> TRA_T] =  "T_CL_COS_I " ;
    nom_t[T_CL_SIN_I >> TRA_T] =  "T_CL_SIN_I " ;

    nom_p[P_COSSIN >> TRA_P] =    "P_COSSIN   " ;
    nom_p[P_COSSIN_P >> TRA_P] =  "P_COSSIN_P " ;
    nom_p[P_COSSIN_I >> TRA_P] =  "P_COSSIN_I " ;


  } // End of first call operations
  

    // Intro - Nombre de zones
    int nzone = bi.nzone ;
    o << "Bases of spectral expansions: " << incindent ;
    for (int l=0 ; l<nzone ; l++) {
	int base_r = (bi.b[l] & MSQ_R) >> TRA_R ;
	int base_t = (bi.b[l] & MSQ_T) >> TRA_T  ;
	int base_p = (bi.b[l] & MSQ_P) >> TRA_P ;
	o << iendl ;
	o << "Domain #" << l << " : r: " << nom_r[base_r] 
	    << ",  theta: " << nom_t[base_t] 
	    << ",  phi: " << nom_p[base_p] ;
    }
    o << decindent << decindent << iendl ;

    //Termine
    return o ;
}
    
			//----------------------//
			// Manipulation de base //
			//----------------------//

void Base_val::set_base_nondef() {
    for (int i=0 ; i<nzone ; i++) {
	b[i] = NONDEF ;
    }
}
			//----------------------//
			// operateur logique    //
			//----------------------//

bool Base_val::operator== (const Base_val& c2) const {

    return (*b == *c2.b) ;
}
