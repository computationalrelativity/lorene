/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2001 Eric Gourgoulhon
 *   Copyright (c) 2001 Keisuke Taniguchi
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


char base_val_mult_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2014/10/06 15:12:56  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2013/01/11 08:20:11  j_novak
 * New radial spectral bases with Legendre polynomials (R_LEG, R_LEGP, R_LEGI).
 *
 * Revision 1.8  2009/10/23 12:55:16  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.7  2009/10/08 16:20:13  j_novak
 * Addition of new bases T_COS and T_SIN.
 *
 * Revision 1.6  2008/08/27 08:46:30  jl_cornou
 * Added R_JACO02 base (Jacobi(0,2) polynomials)
 *
 * Revision 1.5  2004/11/23 15:08:00  m_forot
 * Added the bases for the cases without any equatorial symmetry
 * (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.4  2002/10/16 14:36:30  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/08/02 15:07:41  j_novak
 * Member function determinant has been added to the class Metrique.
 * A better handling of spectral bases is now implemented for the class Tenseur.
 *
 * Revision 1.2  2002/02/07 14:55:07  e_gourgoulhon
 * Add more cases in theta and phi
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2001/08/29  09:31:00  keisuke
 * Addition of the cases T_COSSIN_SP * T_COSSIN_SP,
 *  T_COSSIN_SI * T_COSSIN_SI, etc.
 *
 * Revision 2.2  2001/08/27  14:59:27  keisuke
 * Ajout du cas T_COSSIN_CP * T_COSSIN_SI
 *
 * Revision 2.1  2001/08/27  13:40:18  eric
 * Ajout du cas T_COSSIN_CP * T_COSSIN_SP
 *
 * Revision 2.0  1999/10/26  14:42:47  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Fichier includes
#include <cstdlib>
#include <cstdio>
#include <cassert>

#include "headcpp.h"
#include "type_parite.h"
#include "base_val.h"

/*
 * Routine calculant le produit de deux bases spectrales en utilisant en fait 
 * le produit des symetries par rapport au plan z=0
 * 
 * Si le resultat n'est pas defini le resultat est dans etat == ETATNONDEF
 * 
 */

Base_val operator* (const Base_val& b1, const Base_val& b2) {
    
    assert (b1.nzone == b2.nzone) ;
    
    Base_val res(b1.nzone) ;
    
    int base, indic_r, indic_t, indic_p ;
    int b1_r, b2_r, b1_t, b2_t, b1_p, b2_p ; // Confort ;
    
    int indic_total = 1 ;
    
    //Boucle sur les zones :
    for (int l=0 ; l<b1.nzone ; l++) {
	
	indic_r = -1 ;
	indic_t = -1 ;
	indic_p = -1 ;
	
	b1_r = b1.b[l] & MSQ_R ;
	b1_t = b1.b[l] & MSQ_T ;
	b1_p = b1.b[l] & MSQ_P ;
	b2_r = b2.b[l] & MSQ_R ;
	b2_t = b2.b[l] & MSQ_T ;
	b2_p = b2.b[l] & MSQ_P ;
	
	base = 0 ;
	
	switch (b1_p) {
	    case P_COSSIN :
		switch (b2_p) {
		    case P_COSSIN :
			base = P_COSSIN ;
			indic_p = 1 ;
			break ;
		    default :
			break ;
		}
		break ;
		
	    case P_COSSIN_P :
		switch (b2_p) {
		    case P_COSSIN_P :
			base = P_COSSIN_P ;
			indic_p = 1 ;
			break ;

                    case P_COSSIN_I :
			base = P_COSSIN_I ;
			indic_p = 1 ;
			break ;

                    default :
			break ;
		}
		break ;

	    case P_COSSIN_I :
		switch (b2_p) {
		    case P_COSSIN_P :
			base = P_COSSIN_I ;
			indic_p = 1 ;
			break ;

                    case P_COSSIN_I :
			base = P_COSSIN_P ;
			indic_p = 1 ;
			break ;

                    default :
			break ;
		}
		break ;

	    default :
		break ;
	}

	switch (b1_t) {
	    
	    case T_COSSIN_CP :
		switch (b2_t) {
		    case T_COSSIN_CP :
			base = base | T_COSSIN_CP ;
			indic_t = 1 ;
			break ;
			
		    case T_COSSIN_CI :
			base = base | T_COSSIN_CI ;
			indic_t = 1 ;
			break ;
			
		    case T_COSSIN_SP :
			base = base | T_COSSIN_SP ;
			indic_t = 1 ;
			break ;

		    case T_COSSIN_SI :
			base = base | T_COSSIN_SI ;
			indic_t = 1 ;
			break ;
			
		    default :
			break ;
		}
		break ;
		
	    case T_COSSIN_CI : 
		switch (b2_t) {
		    case T_COSSIN_CP :
			base = base | T_COSSIN_CI ;
			indic_t = 1 ;
			break ;
			
		    case T_COSSIN_CI :
			base = base | T_COSSIN_CP ;
			indic_t = 1 ;
			break ;
			
		    case T_COSSIN_SP :
			base = base | T_COSSIN_SI ;
			indic_t = 1 ;
			break ;

		    case T_COSSIN_SI :
			base = base | T_COSSIN_SP ;
			indic_t = 1 ;
			break ;

		    default :
			break ; 
		}
		break ;
		
	    case T_COSSIN_SP :
		switch (b2_t) {
		    case T_COSSIN_CP :
			base = base | T_COSSIN_SP ;
			indic_t = 1 ;
			break ;

		    case T_COSSIN_CI :
			base = base | T_COSSIN_SI ;
			indic_t = 1 ;
			break ;

		    case T_COSSIN_SP :
			base = base | T_COSSIN_CP ;
			indic_t = 1 ;
			break ;
			
		    case T_COSSIN_SI :
			base = base | T_COSSIN_CI ;
			indic_t = 1 ;
			break ;

		    default :
			break ;
		}
		break ;

	    case T_COSSIN_SI :
		switch (b2_t) {
		    case T_COSSIN_CP :
			base = base | T_COSSIN_SI ;
			indic_t = 1 ;
			break ;

		    case T_COSSIN_CI :
			base = base | T_COSSIN_SP ;
			indic_t = 1 ;
			break ;

		    case T_COSSIN_SP :
			base = base | T_COSSIN_CI ;
			indic_t = 1 ;
			break ;

		    case T_COSSIN_SI :
			base = base | T_COSSIN_CP ;
			indic_t = 1 ;
			break ;
			
		    default :
			break ;
		}
		break ;
		
	    case T_COS_P :
		switch (b2_t) {
		    case T_COS_P :
			base = base | T_COS_P ;
			indic_t = 1 ;
			break ;

		    case T_COS_I :
			base = base | T_COS_I ;
			indic_t = 1 ;
			break ;

		    case T_SIN_I :
			base = base | T_SIN_I ;
			indic_t = 1 ;
			break ;

		    case T_SIN_P :
			base = base | T_SIN_P ;
			indic_t = 1 ;
			break ;

                    default :
			break ;
		}
		break ;

	    case T_COS_I :
		switch (b2_t) {
		    case T_COS_P :
			base = base | T_COS_I ;
			indic_t = 1 ;
			break ;
		    
		    case T_COS_I :
			base = base | T_COS_P ;
			indic_t = 1 ;
			break ;

		    case T_SIN_I :
			base = base | T_SIN_P ;
			indic_t = 1 ;
			break ;

		    case T_SIN_P :
			base = base | T_SIN_I ;
			indic_t = 1 ;
			break ;

		    default :
			break ;
		}
		break ;
	    
	    case T_SIN_P :
		switch (b2_t) {
		    case T_SIN_P :
			base = base | T_COS_P ;
			indic_t = 1 ;
			break ;

		    case T_COS_P :
			base = base | T_SIN_P ;
			indic_t = 1 ;
			break ;

		    case T_COS_I :
			base = base | T_SIN_I ;
			indic_t = 1 ;
			break ;

                    case T_SIN_I :
			base = base | T_COS_I ;
			indic_t = 1 ;
			break ;

		    default :
			break ;
		}
		break ;

	    case T_SIN_I :
		switch (b2_t) {
		    case T_SIN_I :
			base = base | T_COS_P ;
			indic_t = 1 ;
			break ;

		    case T_COS_I :
			base = base | T_SIN_P ;
			indic_t = 1 ;
			break ;

		    case T_COS_P :
			base = base | T_SIN_I ;
			indic_t = 1 ;
			break ;

		    case T_SIN_P :
			base = base | T_COS_I ;
			indic_t = 1 ;
			break ;

		    default :
			break ;
		}
		break ;

	    case T_COSSIN_C :
		switch (b2_t) {
		    case T_COSSIN_C :
			base = base | T_COSSIN_C ;
			indic_t = 1 ;
			break ;
			
		    case T_COSSIN_S :
			base = base | T_COSSIN_S ;
			indic_t = 1 ;
			break ;
						
		    default :
			break ;
		}
		break ;

	     case T_COSSIN_S :
		switch (b2_t) {
		    case T_COSSIN_C :
			base = base | T_COSSIN_S ;
			indic_t = 1 ;
			break ;
			
		    case T_COSSIN_S :
			base = base | T_COSSIN_C ;
			indic_t = 1 ;
			break ;
						
		    default :
			break ;
		}
		break ;

	    case T_LEG_P :
		switch (b2_t) {
		    case T_LEG_P :
			base = base | T_LEG_P ;
			indic_t = 1 ;
			break ;
		    case T_LEG_I :
			base = base | T_LEG_I ;
			indic_t = 1 ; 
			break ;
		    default :
			break ;
		}
		break ;
	    
	    case T_COS :
		switch (b2_t) {
		    case T_COS :
			base = base | T_COS ;
			indic_t = 1 ;
			break ;
		    
		    case T_SIN :
			base = base | T_SIN ;
			indic_t = 1 ;
			break ;

		    default :
			break ;
		}
		break ;
	    
	    case T_SIN :
		switch (b2_t) {
		    case T_SIN :
			base = base | T_COS ;
			indic_t = 1 ;
			break ;

		    case T_COS :
			base = base | T_SIN ;
			indic_t = 1 ;
			break ;

		    default :
			break ;
		}
		break ;

	    case T_LEG_I :
		switch (b2_t) {
		    case T_LEG_P :
			base = base | T_LEG_I ;
			indic_t = 1 ;
			break ;
		    case T_LEG_I :
			base = base | T_LEG_P ;
			indic_t = 1 ; 
			break ;
		    default :
			break ;
		}
		break ; 
	    
	  
	    case T_LEG :
		switch (b2_t) {
		    case T_LEG :
			base = base | T_LEG ;
			indic_t = 1 ;
			break ;
		   
		    default :
			break ;
		}
		break ; 
	    
	    case T_LEG_MP :
		switch (b2_t) {
		    case T_LEG_MP :
			base = base | T_LEG_MP ;
			indic_t = 1 ;
			break ;
		   
		    case T_LEG_MI :
			base = base | T_LEG_MI ;
			indic_t = 1 ;
			break ;
		   
		    default :
			break ;
		}
		break ; 
	    
	    case T_LEG_MI :
		switch (b2_t) {
		    case T_LEG_MP :
			base = base | T_LEG_MI ;
			indic_t = 1 ;
			break ;

		    case T_LEG_MI :
			base = base | T_LEG_MP ;
			indic_t = 1 ;
			break ;
		   
		    default :
			break ;
		}
		break ; 
	    
	    
	    default :
		break ;
	}
	
	switch (b1_r) {
	    
	    case R_CHEB :
		switch (b2_r) {
		    case R_CHEB :
			base = base | R_CHEB ;
			indic_r = 1 ;
			break ;
		    
		    default :
			break ;
		}
		break ;

	    case R_LEG :
		switch (b2_r) {
		    case R_LEG :
			base = base | R_LEG ;
			indic_r = 1 ;
			break ;
		    
		    default :
			break ;
		}
		break ;

	    case R_JACO02 :
		switch (b2_r) {
		    case R_JACO02 :
			base = base | R_JACO02 ;
			indic_r = 1 ;
			break ;
		    
		    default :
			break ;
		}
		break ;
		
	    case R_CHEBU :
		switch (b2_r) {
		    case R_CHEBU :
			base = base | R_CHEBU ;
			indic_r = 1 ;
			break ;
			
		    default :
			break ;
		}
		break ;
		
	    case R_CHEBPIM_P :
		switch (b2_r) {
		    case R_CHEBPIM_P :
			base = base | R_CHEBPIM_P ;
			indic_r = 1 ;
			break ;
			
		    case R_CHEBPIM_I :
			base = base | R_CHEBPIM_I ;
			indic_r = 1 ;
			break ;
		    
		    default :
			break ;
		}
		break ;
		
	    case R_CHEBPIM_I :
		switch (b2_r) {
		    case R_CHEBPIM_P :
			base = base | R_CHEBPIM_I ;
			indic_r = 1 ;
			break ;
			
		    case R_CHEBPIM_I :
			base = base | R_CHEBPIM_P ;
			indic_r = 1 ;
			break ;
		    
		    default :
			break ;
		}
		break ;

	     case R_CHEBPI_I :
		switch (b2_r) {
		    case R_CHEBPI_P :
			base = base | R_CHEBPI_I ;
			indic_r = 1 ;
			break ;
			
		    case R_CHEBPI_I :
			base = base | R_CHEBPI_P ;
			indic_r = 1 ;
			break ;
		    
		    default :
			break ;
		}
		break ;

	     case R_CHEBPI_P :
		switch (b2_r) {
		    case R_CHEBPI_P :
			base = base | R_CHEBPI_P ;
			indic_r = 1 ;
			break ;
			
		    case R_CHEBPI_I :
			base = base | R_CHEBPI_I ;
			indic_r = 1 ;
			break ;
		    
		    default :
			break ;
		}
		break ;
		
	    case R_CHEBP :
		switch (b2_r) {
		    case R_CHEBP :
			base = base | R_CHEBP ;
			indic_r = 1 ;
			break ;
			
		    case R_CHEBI : 
			base = base | R_CHEBI ;
			indic_r = 1 ;
			break ;
			
		    default :
			break ;
		}
		break ;
		
	    case R_CHEBI :
		switch (b2_r) {
		    case R_CHEBP :
			base = base | R_CHEBI ;
			indic_r = 1 ;
			break ;
			
		    case R_CHEBI : 
			base = base | R_CHEBP ;
			indic_r = 1 ;
			break ;
			
		    default :
			break ;
		}
		break ;
	 
	    case R_LEGP :
		switch (b2_r) {
		    case R_LEGP :
			base = base | R_LEGP ;
			indic_r = 1 ;
			break ;
			
		    case R_LEGI : 
			base = base | R_LEGI ;
			indic_r = 1 ;
			break ;
			
		    default :
			break ;
		}
		break ;
		
	    case R_LEGI :
		switch (b2_r) {
		    case R_LEGP :
			base = base | R_LEGI ;
			indic_r = 1 ;
			break ;
			
		    case R_LEGI : 
			base = base | R_LEGP ;
			indic_r = 1 ;
			break ;
			
		    default :
			break ;
		}
		break ;
	 
	 default :
	    break ;   
	}
    
	if (indic_r*indic_t*indic_p == -1)
	    indic_total = -1 ;
	
	res.b[l] = base ;
    }
    
    if (indic_total == -1)
	res.set_base_nondef() ;
	
    return res ;    
}
