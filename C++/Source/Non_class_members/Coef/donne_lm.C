/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Jerome Novak
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


char donne_lm_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.4  2004/11/23 15:13:50  m_forot
 * Added the bases for the cases without any equatorial symmetry
 * (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.3  2003/09/16 12:11:59  j_novak
 * Added the base T_LEG_II.
 *
 * Revision 1.2  2002/10/16 14:36:54  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 3.0  2000/10/09  09:15:01  novak
 * correction pour les cas de bases en r alternees
 *
 * Revision 2.8  2000/10/04  14:55:33  eric
 * Ajout des bases T_LEG_IP et T_LEG_PI.
 *
 * Revision 2.7  1999/12/16  16:41:33  phil
 * *** empty log message ***
 *
 * Revision 2.6  1999/12/16  16:23:39  phil
 * vire un assert
 *
 * Revision 2.5  1999/12/16  16:21:38  phil
 * correction cas nt = 1
 *
 * Revision 2.4  1999/09/16  12:06:11  phil
 * correction des cas antisymetriques en z=0
 *
 * Revision 2.3  1999/09/14  17:52:59  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/09/14  17:44:16  phil
 * *** empty log message ***
 *
 * Revision 2.1  1999/04/13  13:50:01  phil
 * *** empty log message ***
 *
 * Revision 2.0  1999/04/13  13:31:01  phil
 * *** empty log message ***
 *
 * Revision 1.1  1999/04/13  13:30:28  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Entetes C
#include <stdlib.h>

// Entete Lorene
#include "headcpp.h"
#include "type_parite.h"
#include "base_val.h"

/*
 * Fonction affection les nombres l_quant, m_quant et la base en r
 * 
 * ENTREES :	  nz : le nombre de zones
 *		  zone : la zone de travail  
 *		  j et k : indices en theta et phi, respectivement
 *		  base : la base de developpement
 * 
 * SORTIES : les variables m_quant, l_quant et base_r.
 * 
 */

//-----------------------------------------------------------------
// Developpement en P_COSSIN pour phi et T_LEG en theta
//-------------------------------------------------------------------

void donne_lm_nonsymTP (int j, int k, int &m_quant, int &l_quant) {

    m_quant = (k%2 == 0) ? k/2 : (k-1)/2;
    l_quant = j ;
    
}


	 //-----------------------------------------------------------------
	// Developpement en P_COSSIN pour phi et T_LEG_P en theta
       //-------------------------------------------------------------------

void donne_lm_nonsym (int j, int k, int &m_quant, int &l_quant) {

    m_quant = (k%2 == 0) ? k/2 : (k-1)/2;
    l_quant = (m_quant%2 == 0) ? 2*j : 2*j+1 ;
    
}

	 //-----------------------------------------------------------------
	// Developpement en P_COSSIN pour phi et T_LEG_I en theta
       //-------------------------------------------------------------------

void donne_lm_nonsym_anti (int j, int k, int &m_quant, int &l_quant) {

    m_quant = (k%2 == 0) ? k/2 : (k-1)/2;
    l_quant = (m_quant%2 == 1) ? 2*j : 2*j+1 ;
    
}

	 //------------------------------------------------------
	// Developpement en P_COSSIN_P pour phi et T_LEG_PP en theta
       //-------------------------------------------------------

void donne_lm_sym (int j, int k, int &m_quant, int &l_quant) {

    m_quant = (k%2 == 0) ? k : k-1;
    l_quant = 2*j ;
    
}

    
	 //-------------------------------------------------------
	// Developpement en P_COSSIN_P pour phi et T_LEG_IP en theta
       //---------------------------------------------------------

void donne_lm_t_leg_ip (int j, int k, int &m_quant, int &l_quant) {

    m_quant = (k%2 == 0) ? k : k-1 ;
    l_quant = 2*j+1 ;
    
}


	 //-------------------------------------------------------
	// Developpement en P_COSSIN_I pour phi et T_LEG_PI en theta
       //---------------------------------------------------------

void donne_lm_t_leg_pi (int j, int k, int &m_quant, int &l_quant) {

    if (k<=2) {
	m_quant = 1 ; 
    }
    else{
	m_quant = (k%2 == 0) ? k-1 : k ; 
    }
    
    l_quant = 2*j+1 ;
    
}

	 //-------------------------------------------------------
	// Developpement en P_COSSIN_I pour phi et T_LEG_II en theta
       //---------------------------------------------------------

void donne_lm_t_leg_ii (int j, int k, int &m_quant, int &l_quant) {

    if (k<=2) {
	m_quant = 1 ; 
    }
    else{
	m_quant = (k%2 == 0) ? k-1 : k ; 
    }
    
    l_quant = 2*j ;
    
}



		//-----------------------------
	       //       La fonction 
	      //-------------------------------
	      
void donne_lm (int nz, int zone, int j, int k, Base_val base,
		 int &m_quant, int &l_quant, int& base_r) {
    
    //verifications :
    assert (zone >= 0) ;
    assert (zone < nz) ;
    
    int base_t = (base.b[zone] & MSQ_T) ;
    int base_p = (base.b[zone] & MSQ_P) ;
    base_r = (base.b[zone] & MSQ_R) ;
    
    switch (base_p) {
	case P_COSSIN : 
	    // cas sym ou antisym en z=0 ...
	    switch (base_t) {

	      	case T_LEG :
		  donne_lm_nonsymTP (j, k, m_quant, l_quant) ;
		  break ;  

		case T_LEG_P :
		  donne_lm_nonsym (j, k, m_quant, l_quant) ;
		  break ;  
	    
		case T_LEG_I :
		  donne_lm_nonsym_anti (j, k, m_quant, l_quant) ;
		  break ;  
	    
		default :
		    cout << "donne_lm : cas inconnu ..." << endl ;
		    abort() ;
		    break ;
	    }
	    break ;

	case P_COSSIN_P :
	    switch (base_t) {

		case T_LEG_PP :
		    donne_lm_sym (j, k, m_quant, l_quant) ;
		    break ; 

		case T_LEG_IP :
		  donne_lm_t_leg_ip (j, k, m_quant, l_quant); 
		    break ; 

		default :
		    cout << "donne_lm : cas inconnu ..." << endl ;
		    abort() ;
		    break ;
	    }
	    break ; 

	case P_COSSIN_I :
	    switch (base_t) {

		case T_LEG_PI :
		    donne_lm_t_leg_pi (j, k, m_quant, l_quant) ;
		    break ; 

		case T_LEG_II :
		    donne_lm_t_leg_ii (j, k, m_quant, l_quant) ;
		    break ; 

		default :
		    cout << "donne_lm : cas inconnu ..." << endl ;
		    abort() ;
		    break ;

	    }
	    break ;
	
	default :
	    cout << "donne_lm : cas inconnu ..." << endl ;
	    abort() ;
	    break ;
    }
    switch (base_r) {

        case R_CHEBPI_P :
	  base_r = (l_quant%2 == 0) ? R_CHEBP : R_CHEBI ;
	  break ;
	  
        case R_CHEBPI_I :
	  base_r = (l_quant%2 == 1) ? R_CHEBP : R_CHEBI ;
	  break ;
	  
        case R_CHEBPIM_P :
	  base_r = (m_quant%2 == 0) ? R_CHEBP : R_CHEBI ;
	  break ;
	  
        case R_CHEBPIM_I :
	  base_r = (m_quant%2 == 1) ? R_CHEBP : R_CHEBI ;
	  break ;
	  
    }
}
