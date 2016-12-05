/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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
 * $Id$
 * $Log$
 * Revision 1.9  2016/12/05 16:18:02  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.8  2014/10/13 08:53:14  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:16:04  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2009/10/23 12:54:47  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.5  2009/10/13 19:45:01  j_novak
 * New base T_LEG_MP.
 *
 * Revision 1.4  2005/02/16 15:19:55  m_forot
 * Add the case T_LEG
 *
 * Revision 1.3  2003/09/16 12:11:59  j_novak
 * Added the base T_LEG_II.
 *
 * Revision 1.2  2002/10/16 14:36:57  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.8  2000/10/04  14:56:34  eric
 * nullite_plm_nonsym_anti : borne_sup est mise toujours egale a nt-2
 *   (et non plus a nt-1 dans le cas m pair).
 * Ajout des bases T_LEG_IP et T_LEG_PI (deja dans la version 2.7).
 *
 * Revision 2.7  2000/10/03  14:20:09  eric
 * *** empty log message ***
 *
 * Revision 2.6  1999/12/16  16:41:27  phil
 * *** empty log message ***
 *
 * Revision 2.5  1999/12/16  16:16:45  phil
 * correction cas nt = 1
 *
 * Revision 2.4  1999/09/16  12:05:51  phil
 * correction des cas antisym en z=0
 *
 * Revision 2.3  1999/09/14  17:52:47  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/09/14  17:41:48  phil
 * On commence l'ajout des cas antisymetriques en z=0
 *
 * Revision 2.1  1999/04/13  13:49:10  phil
 * *** empty log message ***
 *
 * Revision 2.0  1999/04/13  13:31:15  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Entetes C
#include <cstdlib>

// Entete Lorene
#include "headcpp.h"
#include "type_parite.h"
#include "base_val.h"


// fonction testant la nullite des fonctions de developpements
// j indice en theta -- nt nbre de points en theta
// k indice en phi   -- np nbre de points en phi

	 //-------------------------------------------------------
	// Developpement en P_COSSIN pour phi et T_LEG en theta
       //---------------------------------------------------------

namespace Lorene {
int nullite_plm_t_leg (int j, int nt, int k, int np) {

    int m = (k%2 == 0) ? k/2 : (k-1)/2 ;
    int borne_sup = nt-1 ;
    int borne_inf = m ;
    if ((j<borne_inf) || (j>borne_sup) || (k==1) || (k>np)) 
    return 0 ; else return 1 ; 
}


	 //-------------------------------------------------------
	// Developpement en P_COSSIN pour phi et T_LEG_P en theta
       //---------------------------------------------------------

int nullite_plm_nonsym (int j, int nt, int k, int np) {

    int m = (k%2 == 0) ? k/2 : (k-1)/2 ;
    int borne_sup = (m%2 == 0) ? nt-1 : nt-2 ;
    int borne_inf = (m%2 == 0) ? m/2 : (m-1)/2 ;
    if ((j<borne_inf) || (j>borne_sup) || (k==1) || (k>np)) 
    return 0 ; else return 1 ; 
}

	 //-------------------------------------------------------
	// Developpement en P_COSSIN pour phi et T_LEG_I en theta
       //---------------------------------------------------------

int nullite_plm_nonsym_anti (int j, int nt, int k, int np) {

    int m = (k%2 == 0) ? k/2 : (k-1)/2 ;

    int borne_sup = nt-2 ;

    int borne_inf = (m%2 == 0) ? m/2 : (m-1)/2 ;
    if ((j<borne_inf) || (j>borne_sup) || (k==1) || (k>np)) 
    return 0 ; else return 1 ; 
}




	 //------------------------------------------------------
	// Developpement en P_COSSIN_P pour phi et T_LEG_PP en theta
       //------------------------------------------------------

int nullite_plm_sym (int j, int nt, int k, int np) {
    
    int m = (k%2 == 0) ? k : k-1 ;
    int borne_inf = m/2 ;
    int borne_sup = nt-1 ;
    if ((j<borne_inf) || (j>borne_sup) || (k==1) || (k>np)) 
    return 0 ; else return 1 ; 
}


	 //-------------------------------------------------------
	// Developpement en P_COSSIN_P pour phi et T_LEG_IP en theta
       //---------------------------------------------------------

int nullite_plm_t_leg_ip(int j, int nt, int k, int np) {

    int m = (k%2 == 0) ? k : k-1 ;
    int borne_sup =  nt-2 ;
    int borne_inf =  m/2 ;
    if ((j<borne_inf) || (j>borne_sup) || (k==1) || (k>np)) 
    return 0 ; else return 1 ; 
}


	 //-------------------------------------------------------
	// Developpement en P_COSSIN_I pour phi et T_LEG_PI en theta
       //---------------------------------------------------------

int nullite_plm_t_leg_pi(int j, int nt, int k, int np) {

    int m ;
    if (k<=2) {
	m = 1 ; 
    }
    else{
	m = (k%2 == 0) ? k-1 : k ; 
    }
    
    int borne_sup =  nt-2 ;
    int borne_inf =  (m-1)/2 ;
    if ((j<borne_inf) || (j>borne_sup) || (k==1) || (k>np)) 
    return 0 ; else return 1 ; 
}

	 //-------------------------------------------------------
	// Developpement en P_COSSIN_I pour phi et T_LEG_II en theta
       //---------------------------------------------------------

int nullite_plm_t_leg_ii(int j, int nt, int k, int np) {

    int m ;
    if (k<=2) {
	m = 1 ; 
    }
    else{
	m = (k%2 == 0) ? k-1 : k ; 
    }
    
    int borne_sup =  nt-2 ;
    int borne_inf =  (m+1)/2 ;
    if ((j<borne_inf) || (j>borne_sup) || (k==1) || (k>np)) 
    return 0 ; else return 1 ; 
}

	 //----------------------------------------------------------
	// Developpement en P_COSSIN_P pour phi et T_LEG_MP en theta
       //------------------------------------------------------------

int nullite_plm_t_leg_mp (int j, int nt, int k, int np) {
    
    int m = (k%2 == 0) ? k : k-1 ;
    int borne_inf = m ;
    int borne_sup = nt-1 ;
    if ((j<borne_inf) || (j>borne_sup) || (k==1) || (k>np)) 
    return 0 ; else return 1 ; 
}


	 //----------------------------------------------------------
	// Developpement en P_COSSIN_P pour phi et T_LEG_MI en theta
       //------------------------------------------------------------

int nullite_plm_t_leg_mi (int j, int nt, int k, int np) {
    
    int m = 2*( (k-1) / 2) + 1 ;
    int borne_inf = m ;
    int borne_sup = nt-1 ;
    if ((j<borne_inf) || (j>borne_sup) || (k==1) || (k>np)) 
    return 0 ; else return 1 ; 
}




		//-----------------------------
	       //       La fonction 
	      //-------------------------------
	      
int nullite_plm (int j, int nt, int k, int np, Base_val base) {
    
    // on recupere les bases angulaires dans le noyau :
    // elles doivent etre identiques dans toutes les zones.
    
    int base_t = (base.b[0] & MSQ_T) ;
    int base_p = (base.b[0] & MSQ_P) ;
    int result ;
   
    switch (base_p) {
	case P_COSSIN :
	    // cas sym ou antisym en z=0 ...
	    switch (base_t) {
		case T_LEG_P :
		  result = nullite_plm_nonsym (j, nt, k, np) ;
		  break ;  
	    
	    
		case T_LEG_I :
		  result = nullite_plm_nonsym_anti (j, nt, k, np) ;
		  break ;  

		case T_LEG :
		  result = nullite_plm_t_leg (j, nt, k, np) ;
		  break ;  	  
	    
		default :
		    cout << "nullite_plm : cas inconnu ..." << endl ;
		    abort() ;
	    }
	    break ;

	case P_COSSIN_P :
	    switch (base_t) {
		case T_LEG_PP :
		  result = nullite_plm_sym (j, nt, k, np) ;
		  break ;  
	    
	    
		case T_LEG_IP :
		  result = nullite_plm_t_leg_ip (j, nt, k, np) ;
		  break ;  
	    
		case T_LEG_MP :
		  result = nullite_plm_t_leg_mp (j, nt, k, np) ;
		  break ;  
	    
		default :
		    cout << "nullite_plm : cas inconnu ..." << endl ;
		    abort() ;
	    }
	    break ;
	
	case P_COSSIN_I :
	    switch (base_t) {
		case T_LEG_PI :
		  result = nullite_plm_t_leg_pi (j, nt, k, np) ;
		  break ;  
	    
		case T_LEG_II :
		  result = nullite_plm_t_leg_ii (j, nt, k, np) ;
		  break ;  
	    
		case T_LEG_MI :
		  result = nullite_plm_t_leg_mi (j, nt, k, np) ;
		  break ;  
	    
		default :
		    cout << "nullite_plm : cas inconnu ..." << endl ;
		    abort() ;
	    }
	    break ;
	
	default :
	    cout << "nullite_plm : cas inconnu ..." << endl ;
	    abort() ;
    }
	
    return result ;
}
}
