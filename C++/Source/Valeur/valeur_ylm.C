/*
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


char valeur_ylm_C[] = "$Header$" ;


/*
 * Fonction membre de la classe Valeur qui calcule les coefficients
 * de la decomposition spectrale en harmoniques spheriques  
 * a partir des coefficients en cos(l*theta) / sin(l*theta) 
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/09/17 12:30:22  j_novak
 * New checks for changing to T_LEG* bases.
 *
 * Revision 1.3  2003/09/16 08:54:09  j_novak
 * Addition of the T_LEG_II base (odd in theta, only for odd m) and the
 * transformation functions to and from the T_SIN_P base.
 *
 * Revision 1.2  2002/10/16 14:37:16  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.12  2000/09/29  16:10:21  eric
 * Ajout des bases T_LEG_IP et T_LEG_PI.
 *
 * Revision 2.11  2000/03/31  15:58:40  phil
 * changement des bases meme si etat est zero
 *
 * Revision 2.10  1999/12/22  16:25:26  eric
 * Traitement du cas ETATZERO
 * Test sur c_cf avant d'appeler coef().
 *
 * Revision 2.9  1999/12/16  16:41:38  phil
 * *** empty log message ***
 *
 * Revision 2.8  1999/12/16  16:09:56  phil
 * correction cas nt=1
 *
 * Revision 2.7  1999/11/30  12:46:33  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.6  1999/11/19  09:25:57  eric
 * *** empty log message ***
 *
 * Revision 2.5  1999/10/18  14:12:09  eric
 * Les bases sont desormais membres des Mtbl_cf.
 *
 * Revision 2.4  1999/04/14  10:19:46  phil
 * *** empty log message ***
 *
 * Revision 2.3  1999/04/14  09:52:28  phil
 * remplacement de malloc en new
 *
 * Revision 2.2  1999/04/14  09:36:03  phil
 * Changement liberation memoire : free -> delete
 *
 * Revision 2.1  1999/04/13  16:44:07  phil
 * *** empty log message ***
 *
 * Revision 2.0  1999/04/13  16:36:57  phil
 * *** empty log message ***
 *
 * Revision 1.1  1999/04/13  16:36:14  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */

// headers du C
#include <assert.h>
#include <stdlib.h>

// headers Lorene
#include "type_parite.h"
#include "valeur.h"
#include "proto.h"

void ylm_pasprevu(const int*, const double*, double*) ;

//******************************************************************************
void Valeur::ylm() {

    static void (*chbase_t[MAX_BASE])(const int*, const double*,
					        double*) ;
    static int nouv_base_t[MAX_BASE] ;   
    static int premier_appel = 1 ;
    
    int deg[3] ;
    int i, l ;

    if (premier_appel==1) {
	premier_appel = 0 ;

	for (i=0; i<MAX_BASE; i++) {
	    chbase_t[i] = ylm_pasprevu ;
	    nouv_base_t[i] = NONDEF ; 
	}
		
	chbase_t[T_COSSIN_CP >> TRA_T] = chb_cossincp_legp ;
	nouv_base_t[T_COSSIN_CP >> TRA_T] = T_LEG_P ;

	chbase_t[T_COSSIN_CI >> TRA_T] = chb_cossinci_legi ;
	nouv_base_t[T_COSSIN_CI >> TRA_T] = T_LEG_I ;

	chbase_t[T_COS_P >> TRA_T] = chb_cosp_legpp ;
	nouv_base_t[T_COS_P >> TRA_T] = T_LEG_PP ;

	chbase_t[T_COS_I >> TRA_T] = chb_cosi_legip ;
	nouv_base_t[T_COS_I >> TRA_T] = T_LEG_IP ;

	chbase_t[T_SIN_I >> TRA_T] = chb_sini_legpi ;
	nouv_base_t[T_SIN_I >> TRA_T] = T_LEG_PI ;

	chbase_t[T_SIN_P >> TRA_T] = chb_sinp_legii ;
	nouv_base_t[T_SIN_P >> TRA_T] = T_LEG_II ;
    }

//---------------------------------------------------------------------------
// fin des operation de premier appel 
//---------------------------------------------------------------------------

    // Tout null ?
    int nzone = get_mg()->get_nzone() ;
    if (etat == ETATZERO) {
	for (int l=0 ; l<nzone ; l++) {
	    int vbase_r = base.b[l] & MSQ_R  ;
	    int vbase_t = base.b[l] & MSQ_T  ;
	    int vbase_p = base.b[l] & MSQ_P  ;
	    int vbase_t_tra = vbase_t >> TRA_T ;
	    base.b[l] = ( vbase_p | nouv_base_t[vbase_t_tra] ) | vbase_r ;
	}
	return ;
    }
    
    // Protection
    assert(etat != ETATNONDEF) ;
        
    if (c_cf == 0x0) {
	coef() ;	 // The coefficients are required
    }

// Boucle sur les differentes zones
	
	for (l=0; l<nzone; l++) {
	
// On recupere les anciennes bases en r, phi et theta : 
	    int vbase_r = base.b[l] & MSQ_R  ;
	    int vbase_t = base.b[l] & MSQ_T  ;
	    int vbase_p = base.b[l] & MSQ_P  ;

	    if ((vbase_t != T_LEG_P) && (vbase_t != T_LEG_IP) &&
	    (vbase_t != T_LEG_PP) && (vbase_t != T_LEG_I) &&
	    (vbase_t != T_LEG_II) && (vbase_t != T_LEG_PI) ) 
		{ // cas ou le calcul est necessaire

		int vbase_t_tra = vbase_t >> TRA_T ;
		assert(vbase_t_tra < MAX_BASE) ; 
		bool pair = ( (vbase_t == T_COS_P) || (vbase_t == T_COS_I)) ;
		bool impair = ( (vbase_t == T_SIN_P) || (vbase_t == T_SIN_I)) ;

// Nouvelle base : 
		base.b[l] = ( vbase_p | nouv_base_t[vbase_t_tra] ) | vbase_r ;
		if (get_mg()->get_nt(l)==1)
		    continue ;
//... tbl contenant les coefficients dans la zone l : 
		Tbl* cf =  c_cf->t[l]  ;

		if (cf->get_etat() == ETATZERO) continue ; // On ne fait rien si le tbl = 0 

		deg[0] = get_mg()->get_np(l) ;	   // nb. de degres de liberte en phi
		deg[1] = get_mg()->get_nt(l)  ;	   // nb. de degres de liberte en theta
		deg[2] = get_mg()->get_nr(l) ;	   // nb. de degres de liberte en r

//... resultat du calcul :
		double* resu = new double [cf->get_taille()] ;
	    	    	    
// Transformation en theta:
//-------------------------

		if ((pair && (vbase_p == P_COSSIN_I)) ||
		    (impair && (vbase_p == P_COSSIN_P)) )
		  ylm_pasprevu(deg, (cf->t), resu ) ;
		else
		  chbase_t[vbase_t_tra](deg, (cf->t), resu ) ;
	    
// On branche le tbl contenant les coef sur resu : 
		delete [] cf->t ;	// les anciens coef. sont oublies
		cf->t = resu ;	// nouveaux coef.

	    }	// fin du cas ou la transformation devait etre effectuee
	    
    }  // fin de la boucle sur les differentes zones

    // On met les bonnes bases dans c_cf : 
    c_cf->base = base ; 

}

//******************************************************************************

void ylm_pasprevu(const int* , const double*, double* ) {

    cout << 
     "Valeur::ylm: le changement de base demande n'est pas implemente !" 
     << endl ;
    abort() ; 
}

