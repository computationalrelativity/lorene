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


char valeur_ylm_i_C[] = "$Header$" ;

/*
 * Fonction membre de la classe Valeur qui calcule les coefficients
 * de la decomposition en cos(l*theta) / sin(l*theta)  
 * a partir des coefficients de la decomp. en harmoniques spheriques 
 * 
 */


/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:37:16  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.11  2000/09/29  16:10:37  eric
 * Ajout des bases T_LEG_IP et T_LEG_PI.
 *
 * Revision 2.10  2000/03/31  15:58:26  phil
 * changement des bases meme si etat est zero
 *
 * Revision 2.9  1999/12/22  16:25:52  eric
 * Traitement du cas ETATZERO
 * Test sur c_cf avant d'appeler coef().
 *
 * Revision 2.8  1999/12/16  16:41:43  phil
 * *** empty log message ***
 *
 * Revision 2.7  1999/12/16  16:10:07  phil
 * correction cas nt = 1.
 *
 * Revision 2.6  1999/11/30  12:46:40  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.5  1999/11/22  11:36:14  eric
 * Suppression de #include "tenseur.h"
 *
 * Revision 2.4  1999/10/18  14:12:39  eric
 * Les bases sont desormais membres des Mtbl_cf.
 *
 * Revision 2.3  1999/04/14  10:19:52  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/04/14  09:53:21  phil
 * remplacement de malloc en new
 *
 * Revision 2.1  1999/04/14  09:38:47  phil
 * Changement liberation memoire : free -> delete
 *
 * Revision 2.0  1999/04/13  16:47:28  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers C
#include <assert.h>
#include <stdlib.h>

// headers Lorene
#include "valeur.h"
#include "proto.h"

void ylm_i_pasprevu(const int*, const double*, double*) ;

void Valeur::ylm_i() {

    static void (*chbase_t[MAX_BASE])(const int*, const double*,
					        double*) ;
    static int nouv_base_t[MAX_BASE] ;   
    static int premier_appel = 1 ;
    
    int deg[3] ;
    int i, l ;

    if (premier_appel==1) {
	premier_appel = 0 ;

	for (i=0; i<MAX_BASE; i++) {
	    chbase_t[i] = ylm_i_pasprevu ;
	    nouv_base_t[i] = NONDEF ; 
	}
		
	chbase_t[T_LEG_P >> TRA_T] = chb_legp_cossincp ;
	nouv_base_t[T_LEG_P >> TRA_T] = T_COSSIN_CP  ;

	chbase_t[T_LEG_I >> TRA_T] = chb_legi_cossinci ;
	nouv_base_t[T_LEG_I >> TRA_T] = T_COSSIN_CI  ;

	chbase_t[T_LEG_PP >> TRA_T] = chb_legpp_cosp ;
	nouv_base_t[T_LEG_PP >> TRA_T] = T_COS_P  ;

	chbase_t[T_LEG_IP >> TRA_T] = chb_legip_cosi ;
	nouv_base_t[T_LEG_IP >> TRA_T] = T_COS_I  ;

	chbase_t[T_LEG_PI >> TRA_T] = chb_legpi_sini ;
	nouv_base_t[T_LEG_PI >> TRA_T] = T_SIN_I  ;
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
	    
	    int vbase_t_tra =  vbase_t  >> TRA_T ;
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
	    
	    if ((vbase_t != T_COSSIN_CP) && 
	    (vbase_t != T_COS_P) && (vbase_t != T_COSSIN_CI) ) 
					    { // cas ou le calcul est necessaire
	    
		int vbase_t_tra =  vbase_t  >> TRA_T ;
		assert(vbase_t_tra < MAX_BASE) ; 

// Nouvelle base : 
		base.b[l] = ( vbase_p | nouv_base_t[vbase_t_tra] ) | vbase_r ;
		if (get_mg()->get_nt(l)==1)
		    continue ;
//... tbl contenant les coefficients dans la zone l : 
		Tbl* cf =  c_cf->t[l]  ;

		if (cf->get_etat() == ETATZERO) continue ; // On ne fait rien si le tbl = 0 

//... resultat du calcul :
		double* resu = new double [cf->get_taille()] ;

		deg[0] = get_mg()->get_np(l) ;	   // nb. de degres de liberte en phi
		deg[1] = get_mg()->get_nt(l) ;	   // nb. de degres de liberte en theta
		deg[2] = get_mg()->get_nr(l) ;	   // nb. de degres de liberte en r
	    	    
	    
// Transformation en theta:
//-------------------------

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

void ylm_i_pasprevu(const int*, const double*, double*) {

    cout << 
     "Valeur::ylm_i: le changement de base demande n'est pas implemente !" 
     << endl ;
    abort() ; 
}


