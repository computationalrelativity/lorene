/*
 *  Computations of the values at the collocation points from the spectral
 *   coefficients
 *
 */

/*
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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


char valeur_coef_i_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.15  2013/06/06 15:31:33  j_novak
 * Functions to compute Legendre coefficients (not fully tested yet).
 *
 * Revision 1.14  2012/01/17 15:08:02  j_penner
 * using MAX_BASE_2 for the phi coordinate
 *
 * Revision 1.13  2009/10/23 12:56:29  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.12  2009/10/13 13:49:58  j_novak
 * New base T_LEG_MP.
 *
 * Revision 1.11  2009/10/08 16:23:14  j_novak
 * Addition of new bases T_COS and T_SIN.
 *
 * Revision 1.10  2008/10/07 15:01:58  j_novak
 * The case nt=1 is now treated separately.
 *
 * Revision 1.9  2007/12/11 15:28:25  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.8  2004/11/23 15:17:19  m_forot
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.7  2004/08/24 09:14:52  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.6  2003/10/13 20:51:25  e_gourgoulhon
 * Replaced malloc by new
 *
 * Revision 1.5  2003/09/17 12:30:22  j_novak
 * New checks for changing to T_LEG* bases.
 *
 * Revision 1.4  2003/09/16 13:07:41  j_novak
 * New files for coefficient trnasformation to/from the T_LEG_II base.
 *
 * Revision 1.3  2002/11/12 17:44:35  j_novak
 * Added transformation function for T_COS basis.
 *
 * Revision 1.2  2002/10/16 14:37:15  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.9  2000/10/04  14:41:40  eric
 * Ajout des bases T_LEG_IP et T_LEG_PI
 *
 * Revision 2.8  2000/09/07  15:14:44  eric
 * Ajout de la base P_COSSIN_I
 *
 * Revision 2.7  2000/08/16  10:33:04  eric
 * Suppression de Mtbl::dzpuis.
 *
 * Revision 2.6  1999/12/16  16:41:48  phil
 * *** empty log message ***
 *
 * Revision 2.5  1999/11/30  12:43:03  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.4  1999/10/28  07:44:20  eric
 * Modif commentaires.
 *
 * Revision 2.3  1999/10/13  15:51:33  eric
 * Anglisation.
 *
 * Revision 2.2  1999/06/22  15:10:02  phil
 * ajout de dzpuis
 *
 * Revision 2.1  1999/02/22  15:40:25  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

#include <math.h>

// Header Lorene
#include "valeur.h"
#include "proto.h"

void c_est_pas_fait(char * ) ;

void ipasprevu_r(const int*, const int*, double*, const int*, double*) ;
void ipasprevu_t(const int*, const int*, double*, const int*, double*) ;
void ipasprevu_p(const int* , const int* , const int* , double* , double* ) ;
void ibase_non_def_r(const int*, const int*, double*, const int*, double*) ;
void ibase_non_def_t(const int*, const int*, double*, const int*, double*) ;
void ibase_non_def_p(const int* , const int* , const int* , double* , double* ) ;

void Valeur::coef_i() const {
    
    // Variables statiques
    static void (*invcf_r[MAX_BASE])(const int*, const int*, double*, const int*, double*) ;
    static void (*invcf_t[MAX_BASE])(const int*, const int*, double*, const int*, double*) ;
    static void (*invcf_p[MAX_BASE_2])(const int* , const int* , const int*, double* , double* ) ;
    static int premier_appel = 1 ;

    // Premier appel
    if (premier_appel) {
	premier_appel = 0 ;

	for (int i=0; i<MAX_BASE; i++) {
	    invcf_r[i] = ipasprevu_r ;
	    invcf_t[i] = ipasprevu_t ;
	    if(i%2==0){
	    invcf_p[i/2] = ipasprevu_p ; // saves a loop
	    }
	}	

	invcf_r[NONDEF] = ibase_non_def_r ;
	invcf_r[R_CHEB >> TRA_R] = circheb ;	    
	invcf_r[R_CHEBU >> TRA_R] = circheb ;	    
	invcf_r[R_CHEBP >> TRA_R] = circhebp ;	    
	invcf_r[R_CHEBI >> TRA_R] = circhebi ;	    
	invcf_r[R_CHEBPIM_P >> TRA_R] = circhebpimp ;	    
	invcf_r[R_CHEBPIM_I >> TRA_R] = circhebpimi ;	    
	invcf_r[R_CHEBPI_P >> TRA_R] = circhebpip ;	    
	invcf_r[R_CHEBPI_I >> TRA_R] = circhebpii ;	    
	invcf_r[R_LEG >> TRA_R] = cirleg ;	    
	invcf_r[R_LEGP >> TRA_R] = cirlegp ;	    
	invcf_r[R_LEGI >> TRA_R] = cirlegi ;	    
	invcf_r[R_JACO02 >> TRA_R] = cirjaco02 ;

	invcf_t[NONDEF] = ibase_non_def_t ;
	invcf_t[T_COS >> TRA_T] = citcos ;
	invcf_t[T_SIN >> TRA_T] = citsin ;
	invcf_t[T_COS_P >> TRA_T] = citcosp ;
	invcf_t[T_COS_I >> TRA_T] = citcosi ;
	invcf_t[T_SIN_P >> TRA_T] = citsinp ;
	invcf_t[T_SIN_I >> TRA_T] = citsini ;
	invcf_t[T_COSSIN_CP >> TRA_T] = citcossincp ;
	invcf_t[T_COSSIN_SI >> TRA_T] = citcossinsi ;
	invcf_t[T_COSSIN_SP >> TRA_T] = citcossinsp ;
	invcf_t[T_COSSIN_CI >> TRA_T] = citcossinci ;
	invcf_t[T_COSSIN_S >> TRA_T] = citcossins ;
	invcf_t[T_COSSIN_C >> TRA_T] = citcossinc ;
	invcf_t[T_LEG_P >> TRA_T] = citlegp ;
	invcf_t[T_LEG_PP >> TRA_T] = citlegpp ;
	invcf_t[T_LEG_I >> TRA_T] = citlegi ;
	invcf_t[T_LEG_IP >> TRA_T] = citlegip ;
	invcf_t[T_LEG_PI >> TRA_T] = citlegpi ;
	invcf_t[T_LEG_II >> TRA_T] = citlegii ;	
	invcf_t[T_LEG_MP >> TRA_T] = citlegmp ;
	invcf_t[T_LEG_MI >> TRA_T] = citlegmi ;
	invcf_t[T_LEG >> TRA_T] = citleg ;

	invcf_p[NONDEF] = ibase_non_def_p ;
	invcf_p[P_COSSIN >> TRA_P] = cipcossin ;
	invcf_p[P_COSSIN_P >> TRA_P] = cipcossin ;
	invcf_p[P_COSSIN_I >> TRA_P] = cipcossini ;

    }  // fin des operation de premier appel

	    //------------------//
	    // DEBUT DU CALCUL  //
	    //------------------//	    
	   
    // Tout null ?
    if (etat == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(etat != ETATNONDEF) ;
        
    // Peut-etre rien a faire
    if (c != 0x0) {
	return ;
    }
    
    // Il faut bosser
    assert(c_cf != 0x0) ;		// ..si on peut
    assert(c_cf->base == base) ;		// Consistence des bases

    c = new Mtbl(mg) ;
    c->set_etat_qcq() ;
    
    // Boucles sur les zones
    int nz = mg->get_nzone() ;
    for (int l=0; l<nz; l++) {
	
	// Initialisation des valeurs de this->c_cf avec celle de this->c :
	Tbl* f =  (c->t)[l]  ;
	const Tbl* cf =  (c_cf->t)[l]  ;

	if (cf->get_etat() == ETATZERO) {
	    f->set_etat_zero() ;
	    continue ; // on ne fait rien si le tbl(cf) = 0  
	}

	f->set_etat_qcq() ;

	int np = f->get_dim(2) ;
	int nt = f->get_dim(1) ;
	int nr = f->get_dim(0) ;

	int np_c = cf->get_dim(2) ;
	int nt_c = cf->get_dim(1) ;
	int nr_c = cf->get_dim(0) ;	    

	// Attention a ce qui suit... (deg et dim)
	int deg[3] ;
	deg[0] = np ;
	deg[1] = nt ;
	deg[2] = nr ;

	int dimc[3] ;
	dimc[0] = np_c ;
	dimc[1] = nt_c ;
	dimc[2] = nr_c ;

	// Allocation de l'espace memoire pour le tableau de travail trav
	int ntot = cf->get_taille() ;
	double* trav = new double[ntot] ; 

	// On recupere les bases en r, theta et phi : 
	int base_r = ( base.b[l] & MSQ_R ) >> TRA_R ;
	int base_t = ( base.b[l] & MSQ_T ) >> TRA_T ;
	int base_p = ( base.b[l] & MSQ_P ) >> TRA_P ;
	int vbase_t = base.b[l] & MSQ_T ;
	int vbase_p = base.b[l] & MSQ_P ;

	assert(base_r < MAX_BASE) ; 
	assert(base_t < MAX_BASE) ; 
	assert(base_p < MAX_BASE_2) ; 

	// Transformation inverse en r:
	if ( nr == 1 ) {
	    for (int i=0; i<ntot; i++) {
		trav[i] = cf->t[i] ;	// simple recopie cf --> trav
	    }	    
	}
	else {
	    invcf_r[base_r]( deg, dimc, (cf->t), dimc, trav ) ;
	}
	
	// Partie angulaire
	if ( np == 1) {
	    if (nt==1) {
		for (int i=0 ; i<f->get_taille() ; i++)
		    f->t[i] = trav[i] ;
		if ((vbase_t == T_LEG_PP) || (vbase_t == T_LEG_PI) || 
		    (vbase_t == T_LEG_IP) || (vbase_t == T_LEG_II) ||
		    (vbase_t == T_LEG_P) || (vbase_t == T_LEG_I) ||
		    (vbase_t == T_LEG_MP) || (vbase_t == T_LEG_MI) ||
		    (vbase_t == T_LEG) ) {
		    
		    *f /=sqrt(2.) ;		
		}
	    }
	    
	    else {
	    bool pair = ( (vbase_t == T_LEG_PP) || (vbase_t == T_LEG_IP)
			  || (vbase_t == T_LEG_MP) ) ;
	    bool impair = ( (vbase_t == T_LEG_PI) || (vbase_t == T_LEG_II)
			  || (vbase_t == T_LEG_MI)) ;
	    
	    if ((pair && (vbase_p == P_COSSIN_I)) ||
		(impair && (vbase_p == P_COSSIN_P)) )
	      ipasprevu_t(deg, dimc, trav, deg, (f->t) ) ;
	    else	    
	      invcf_t[base_t]( deg, dimc, trav, deg, (f->t) ) ;
	  }
	}
	else {
	  // Cas 3-D
	  //  ...... Transformation inverse en theta:
	  bool pair = ( (vbase_t == T_LEG_PP) || (vbase_t == T_LEG_IP)
			|| (vbase_t == T_LEG_MP) ) ;
	  bool impair = ( (vbase_t == T_LEG_PI) || (vbase_t == T_LEG_II)
			  || (vbase_t == T_LEG_MI) ) ;
	  
	  if ((pair && (vbase_p == P_COSSIN_I)) ||
	      (impair && (vbase_p == P_COSSIN_P)) )
	    ipasprevu_t(deg, dimc, trav, dimc, trav ) ;
	  else	    
	    invcf_t[base_t]( deg, dimc, trav, dimc, trav ) ;
	  //  ...... Transformation inverse en phi:
	  invcf_p[base_p]( deg, dimc, deg, trav, (f->t) ) ;
	}
	// Menage
	delete [] trav ;
    }  // fin de la boucle sur les differentes zones

}

		    //------------------------//
		    // Les machins pas prevus //
		    //------------------------//

void ipasprevu_r(const int*, const int*, double*, const int*, double*) {
    cout << "Valeur::coef_i: the required expansion basis in r " << endl ;
    cout << "  is not implemented !" << endl ;
    abort() ; 
}

void ipasprevu_t(const int*, const int*, double*, const int*, double* ) {
    cout << "Valeur::coef_i: the required expansion basis in theta " << endl ;
    cout << "  is not implemented !" << endl ;
    abort() ; 
}

void ipasprevu_p(const int*, const int*, const int*, double*, double* ) {
    cout << "Valeur::coef_i: the required expansion basis in phi " << endl ;
    cout << "  is not implemented !" << endl ;
    abort() ; 
}

void ibase_non_def_r(const int*, const int*, double*, const int*, double*) {
    cout << "Valeur::coef_i: the expansion basis in r is undefined !" << endl ;
    abort() ; 
}

void ibase_non_def_t(const int*, const int*, double*, const int*, double*) {
    cout << "Valeur::coef_i: the expansion basis in theta is undefined !" 
	 << endl ;
    abort() ; 
}

void ibase_non_def_p(const int*, const int*, const int*, double*, double*) {
    cout << "Valeur::coef_i: the expansion basis in phi is undefined !" << endl ;
    abort() ; 
}

