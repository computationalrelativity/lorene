/*
 * Computation of the spectral coefficients.
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


char valeur_coef_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.14  2012/01/17 15:07:57  j_penner
 * using MAX_BASE_2 for the phi coordinate
 *
 * Revision 1.13  2009/10/23 12:56:29  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.12  2009/10/13 13:49:58  j_novak
 * New base T_LEG_MP.
 *
 * Revision 1.11  2008/10/07 15:01:58  j_novak
 * The case nt=1 is now treated separately.
 *
 * Revision 1.10  2008/05/24 15:09:02  j_novak
 * Getting back to previous version, the new one was an error.
 *
 * Revision 1.9  2008/05/24 15:05:22  j_novak
 * New method Scalar::match_tau to match the output of an explicit time-marching scheme with the tau method.
 *
 * Revision 1.8  2008/02/18 13:53:51  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.7  2007/12/11 15:28:25  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.6  2004/11/23 15:17:19  m_forot
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
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
 * Revision 2.10  2000/10/04  14:41:26  eric
 * Ajout des bases T_LEG_IP et T_LEG_PI
 *
 * Revision 2.9  2000/09/29  16:09:25  eric
 * Mise a zero des coefficients k=1 et k=2 dans le cas np=1.
 *
 * Revision 2.8  2000/09/07  15:14:30  eric
 * Ajout de la base P_COSSIN_I
 *
 * Revision 2.7  2000/08/16  10:32:53  eric
 * Suppression de Mtbl::dzpuis.
 * >> .
 *
 * Revision 2.6  1999/11/30  12:42:29  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.5  1999/11/24  16:06:13  eric
 * Ajout du test de l'admissibilite FFT des nombres de degres de liberte.
 *
 * Revision 2.4  1999/10/28  07:43:14  eric
 * Modif commentaires.
 *
 * Revision 2.3  1999/10/13  15:51:07  eric
 * Ajout de la base dans l'appel au constructeur de Mtbl_cf.
 *
 * Revision 2.2  1999/06/22  14:21:23  phil
 * Ajout de dzpuis
 *
 * Revision 2.1  1999/03/01  14:55:12  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/22  15:40:37  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
#include<cmath>

// Header Lorene
#include "mtbl.h"
#include "mtbl_cf.h"
#include "valeur.h"
#include "proto.h"

// Prototypage local
void pasprevu_r(const int*, const int*, double*, const int*, double*) ;
void pasprevu_t(const int*, const int*, double*, const int*, double*) ;
void pasprevu_p(const int* ,const int* ,  double* ) ;

void base_non_def_r(const int*, const int*, double*, const int*, double*) ;
void base_non_def_t(const int*, const int*, double*, const int*, double*) ;
void base_non_def_p(const int* ,const int* ,  double* ) ;

bool admissible_fft(int ) ; 

void Valeur::coef() const {
    
    // Variables statiques
    static void (*coef_r[MAX_BASE])(const int*, const int*, double*, const int*, double*) ;
    static void (*coef_t[MAX_BASE])(const int*, const int*, double*, const int*, double*) ;
    static void (*coef_p[MAX_BASE_2])(const int* ,const int* ,  double* ) ;
    static int premier_appel = 1 ;
    
    // Premier appel
    if (premier_appel) {
	premier_appel = 0 ;

	for (int i=0; i<MAX_BASE; i++) {
	    coef_r[i] = pasprevu_r ;
	    coef_t[i] = pasprevu_t ;
	    if(i%2 == 0){
	    coef_p[i/2] = pasprevu_p ;
	    }
	}	

	coef_r[NONDEF] = base_non_def_r ;
	coef_r[R_CHEB >> TRA_R] = cfrcheb ;	    
	coef_r[R_CHEBU >> TRA_R] = cfrcheb ;	    
	coef_r[R_CHEBP >> TRA_R] = cfrchebp ;	    
	coef_r[R_CHEBI >> TRA_R] = cfrchebi ;	    
	coef_r[R_CHEBPIM_P >> TRA_R] = cfrchebpimp ;	    
	coef_r[R_CHEBPIM_I >> TRA_R] = cfrchebpimi ;	    
	coef_r[R_CHEBPI_P >> TRA_R] = cfrchebpip ;	    
	coef_r[R_CHEBPI_I >> TRA_R] = cfrchebpii ;
	coef_r[R_JACO02 >> TRA_R] = cfrjaco02 ;

	coef_t[NONDEF] = base_non_def_t ;
	coef_t[T_COS >> TRA_T] = cftcos ;
	coef_t[T_SIN >> TRA_T] = cftsin ;
	coef_t[T_COS_P >> TRA_T] = cftcosp ;
	coef_t[T_COS_I >> TRA_T] = cftcosi ;
	coef_t[T_SIN_P >> TRA_T] = cftsinp ;
	coef_t[T_SIN_I >> TRA_T] = cftsini ;
	coef_t[T_COSSIN_CP >> TRA_T] = cftcossincp ;
	coef_t[T_COSSIN_SI >> TRA_T] = cftcossinsi ;
	coef_t[T_COSSIN_SP >> TRA_T] = cftcossinsp ;
	coef_t[T_COSSIN_CI >> TRA_T] = cftcossinci ;
	coef_t[T_COSSIN_S >> TRA_T] = cftcossins ;
	coef_t[T_COSSIN_C >> TRA_T] = cftcossinc ;
	coef_t[T_LEG_P >> TRA_T] = cftlegp ;
	coef_t[T_LEG_PP >> TRA_T] = cftlegpp ;
	coef_t[T_LEG_I >> TRA_T] = cftlegi ;
	coef_t[T_LEG_IP >> TRA_T] = cftlegip ;
	coef_t[T_LEG_PI >> TRA_T] = cftlegpi ;
	coef_t[T_LEG_II >> TRA_T] = cftlegii ;
	coef_t[T_LEG_MP >> TRA_T] = cftlegmp ;
	coef_t[T_LEG_MI >> TRA_T] = cftlegmi ;
	coef_t[T_LEG >> TRA_T] = cftleg ;

	coef_p[NONDEF] = base_non_def_p ;
	coef_p[P_COSSIN >> TRA_P] = cfpcossin ;
	coef_p[P_COSSIN_P >> TRA_P] = cfpcossin ;
	coef_p[P_COSSIN_I >> TRA_P] = cfpcossini ;

    }  // fin des operation de premier appel


	    //-------------------//
	    //  DEBUT DU CALCUL  //
	    //-------------------//

    // Tout null ?
    if (etat == ETATZERO) {
	return ;
    }
    
    // Protection
    assert(etat != ETATNONDEF) ;
        
    // Peut-etre rien a faire
    if (c_cf != 0x0) {
	return ;
    }
    
    // Il faut bosser
    assert(c != 0x0) ;		// ..si on peut
    c_cf = new Mtbl_cf(mg, base) ;
    c_cf->set_etat_qcq() ;
    
    // Boucles sur les zones
    int nz = mg->get_nzone() ;
    for (int l=0; l<nz; l++) {
	
	// Initialisation des valeurs de this->c_cf avec celle de this->c :
	const Tbl* f =  (c->t)[l]  ;
	Tbl* cf =  (c_cf->t)[l]  ;

	if (f->get_etat() == ETATZERO) {
	    cf->set_etat_zero() ;
	    continue ; // on ne fait rien si le tbl = 0  
	}

	cf->set_etat_qcq() ;
	    
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

	int dim[3] ;
	dim[0] = np_c ;
	dim[1] = nt_c ;
	dim[2] = nr_c ;

	int nrnt = nr * nt ;
	int nrnt_c = nr_c * nt_c ; 
	    
	for (int i=0; i<np ; i++) {
	    for (int j=0; j<nt ; j++) {
		for (int k=0; k<nr ; k++) {
		    int index = nrnt * i + nr * j + k ;
		    int index_c = nrnt_c * i + nr_c * j + k ;			
		    (cf->t)[index_c] = (f->t)[index] ;
		}
	    }
	}

	// On recupere les bases en r, theta et phi : 
	int base_r = ( base.b[l] & MSQ_R ) >> TRA_R ;
	int base_t = ( base.b[l] & MSQ_T ) >> TRA_T ;
	int base_p = ( base.b[l] & MSQ_P ) >> TRA_P ;
	int vbase_t = base.b[l] & MSQ_T ;
	int vbase_p = base.b[l] & MSQ_P ;

	assert(base_r < MAX_BASE) ; 
	assert(base_t < MAX_BASE) ; 
	assert(base_p < MAX_BASE) ; 

	// Transformation en phi 
	// ---------------------
	if ( np > 1 ) {
	    assert( admissible_fft(np) ) ; 
	    
	    coef_p[base_p]( deg,  dim, (cf->t) ) ;
	}
	else{	// Cas np=1 : mise a zero des coefficients k=1 et k=2 :
	    for (int i=nrnt; i<3*nrnt; i++) {
		cf->t[i] = 0 ; 
	    }	    
	}

	// Transformation en theta:
	// ------------------------

	if ( nt > 1 ) {
	    assert( admissible_fft(nt-1) ) ; 
	    bool pair = ( (vbase_t == T_LEG_PP) || (vbase_t == T_LEG_IP)
			  || (vbase_t == T_LEG_MP) ) ;
	    bool impair = ( (vbase_t == T_LEG_PI) || (vbase_t == T_LEG_II)
			    || (vbase_t == T_LEG_MI) ) ;

	    if ((pair && (vbase_p == P_COSSIN_I)) ||
		(impair && (vbase_p == P_COSSIN_P)) )
		  pasprevu_t(deg, dim, (cf->t), dim, (cf->t) ) ;
	    else	    
	      coef_t[base_t](deg, dim, (cf->t), dim, (cf->t)) ;
	}
	else {
	    if ((vbase_t == T_LEG_PP) || (vbase_t == T_LEG_PI) || 
		(vbase_t == T_LEG_IP) || (vbase_t == T_LEG_II) ||
		(vbase_t == T_LEG_P) || (vbase_t == T_LEG_I) ||
		(vbase_t == T_LEG) || (vbase_t == T_LEG_MP)
		|| (vbase_t == T_LEG_MI) ) {
	       
		      *c_cf->t[l] *=sqrt(2.) ;		
	    }
	}
	

	// Transformation en r:
	// --------------------
	if ( nr > 1 ) {
	    assert( admissible_fft(nr-1) ) ; 
	    coef_r[base_r](deg, dim, (cf->t), dim, (cf->t)) ;
	}	
	   
    }  // fin de la boucle sur les differentes zones
    
}

		    //------------------------//
		    // Les machins pas prevus //
		    //------------------------//

void pasprevu_r(const int*, const int*, double*, const int*, double*) {
    cout << "Valeur::coef: the required expansion basis in r " << endl ;
    cout << "  is not implemented !" << endl ;
    abort() ;
}

void pasprevu_t(const int*, const int*, double*, const int*, double*) {
    cout << "Valeur::coef: the required expansion basis in theta " << endl ;
    cout << "  is not implemented !" << endl ;
    abort() ;
}

void pasprevu_p(const int*, const int*, double*) {
    cout << "Valeur::coef: the required expansion basis in phi " << endl ;
    cout << "  is not implemented !" << endl ;
    abort() ;
}

void base_non_def_r(const int*, const int*, double*, const int*, double*) {
    cout << "Valeur::coef: the expansion basis in r is undefined !" << endl ;
    abort() ;
}

void base_non_def_t(const int*, const int*, double*, const int*, double*) {
    cout << "Valeur::coef: the expansion basis in theta is undefined !" << endl ;
    abort() ;
}

void base_non_def_p(const int*, const int*, double*) {
    cout << "Valeur::coef: the expansion basis in phi is undefined !" << endl ;
    abort() ;
}

