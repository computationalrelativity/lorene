/*
 *  Methods of the class tenseur for solving vectorial Poisson equations
 *   with a falloff condition at the outer boundary
 *
 *    (see file tenseur.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Joshua A. Faber
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char tenseur_pde_falloff_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/11/30 20:55:33  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "cmp.h"
#include "param.h"
#include "tenseur.h"

		    //-----------------------------------//
		    //      Vectorial Poisson equation	 //
		    //-----------------------------------//

// Version avec parametres
// -----------------------
void Tenseur::poisson_vect_falloff(double lambda, Param& para, Tenseur& shift
			    , Tenseur& vecteur, Tenseur& scalaire, int k_falloff) const {
    assert (lambda != -1) ;
    
    // Verifications d'usage ...
    assert (valence == 1) ;
    assert (shift.get_valence() == 1) ;
    assert (shift.get_type_indice(0) == type_indice(0)) ;
    assert (vecteur.get_valence() == 1) ;
    assert (vecteur.get_type_indice(0) == type_indice(0)) ;
    assert (scalaire.get_valence() == 0) ;
    assert (etat != ETATNONDEF) ;

    // Nothing to do if the source is zero
    if (etat == ETATZERO) {

	shift.set_etat_zero() ; 

	vecteur.set_etat_qcq() ;
	for (int i=0; i<3; i++) {
	    vecteur.set(i) = 0 ; 
	}

	scalaire.set_etat_qcq() ;
	scalaire.set() = 0 ;  

	return ; 
    }

    // On construit le tableau contenant le terme P_i ...
    for (int i=0 ; i<3 ; i++) {
	Param* par = mp->donne_para_poisson_vect(para, i) ; 

	(*this)(i).poisson_falloff(*par, vecteur.set(i),k_falloff) ;

	if (par != 0x0)
	  delete par ; 
    }
    vecteur.set_triad( *triad ) ; 
    
    // Equation de Poisson scalaire :
    Tenseur source_scal (-skxk(*this)) ;
      
    Param* par = mp->donne_para_poisson_vect(para, 3) ; 

    source_scal().poisson_falloff(*par, scalaire.set(), k_falloff) ;
    
    if (par !=0x0)
      delete par ; 

    // On construit le tableau contenant le terme d xsi / d x_i ...
    Tenseur auxiliaire(scalaire) ;
    Tenseur dxsi (auxiliaire.gradient()) ;
 
    // On construit le tableau contenant le terme x_k d P_k / d x_i
    Tenseur dp (skxk(vecteur.gradient())) ;
    
    // Il ne reste plus qu'a tout ranger dans P :
    // The final computation is done component by component because
    // d_khi and x_d_w are covariant comp. whereas w_shift is
    // contravariant

    shift.set_etat_qcq() ; 

    for (int i=0 ; i<3 ; i++)
	shift.set(i) = (lambda+2)/2/(lambda+1) * vecteur(i) 
			    - (lambda/2/(lambda+1)) * (dxsi(i) + dp(i)) ;   
			    
    shift.set_triad( *(vecteur.triad) ) ; 

}


// Version sans parametres
// -----------------------
Tenseur Tenseur::poisson_vect_falloff(double lambda, Tenseur& vecteur, 
				    Tenseur& scalaire, int k_falloff) const {
      
    Param bidon ;
    Tenseur resu(*mp, valence, type_indice, triad, metric, poids) ;
    resu.set_etat_qcq() ;
    poisson_vect_falloff(lambda, bidon, resu, vecteur, scalaire, k_falloff) ;
    return resu ;
}

