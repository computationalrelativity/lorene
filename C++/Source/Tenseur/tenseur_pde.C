/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char tenseur_pde_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/08/07 16:14:11  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.2  2002/07/09 16:46:23  p_grandclement
 * The Param in the case of an affine mapping is now 0x0 and not deleted
 * (I wonder why it was working before)
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.13  2000/10/04  14:58:32  eric
 * Ajout de shift.set_etat_qcq() avant l'affection de shift.
 *
 * Revision 2.12  2000/09/27  15:07:22  eric
 * Traitement source nulle dans poisson_vect.
 *
 * Revision 2.11  2000/05/22  15:48:18  phil
 * modification de Oohara pour passer avec dzpuis == 2
 * 		pardon 					3
 *
 * Revision 2.10  2000/05/22  15:00:46  phil
 * Modification de la methode de Shibata :
 * on doit passer une source en r^4 et l equation scalaire est alors resolue en utilisant l'algotihme avec dzpuis == 3
 *
 * Revision 2.9  2000/03/10  15:51:42  eric
 * Traitement dzpuis de source_scal.
 *
 * Revision 2.8  2000/03/08  10:21:05  eric
 * Appel de delete sur le Param* retourne par Map_et::donne_para_poisson_vect[
 * lorsqu'il n'est plus utilise (correction Memory leak).
 *
 * Revision 2.7  2000/03/07  16:53:42  eric
 * *** empty log message ***
 *
 * Revision 2.6  2000/03/07  15:43:32  phil
 * gestion des cas dzpuis ==4
 *
 * Revision 2.5  2000/02/21  12:55:09  eric
 * Traitement des triades.
 *
 * Revision 2.4  2000/02/16  17:13:05  eric
 * Correction
 *   mp->donne_para_poisson_vect(para, 4)
 * devient
 *   mp->donne_para_poisson_vect(para, 3)
 * dans Tenseur::poisson_vect
 * /
 *
 * Revision 2.3  2000/02/15  10:26:49  phil
 * le calcul de sol n'appelle plus Map::poisson_vect mais est fait dans
 * Tenseur::poisson_vect (respectivement poisson_vect_oohara)
 *
 * Revision 2.2  2000/02/09  19:32:58  eric
 * La triade de decomposition est desormais passee en argument des constructeurs.
 *
 * Revision 2.1  2000/02/09  10:01:32  phil
 * ajout version oohara
 *
 * Revision 2.0  2000/01/21  12:58:57  phil
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
void Tenseur::poisson_vect(double lambda, Param& para, Tenseur& shift
			    , Tenseur& vecteur, Tenseur& scalaire) const {
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

    for (int i=0 ; i<3 ; i++)
	assert ((*this)(i).check_dzpuis(4)) ;

    // On construit le tableau contenant le terme P_i ...
    for (int i=0 ; i<3 ; i++) {
	Param* par = mp->donne_para_poisson_vect(para, i) ; 

	(*this)(i).poisson(*par, vecteur.set(i)) ;

	if (par != 0x0)
	  delete par ; 
    }
    vecteur.set_triad( *triad ) ; 
    
    // Equation de Poisson scalaire :
    Tenseur source_scal (-skxk(*this)) ;
      
    assert (source_scal().check_dzpuis(3)) ; 

    Param* par = mp->donne_para_poisson_vect(para, 3) ; 

    source_scal().poisson(*par, scalaire.set()) ;
    
    if (par !=0x0)
      delete par ; 

    // On construit le tableau contenant le terme d xsi / d x_i ...
    Tenseur auxiliaire(scalaire) ;
    Tenseur dxsi (auxiliaire.gradient()) ;
    dxsi.dec2_dzpuis() ;
 
    // On construit le tableau contenant le terme x_k d P_k / d x_i
    Tenseur dp (skxk(vecteur.gradient())) ;
    dp.dec_dzpuis() ;
    
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
Tenseur Tenseur::poisson_vect(double lambda, Tenseur& vecteur, 
				    Tenseur& scalaire) const {
      
    Param bidon ;
    Tenseur resu(*mp, valence, type_indice, triad, metric, poids) ;
    resu.set_etat_qcq() ;
    poisson_vect(lambda, bidon, resu, vecteur, scalaire) ;
    return resu ;
}


		    //-----------------------------------//
		    //      Vectorial Poisson equation	 //
		    //      using Oohara scheme 	 //
		    //-----------------------------------//
		    
// Version avec parametres
// -----------------------
void Tenseur::poisson_vect_oohara(double lambda, Param& para, Tenseur& shift, 
			    Tenseur& chi) const {
    
     // Ne marche pas pour lambda =-1
    assert (lambda != -1) ;
    
    // Verifications d'usage ...
    assert (valence == 1) ;
    assert (shift.get_valence() == 1) ;
    assert (shift.get_type_indice(0) == type_indice(0)) ;
    assert (chi.get_valence() == 0) ;
    assert (etat != ETATNONDEF) ;

    // Nothing to do if the source is zero
    if (etat == ETATZERO) {
	shift.set_etat_zero() ; 
	chi.set_etat_qcq() ;
	chi.set() = 0 ; 
	return ; 
    }

    for (int i=0 ; i<3 ; i++)
	assert ((*this)(i).check_dzpuis(3) ||
		(*this)(i).check_dzpuis(4)) ;
    

    Tenseur copie(*this) ;
    copie.dec2_dzpuis() ;
    if ((*this)(0).check_dzpuis(4))
	copie.dec2_dzpuis() ;
    else
	copie.dec_dzpuis() ;
    
    Tenseur source_scal(contract(copie.gradient(), 0, 1)/(1.+lambda)) ;
    source_scal.inc2_dzpuis() ;
    
    Param* par = mp->donne_para_poisson_vect(para, 3) ; 
    
    source_scal().poisson(*par, chi.set());
    
    if (par !=0x0)
      delete par ; 
  
    Tenseur source_vect(*this) ;
    if ((*this)(0).check_dzpuis(4))
	source_vect.dec_dzpuis() ;
    
    Tenseur chi_grad (chi.gradient()) ;
    chi_grad.inc_dzpuis() ;
    
    for (int i=0 ; i<3 ; i++)
	source_vect.set(i) -= lambda*chi_grad(i) ;
	
    assert( *(source_vect.triad) == *((chi.gradient()).get_triad()) ) ;
    
    for (int i=0 ; i<3 ; i++) {
	par = mp->donne_para_poisson_vect(para, i) ;

	source_vect(i).poisson(*par, shift.set(i)) ;   

	if (par !=0x0)
	  delete par ; 
    }
    shift.set_triad( *(source_vect.triad) ) ; 

}


// Version sans parametres
// -----------------------
Tenseur Tenseur::poisson_vect_oohara(double lambda, Tenseur& scalaire) const {
      
    Param bidon ;
    Tenseur resu(*mp, valence, type_indice, triad, metric, poids) ;
    resu.set_etat_qcq() ;
    poisson_vect_oohara(lambda, bidon, resu, scalaire) ;
    return resu ;
}

