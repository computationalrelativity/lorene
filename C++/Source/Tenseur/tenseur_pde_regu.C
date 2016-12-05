/*
 *  Method of class Tenseur to solve a vector Poisson equation
 *  by regularizing its source.
 *
 *  (see file tenseur.h for documentation).
 *
 */

/*
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
 *   Copyright (c) 2001 Philippe Grandclement
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
 * Revision 1.5  2016/12/05 16:18:17  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:42  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2003/10/03 15:58:51  j_novak
 * Cleaning of some headers
 *
 * Revision 1.2  2002/08/07 16:14:11  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2001/01/15  11:01:34  phil
 * vire version sans parametres
 *
 * Revision 2.0  2000/10/06  15:34:03  keisuke
 * Initial revision.
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "param.h"
#include "tenseur.h"

		    //-----------------------------------//
		    //      Vectorial Poisson equation	 //
		    //-----------------------------------//

// Version avec parametres
// -----------------------
namespace Lorene {
void Tenseur::poisson_vect_regu(int k_div, int nzet, double unsgam1,
				double lambda, Param& para, Tenseur& shift,
				Tenseur& vecteur, Tenseur& scalaire) const {
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

    Tenseur vecteur_regu(*mp, 1, CON, mp->get_bvect_cart(), metric, poids) ;
    Tenseur vecteur_div(*mp, 1, CON, mp->get_bvect_cart(), metric, poids) ;
    Tenseur dvect_div(*mp, 1, CON, mp->get_bvect_cart(), metric, poids) ;
    Tenseur souvect_regu(*mp, 1, CON, mp->get_bvect_cart(), metric, poids) ;
    Tenseur souvect_div(*mp, 1, CON, mp->get_bvect_cart(), metric, poids) ;

    vecteur_regu.set_etat_qcq() ;
    vecteur_div.set_etat_qcq() ;
    dvect_div.set_etat_qcq() ;
    souvect_regu.set_etat_qcq() ;
    souvect_div.set_etat_qcq() ;

    // On construit le tableau contenant le terme P_i ...

    // Apply only to x and y components because poisson_regular does not
    // work for z component due to the symmetry.
    for (int i=0 ; i<2 ; i++) {
	Param* par = mp->donne_para_poisson_vect(para, i) ; 

	(*this)(i).poisson_regular(k_div, nzet, unsgam1, *par,
				   vecteur.set(i),
				   vecteur_regu.set(i), vecteur_div.set(i),
				   dvect_div,
				   souvect_regu.set(i), souvect_div.set(i)) ;

	delete par ; 
    }

    Param* par = mp->donne_para_poisson_vect(para, 2) ;

    (*this)(2).poisson(*par, vecteur.set(2)) ;

    delete par ;

    vecteur.set_triad( *triad ) ; 
    
    // Equation de Poisson scalaire :
    Tenseur source_scal (-skxk(*this)) ;
      
    assert (source_scal().check_dzpuis(3)) ; 

    par = mp->donne_para_poisson_vect(para, 3) ; 

    Tenseur scalaire_regu(*mp, metric, poids) ;
    Tenseur scalaire_div(*mp, metric, poids) ;
    Tenseur dscal_div(*mp, 1, CON, mp->get_bvect_cart(), metric, poids) ;
    Cmp souscal_regu(mp) ;
    Cmp souscal_div(mp) ;

    scalaire_regu.set_etat_qcq() ;
    scalaire_div.set_etat_qcq() ;
    dscal_div.set_etat_qcq() ;

    souscal_regu.std_base_scal() ;
    souscal_div.std_base_scal() ;

    source_scal().poisson_regular(k_div, nzet, unsgam1, *par,
				  scalaire.set(),
				  scalaire_regu.set(), scalaire_div.set(),
				  dscal_div, souscal_regu, souscal_div) ;
    
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


}
