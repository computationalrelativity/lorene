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


char map_af_lap_C[] = "$Header$" ;



/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:41  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.16  2000/08/16  10:35:41  eric
 * Suppression de Mtbl::dzpuis.
 *
 * Revision 2.15  2000/02/25  09:01:47  eric
 * Remplacement de (uu.get_dzpuis() == 0) par uu.check_dzpuis(0).
 *
 * Revision 2.14  2000/02/08  14:19:53  phil
 * correction annulation derniere zone
 *
 * Revision 2.13  2000/01/27  17:52:16  phil
 * corrections diverses et variees
 *
 * Revision 2.12  2000/01/26  13:10:34  eric
 * Reprototypage complet :
 * le resultat est desormais suppose alloue a l'exterieur de la routine
 * et est passe en argument (Cmp& resu).
 *
 * Revision 2.11  1999/11/30  12:49:53  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.10  1999/11/29  12:57:42  eric
 * REORGANISATION COMPLETE: nouveau prototype : Valeur --> Cmp
 *   Utilisation de la nouvelle arithmetique des Valeur's.
 *
 * Revision 2.9  1999/10/27  15:44:00  eric
 * Suppression du membre Cmp::c.
 *
 * Revision 2.8  1999/10/14  14:27:35  eric
 * Methode const.
 *
 * Revision 2.7  1999/09/06  16:26:03  phil
 * ajout de la version Cmp
 *
 * Revision 2.6  1999/09/06  14:51:20  phil
 * ajout du laplacien
 *
 * Revision 2.5  1999/05/03  15:22:00  phil
 * Correction des bases
 *
 * Revision 2.4  1999/04/28  10:33:02  phil
 * correction du cas zec_mult_r = 4
 *
 * Revision 2.3  1999/04/27  09:22:25  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/04/27  09:17:27  phil
 * corrections diverses et variees ....
 *
 * Revision 2.1  1999/04/26  17:24:21  phil
 * correction de gestion de base
 *
 * Revision 2.0  1999/04/26  16:33:55  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Fichiers include
// ----------------
#include <stdlib.h>
#include <math.h>

#include "map.h"
#include "cmp.h"
#include "type_parite.h"
#include "utilitaires.h"

//******************************************************************************


/*
 *  Calcul du laplacien d'un Cmp f dans le cas ou le mapping
 *  (coordonnees-grille) --> (coordonnees physiques) est affine. 
 *  Le Laplacien est calcule suivant l'equation:
 * 
 *    Lap(f) = d^2 f/dr^2 + 1/r ( 2 df/dr + 1/r Lap_ang(f) )		(1)
 * 
 *  avec 
 * 
 *    Lap_ang(f) := d^2 f/dtheta^2 + 1/tan(theta) df/dtheta 
 *			+ 1/sin^2(theta) d^2 f/dphi^2			(2)
 *
 *  Le laplacien angulaire (2) est calcule par passage aux harmoniques 
 *  spheriques, suivant la formule
 *
 *    Lap_ang( f_{lm} Y_l^m ) = - l(l+1) f_{lm} Y_l^m			(3)
 *   
 *  Dans la zone externe compactifiee (ZEC), la routine calcule soit Lap(f), 
 *   soit r^2 Lap(f), soit r^4 Lap(f) suivant la valeur du drapeau zec_mult_r :
 * 
 *   -- pour zec_mult_r = 0, on calcule Lap(f) suivant la formule
 *
 *        Lap(f) = u^2 [ u^2 d^2 f/du^2 + Lap_ang(f) ] ,  avec u = 1/r  (4)
 * 
 *   -- pour zec_mult_r = 2, on calcule r^2 Lap(f) suivant la formule
 *
 *        r^2 Lap(f) = u^2 d^2 f/du^2 + Lap_ang(f)			(5)
 * 
 *   -- pour zec_mult_r = 4, on calcule 4^2 Lap(f) suivant la formule
 *
 *	  r^4 Lap(f) = d^2 f /du^2 + 1/u^2 Lap_ang(f)			(6)
 *  
 *
 *
 * Entree:
 * ------
 *  const Cmp& uu :	    champ f dont on veut calculer le laplacien
 *			    
 *
 *  int zec_mult_r :	    drapeau indiquant la quantite calculee dans la ZEC:
 *			     zec_mult_r = 0  : lapff = Lap(f)
 *			     zec_mult_r = 2  : lapff = r^2 Lap(f)
 *			     zec_mult_r = 4  : lapff = r^4 Lap(f)
 * Sortie:
 * ------
 *   Cmp& resu :	 Lap(f)
 *
 */

void Map_af::laplacien(const Cmp& uu, int zec_mult_r, Cmp& resu) const {
 
    assert (uu.get_etat() != ETATNONDEF) ; 
    assert (uu.get_mp()->get_mg() == mg) ; 

    // Particular case of null value:

    if (uu.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
	return ; 
    }

    assert( uu.get_etat() == ETATQCQ ) ; 
    assert( uu.check_dzpuis(0) ) ; 
    resu.set_etat_qcq() ;    

    int nz = get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ;        

    // On recupere les bases d'entree : 
    Base_val base_entree =  (uu.va).base  ;
    
    // Separation zones internes / ZEC :
    Cmp uuext(uu) ; 
    
    Valeur ff = uu.va ;

    if (get_mg()->get_type_r(nzm1) == UNSURR) {  // il existe une ZEC
	uuext.annule(0, nzm1-1) ;
	ff.annule(nzm1) ;
    }
    else {
	uuext.set_etat_zero() ;		// pas de ZEC
    }


    //=========================================================================
    // Calcul dans les zones internes (noyau + coquilles)
    //=========================================================================

	//----------------------------
	// Derivations par rapport a x  
	//----------------------------

	Valeur d2ff = ff.d2sdx2() ;		    // d^2/dx^2  partout 
	
	Valeur dff = ff.dsdx() ;		    // d/dx partout

	//---------------------------
	// Calcul de 1/x Lap_ang(f) ---> resultat dans ff
	//---------------------------

	//... Passage en harmoniques spheriques
	ff.coef() ;     
	ff.ylm() ;		

	//... Multiplication par -l*(l+1)  
	ff = ff.lapang() ; 

	//... Division par x:    
	ff = ff.sx() ;		// 1/x	     ds. le noyau
				// Id	     ds. les coquilles
    
	//----------------------------------------------
	// On repasse dans l'espace des configurations 
	//  pour effectuer les multiplications par 
	//  les derivees du chgmt. de var. 
	//----------------------------------------------

	d2ff.coef_i() ;
	dff.coef_i() ;
	ff.coef_i() ;
	
	
	//-----------------------------------------------
	//  Calcul de 1/x ( 2 df/dr + 1/r Lap_ang(f) ) 
	//  Le resultat est mis dans ff
	//-----------------------------------------------

	Base_val sauve_base = dff.base ;
	ff = double(2) * ( dxdr * dff)  + xsr * ff ; 
	ff.base = sauve_base ;
				
	ff = ff.sx() ;	    // 1/x	     ds. le noyau
			    // Id	     ds. les coquilles
	ff.coef_i() ;

	//---------------------------------------------
	// Calcul de Lap(f) suivant l'Eq. (1)
	//  Le resultat est mis dans ff
	//-----------------------------------------------
	
	sauve_base = d2ff.base ;
	ff = dxdr * dxdr * d2ff + xsr * ff ;
	ff.base = sauve_base ;
	

    if (get_mg()->get_type_r(nzm1) == UNSURR) {  // il existe une ZEC

    //=========================================================================
    // Calcul dans la ZEC
    //=========================================================================
 
	    Valeur& ffext = uuext.va ;
 
	    //----------------------------
	    // Derivation par rapport a x  
	    //----------------------------
	    
	    d2ff = ffext.d2sdx2() ;	    // d^2/dx^2  partout 
	    
	    //---------------------------
	    // Calcul de Lap_ang(f) ---> resultat dans ffext
	    //---------------------------

	    //... Passage en harmoniques spheriques    
	    ffext.coef() ; 
	    ffext.ylm() ;		

	    //... Multiplication par -l*(l+1)  
	    
	    ffext = ffext.lapang() ; 
	    
	    switch (zec_mult_r) {
		
		case 0 : {	// calcul de Lap(f) suivant l'Eq. (4)
		    
		    d2ff.mult2_xm1_zec() ;	// multiplication de d^2f/dx^2 
						//  par (x-1)^2
		    
		    sauve_base = ffext.base ;
		    ffext = dxdr * dxdr / (xsr*xsr) * d2ff + ffext ; 
		    ffext.base = sauve_base ;
		    
		    ffext.mult2_xm1_zec() ;	// multiplication par (x-1)^2
		    ffext.coef_i() ; 
		    sauve_base = ffext.base ;
		    ffext = ffext / (xsr*xsr) ; 		    
		    ffext.base = sauve_base ;
		    break ; 
		}
		
		case 2 : {	// calcul de r^2 Lap(f) suivant l'Eq. (5)
		    
		    d2ff.mult2_xm1_zec() ;	// multiplication de d^2f/dx^2 
						//  par (x-1)^2		    
		    sauve_base = ffext.base ;
		    ffext = dxdr*dxdr / (xsr*xsr) * d2ff + ffext ;
		    ffext.base = sauve_base ; 
		    break ; 
		}
		
		case 4 : {	// calcul de r^4 Lap(f) suivant l'Eq. (6)
		    
		    ffext = ffext.sx2() ; // division de Lap_ang(f) par (x-1)^2

		    sauve_base = ffext.base ;
		    ffext = dxdr*dxdr * d2ff + xsr*xsr * ffext ;
		    ffext.base = sauve_base ;
		    break ; 
		}
		
		default : {
		    cout << "Map_af::laplacien : unknown value of zec_mult_r !"
			 << endl << " zec_mult_r = " << zec_mult_r << endl ; 
		    abort() ; 
		}
	    }	// fin des differents cas pour zec_mult_r

	// Resultat final

	ff = ff + ffext ; 

    } // fin du cas ou il existe une ZEC 
    
    // Les bases de sorties sont egales aux bases d'entree:
    ff.base = base_entree ;
    resu = ff ;
    resu.set_dzpuis(zec_mult_r) ; 
    
}




