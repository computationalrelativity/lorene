/*
 *  Member functions of the class Map_radial for various r manipulations
 *  of the Cmp's.
 */

/*
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Philippe Grandclement
 *   Copyright (c) 2001 Jerome Novak
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


char map_radial_r_manip_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2004/01/27 09:33:48  j_novak
 * New method Map_radial::div_r_zec
 *
 * Revision 1.6  2003/11/04 23:00:16  e_gourgoulhon
 * Method div_tant is now defined in file map_radial_th_manip.C.
 *
 * Revision 1.5  2003/10/27 09:02:19  j_novak
 * Corrected a bug in the case of null CED
 *
 * Revision 1.4  2003/10/17 15:07:29  j_novak
 * The order of operations in div_tant() has been changed.
 *
 * Revision 1.3  2003/10/15 10:41:49  e_gourgoulhon
 * Added new method div_tant.
 *
 * Revision 1.2  2002/08/13 08:02:45  j_novak
 * Handling of spherical vector/tensor components added in the classes
 * Mg3d and Tenseur. Minor corrections for the class Metconf.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.17  2001/10/29  15:34:35  novak
 * Ajout de Map_radial::div_r
 *
 * Revision 1.16  2000/08/31 14:19:26  eric
 * *** empty log message ***
 *
 * Revision 1.15  2000/08/31  13:04:05  eric
 * Ajout de la fonction mult_rsint.
 * Reecriture de div_rsint.
 *
 * Revision 1.14  2000/08/31  11:45:19  eric
 * Suppression de assert (!((ci_ext.check_dzpuis(0)))) dans
 *   Map_radial::mult_r
 *
 * Revision 1.13  2000/07/20  14:21:27  eric
 * Ajout de la fonction div_rsint.
 *
 * Revision 1.12  2000/05/22  15:12:27  phil
 * *** empty log message ***
 *
 * Revision 1.11  2000/05/22  14:55:47  phil
 * la fonction mult_r se contente de changer le dzpuis sans faire d'operation !
 * /.
 *
 * Revision 1.10  2000/05/22  14:39:28  phil
 * ajout de inc_dzpuis et dec_dzpuis
 *
 * Revision 1.9  2000/03/31  13:25:12  eric
 * *** empty log message ***
 *
 * Revision 1.8  2000/03/31  13:11:08  eric
 * Map_radial::mult_r : traitement ameliore dans la ZEC.
 *
 * Revision 1.7  1999/12/02  11:27:34  eric
 * *** empty log message ***
 *
 * Revision 1.6  1999/12/02  11:04:05  eric
 * Reorganisation complete de la routine mult_r.
 * Appel de Valeur::mult_x().
 *
 * Revision 1.5  1999/11/30  16:26:26  eric
 * *** empty log message ***
 *
 * Revision 1.4  1999/11/30  15:48:20  eric
 * *** empty log message ***
 *
 * Revision 1.3  1999/11/30  15:31:06  eric
 * *** empty log message ***
 *
 * Revision 1.2  1999/11/30  15:15:02  eric
 * *** empty log message ***
 *
 * Revision 1.1  1999/11/30  14:22:50  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

#include "map.h"
#include "cmp.h"
#include "tensor.h"


			//---------------------------//
			//	    mult_r	     //
			//---------------------------//

void Map_radial::mult_r(Cmp& ci) const {
    
    // Verifications d'usage :
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Nothing to do if the Cmp is null :
    if (ci.get_etat() == ETATZERO) {
	return ; 
    }

    assert(ci.get_etat() == ETATQCQ) ;
        
    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;


    if (mg->get_type_r(nzm1) == UNSURR) {   // Case with ZEC
					    // -------------
    
	// Decomposition inner domains / external domain : 
	// ----------------------------------------------
    
	Cmp ci_ext = ci ; 
	ci_ext.annule(0, nzm1-1) ; 

	ci.annule(nzm1) ; 

	// Inner domains: multiplication by r :
	// ----------------------------------
	//Confort :
	Valeur& val = ci.va ; 
	assert(val.get_mg() == mg) ; 

	val = val.mult_x() ; // Multiplication by xi in the nucleus
		             // Identity in the shells

	Base_val sauve_base = val.base ; 
	val = val / xsr ;    // R/xi in the nucleus
			     // R in the shells
	val.base = sauve_base ; 
	 
	// External domain
	// ---------------
	
	// On change juste le dzpuis !

	//## assert (!((ci_ext.check_dzpuis(0)))) ;
	// On fait just dec_dzpuis ...
	ci_ext.set_dzpuis (ci.get_dzpuis()-1) ;
    
	// Recombination
	// -------------

	ci = ci + ci_ext ; 
	    
    }
    else{   // Case without ZEC
	    //-----------------
	Valeur& uu = ci.va ; 
	uu = uu.mult_x() ;  // Multiplication by xi in the nucleus
			    // Identity in the shells
	
	Base_val sauve_base = uu.base ; 
	uu = uu / xsr ;	    // R/xi in the nucleus
			    // R in the shells 
	uu.base = sauve_base ; 
    }
    
    
}



			//---------------------------//
			//	    mult_r_zec	     //
			//---------------------------//

void Map_radial::mult_r_zec(Cmp& ci) const {
    
    // Verifications d'usage :
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Nothing to do if the Cmp is null :
    if (ci.get_etat() == ETATZERO) {
	return ; 
    }

    assert(ci.get_etat() == ETATQCQ) ;
    
    //Confort :
    Valeur& uu = ci.va ; 
    assert(uu.get_mg() == mg) ; 
    
    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;


    if (mg->get_type_r(nzm1) == UNSURR) {   // Case with ZEC
					    // -------------
    
	// On stocke tout sauf la ZEC

	Valeur val = uu ;
	val.annule(nzm1) ;
   
	// On ne travaile que sur la ZEC :

	Valeur val_ext = uu ;
	val_ext.annule(0, nzm1-1) ; 
	
	val_ext.sxm1_zec() ; // division par (x-1) dans l'espace des coefs.

	Base_val sauve_base = val_ext.base ; 
	val_ext = xsr * val_ext ;
	val_ext.base = sauve_base ; 
    
	// Et on reconstruit le resultat ...
	uu = val + val_ext ; 
	    
    }
    else{   // Case without ZEC
	    //-----------------

	return ; 

    }
    
    
}

			//---------------------------//
			//	    mult_rsint	     //
			//---------------------------//

void Map_radial::mult_rsint(Cmp& ci) const {
    
    assert(ci.get_etat() != ETATNONDEF) ;
    
    if (ci.get_etat() == ETATZERO) {
	return ;			 // Nothing to do if the Cmp is null 
    }

    assert(ci.get_etat() == ETATQCQ) ;
            
    Valeur& val = ci.va ; 
    assert(val.get_mg() == mg) ; 

    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;

    // 1/ Multiplication by sin(theta)
    // -------------------------------
    
    val = val.mult_st() ;	// Multiplication by sin(theta)

    // 2/ Multiplication by r
    // ----------------------    

    Cmp ci_ext(this) ; 

    if (mg->get_type_r(nzm1) == UNSURR) {   // Case with ZEC
					    // -------------
    	// Decomposition inner domains / external domain : 
	// ----------------------------------------------
    
	ci_ext = ci ; 
	ci_ext.annule(0, nzm1-1) ; 
	ci.annule(nzm1) ; 

	// External domain 
	// ---------------
	Valeur& val_ext = ci_ext.va ; 
	assert(val_ext.get_mg() == mg) ; 

	val_ext.sxm1_zec() ;		// Division by (xi-1) 
		
	Base_val sauve_base = val_ext.base ; 
	val_ext = val_ext * xsr ;	// Multiplication by R(xi-1)
	val_ext.base = sauve_base ; 
	
    }
    else{   // Case without ZEC
	    //-----------------

	ci_ext = 0 ; 

    }

    // Inner domains: 
    // -------------
    val = val.mult_x() ;  // Multiplication by xi in the nucleus
		          // Identity in the shells

    Base_val sauve_base = val.base ; 
    val = val / xsr ;	  // Multiplication by R/xi in the nucleus
			  //			  R in the shells
    val.base = sauve_base ; 
	 
    // Recombination
    // -------------

    ci = ci + ci_ext ;     
    
}



			//---------------------------//
			//	    div_rsint	     //
			//---------------------------//

void Map_radial::div_rsint(Cmp& ci) const {
    
    assert(ci.get_etat() != ETATNONDEF) ;
    
    if (ci.get_etat() == ETATZERO) {
	return ;			 // Nothing to do if the Cmp is null 
    }

    assert(ci.get_etat() == ETATQCQ) ;
            
    Valeur& val = ci.va ; 
    assert(val.get_mg() == mg) ; 

    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;

    // 1/ Division by sin(theta)
    // -------------------------
    
    val = val.ssint() ;		// Division by sin(theta)


    // 2/ Division by r
    // ----------------    

    Cmp ci_ext(this) ; 

    if (mg->get_type_r(nzm1) == UNSURR) {   // Case with ZEC
					    // -------------
    
	// Decomposition inner domains / external domain : 
	// ----------------------------------------------
    
	ci_ext = ci ; 
	ci_ext.annule(0, nzm1-1) ; 
	ci.annule(nzm1) ; 

	// External domain
	// ---------------
	Valeur& val_ext = ci_ext.va ; 
	assert(val_ext.get_mg() == mg) ; 

	val_ext = val_ext.mult_x() ;	// Multiplication by (xi-1)
	
	Base_val sauve_base = val_ext.base ; 
	val_ext = val_ext / xsr ;	// Division by (xi-1)/R
	val_ext.base = sauve_base ; 

    }
    else{   // Case without ZEC
	    //-----------------

	ci_ext = 0 ; 

    }

    // Inner domains: 
    // -------------

    val = val.sx() ;	// Division by xi in the nucleus
			// Identity in the shells

    Base_val sauve_base = val.base ; 
    val = val * xsr ;    // Multiplication by xi/R in the nucleus
			 // Multiplication by 1/R  in the shells
    val.base = sauve_base ; 

    // Recombination
    // -------------

    ci = ci + ci_ext ;         
    
}

			//---------------------------//
			//	    div_r	     //
			//---------------------------//

void Map_radial::div_r(Cmp& ci) const {
    
    assert(ci.get_etat() != ETATNONDEF) ;
    
    if (ci.get_etat() == ETATZERO) {
	return ;			 // Nothing to do if the Cmp is null 
    }

    assert(ci.get_etat() == ETATQCQ) ;
            
    Valeur& val = ci.va ; 
    assert(val.get_mg() == mg) ; 

    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;

    Cmp ci_ext(this) ; 

    if (mg->get_type_r(nzm1) == UNSURR) {   // Case with ZEC
					    // -------------
    
	// Decomposition inner domains / external domain : 
	// ----------------------------------------------
    
	ci_ext = ci ; 
	ci_ext.annule(0, nzm1-1) ; 
	ci.annule(nzm1) ; 

	// External domain
	// ---------------
	Valeur& val_ext = ci_ext.va ; 
	assert(val_ext.get_mg() == mg) ; 

	val_ext = val_ext.mult_x() ;	// Multiplication by (xi-1)
	
	Base_val sauve_base = val_ext.base ; 
	val_ext = val_ext / xsr ;	// Division by (xi-1)/R
	val_ext.base = sauve_base ; 

    }
    else{   // Case without ZEC
	    //-----------------

	ci_ext = 0 ; 

    }

    // Inner domains: 
    // -------------

    val = val.sx() ;	// Division by xi in the nucleus
			// Identity in the shells

    Base_val sauve_base = val.base ; 
    val = val * xsr ;    // Multiplication by xi/R in the nucleus
			 // Multiplication by 1/R  in the shells
    val.base = sauve_base ; 

    // Recombination
    // -------------

    ci = ci + ci_ext ;         
    
}


			//---------------------------//
			//	    div_r_zec	     //
			//---------------------------//

void Map_radial::div_r_zec(Scalar& uu) const {
    
    // Verifications d'usage :
    assert(uu.get_etat() != ETATNONDEF) ;
    
    // Nothing to do if the Cmp is null :
    if (uu.get_etat() == ETATZERO) {
	return ; 
    }

    assert(uu.get_etat() == ETATQCQ) ;
    
    //Confort :
    const Valeur& vu = uu.get_spectral_va() ; 
    assert(vu.get_mg() == mg) ; 
    
    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;

    if (mg->get_type_r(nzm1) == UNSURR) {   // Case with ZEC
					    // -------------
    	// On stocke tout sauf la ZEC

	Valeur val = vu ;
	val.annule(nzm1) ;
   
	// On ne travaile que sur la ZEC :

	Valeur val_ext = vu ;
	val_ext.annule(0, nzm1-1) ; 
	
	val_ext.mult_xm1_zec() ; // division par (x-1) dans l'espace des coefs.

	Base_val sauve_base = val_ext.base ; 
	val_ext = val_ext / xsr ;
	val_ext.base = sauve_base ; 
    
	// Et on reconstruit le resultat ...
	uu.set_spectral_va() = val + val_ext ; 
	    
    }
    else{   // Case without ZEC
	    //-----------------

	return ; 

    }
    
    
}


			//---------------------------//
			//	  dec_dzpuis	     //
			//---------------------------//
			
void Map_radial::dec_dzpuis(Cmp& ci) const {

    // Verifications d'usage :
    assert(ci.get_etat() != ETATNONDEF) ;
    
    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;

    // Nothing to do if the Cmp is null or if there is no ZEC:
    if ((ci.get_etat() == ETATZERO) || (mg->get_type_r(nzm1) != UNSURR)) {
      ci.set_dzpuis( ci.get_dzpuis() - 1 ) ; 
      return ; 
    }

    assert(ci.get_etat() == ETATQCQ) ;
    
    Valeur& uu = ci.va ; 
    assert(uu.get_mg() == mg) ; 
    
    
    // Decomposition inner domains / external domain : 
    // ----------------------------------------------
    Valeur uu_ext = uu ; 
    uu_ext.annule(0, nzm1-1) ; 

    uu.annule(nzm1) ; 
    
    // Computation in the external domain (division by r)
    // ----------------------------------
    
    uu_ext.mult_xm1_zec() ;    // Multiplication by (xi-1) in the ZEC

    Base_val sauve_base = uu_ext.base ; 
    uu_ext = uu_ext / xsr ;	    // u^2/(xi-1) in the ZEC
    uu_ext.base = sauve_base ; 

    // Final result:
    // ------------
    uu = uu + uu_ext ; 
 
    ci.set_dzpuis( ci.get_dzpuis() - 1 ) ; 
    
}

			//---------------------------//
			//	  inc_dzpuis	     //
			//---------------------------//
			
void Map_radial::inc_dzpuis(Cmp& ci) const {

    // Verifications d'usage :
    assert(ci.get_etat() != ETATNONDEF) ;
    
    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;

    // Nothing to do if the Cmp is null or if there is no ZEC:
    if ((ci.get_etat() == ETATZERO) || (mg->get_type_r(nzm1) != UNSURR)) {
      ci.set_dzpuis( ci.get_dzpuis() + 1 ) ; 
      return ; 
    }

    assert(ci.get_etat() == ETATQCQ) ;
    
    Valeur& uu = ci.va ; 
    assert(uu.get_mg() == mg) ; 
    
    
    // Decomposition inner domains / external domain : 
    // ----------------------------------------------
    Valeur uu_ext = uu ; 
    uu_ext.annule(0, nzm1-1) ;
     
    uu.annule(nzm1) ; 
    
    // Computation in the external domain (multiplication by r)
    // ----------------------------------
    
    uu_ext.sxm1_zec() ;    // Division by (xi-1) in the ZEC

    Base_val sauve_base = uu_ext.base ; 
    uu_ext = uu_ext * xsr ;	    // (xi-1)/u in the ZEC
    uu_ext.base = sauve_base ; 

    // Final result:
    // ------------
    uu = uu + uu_ext ; 
 
    ci.set_dzpuis( ci.get_dzpuis() + 1 ) ; 

}


			//---------------------------//
			//	  dec2_dzpuis	     //
			//---------------------------//
			
void Map_radial::dec2_dzpuis(Cmp& ci) const {

    // Verifications d'usage :
    assert(ci.get_etat() != ETATNONDEF) ;
    
    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;

    // Nothing to do if the Cmp is null or if there is no ZEC:
    if ((ci.get_etat() == ETATZERO) || (mg->get_type_r(nzm1) != UNSURR)) {
      ci.set_dzpuis( ci.get_dzpuis() - 2 ) ; 
      return ; 
    }

    assert(ci.get_etat() == ETATQCQ) ;
    
    Valeur& uu = ci.va ; 
    assert(uu.get_mg() == mg) ; 
    
    
    // Decomposition inner domains / external domain : 
    // ----------------------------------------------
    Valeur uu_ext = uu ; 
    uu_ext.annule(0, nzm1-1) ; 

    uu.annule(nzm1) ; 
    
    // Computation in the external domain (division by r^2)
    // ----------------------------------
    
    uu_ext.mult2_xm1_zec() ;    // Multiplication by (xi-1)^2 in the ZEC

    Base_val sauve_base = uu_ext.base ; 
    uu_ext = uu_ext / (xsr*xsr) ;	    // u^2/(xi-1)^2 in the ZEC
    uu_ext.base = sauve_base ; 

    // Final result:
    // ------------
    uu = uu + uu_ext ; 
 
    ci.set_dzpuis( ci.get_dzpuis() - 2 ) ; 
    
}

			//---------------------------//
			//	  inc2_dzpuis	     //
			//---------------------------//
			
void Map_radial::inc2_dzpuis(Cmp& ci) const {

    // Verifications d'usage :
    assert(ci.get_etat() != ETATNONDEF) ;
    
    int nz = mg->get_nzone() ;
    int nzm1 = nz-1 ;

    // Nothing to do if the Cmp is null or if there is no ZEC:
    if ((ci.get_etat() == ETATZERO) || (mg->get_type_r(nzm1) != UNSURR)) {
      ci.set_dzpuis( ci.get_dzpuis() + 2 ) ; 
      return ; 
    }

    assert(ci.get_etat() == ETATQCQ) ;
    
    Valeur& uu = ci.va ; 
    assert(uu.get_mg() == mg) ; 
    
    
    // Decomposition inner domains / external domain : 
    // ----------------------------------------------
    Valeur uu_ext = uu ; 
    uu_ext.annule(0, nzm1-1) ;
     
    uu.annule(nzm1) ; 
    
    // Computation in the external domain (multiplication by r^2)
    // ----------------------------------
    
    uu_ext.sxm1_zec() ;    // Division by (xi-1) in the ZEC
    uu_ext.sxm1_zec() ;    // Division by (xi-1) in the ZEC

    Base_val sauve_base = uu_ext.base ; 
    uu_ext = uu_ext * (xsr*xsr) ;	    // (xi-1)^2/u^2 in the ZEC
    uu_ext.base = sauve_base ; 

    // Final result:
    // ------------
    uu = uu + uu_ext ; 
 
    ci.set_dzpuis( ci.get_dzpuis() + 2 ) ; 

}

