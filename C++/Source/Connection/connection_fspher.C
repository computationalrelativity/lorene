/*
 *  Methods of class Connection_fspher.
 *
 *	(see file connection.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
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

char connection_fspher_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.16  2004/01/04 21:00:50  e_gourgoulhon
 * Better handling of tensor symmetries in methods p_derive_cov() and
 * p_divergence() (thanks to the new class Tensor_sym).
 *
 * Revision 1.15  2004/01/01 11:24:04  e_gourgoulhon
 * Full reorganization of method p_derive_cov: the main loop is now
 * on the indices of the *output* tensor (to take into account
 * symmetries in the input and output tensors).
 *
 * Revision 1.14  2003/12/27 14:59:52  e_gourgoulhon
 * -- Method derive_cov() suppressed.
 * -- Change of the position of the derivation index from the first one
 *    to the last one in methods p_derive_cov() and p_divergence().
 *
 * Revision 1.13  2003/11/03 13:37:58  j_novak
 * Still dzpuis...
 *
 * Revision 1.12  2003/11/03 11:14:18  j_novak
 * Treatment of the case dzpuis = 4.
 *
 * Revision 1.11  2003/11/03 10:58:30  j_novak
 * Treatment of the general case for divergence.
 *
 * Revision 1.10  2003/10/22 13:08:03  j_novak
 * Better handling of dzpuis flags
 *
 * Revision 1.9  2003/10/16 15:26:48  e_gourgoulhon
 * Name of method Scalar::div_r_ced() changed to Scalar::div_r_inc2().
 *
 * Revision 1.8  2003/10/16 14:21:36  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.7  2003/10/15 10:46:18  e_gourgoulhon
 * Introduced call to the new method Scalar::div_tant to perform
 * division by tan(theta) in derive_cov.
 *
 * Revision 1.6  2003/10/11 16:45:43  e_gourgoulhon
 * Suppressed the call to Itbl::set_etat_qcq() after
 * the construction of the Itbl's.
 *
 * Revision 1.5  2003/10/11 14:39:50  e_gourgoulhon
 * Suppressed declaration of unusued arguments in some methods.
 *
 * Revision 1.4  2003/10/06 13:58:47  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.3  2003/10/05 21:09:23  e_gourgoulhon
 * Method derive_cov: multiplication by r^2 in the CED.
 *
 * Revision 1.2  2003/10/01 21:49:45  e_gourgoulhon
 * First version of derive_cov --- not tested yet.
 *
 * Revision 1.1  2003/10/01 15:42:49  e_gourgoulhon
 * still ongoing...
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>

// Lorene headers
#include "connection.h"

        //------------------------------//
        //          Constructors        //
        //------------------------------//

// Contructor from a spherical flat-metric-orthonormal basis

Connection_fspher::Connection_fspher(const Map& mpi, const Base_vect_spher& bi) 
  : Connection_flat(mpi, bi) {

}		

// Copy constructor
Connection_fspher::Connection_fspher(const Connection_fspher& ci) 
  : Connection_flat(ci) {

}		

	
        //----------------------------//
        //          Destructor        //
        //----------------------------//

Connection_fspher::~Connection_fspher(){
	
}


        //--------------------------------//
        //      Mutators / assignment     //
        //--------------------------------//


void Connection_fspher::operator=(const Connection_fspher& ) {
	
  cout << "Connection_fspher::operator= : not implemented yet !" << endl ; 
  abort() ; 

}	


        //-----------------------------//
        //    Computational methods    //
        //-----------------------------//

// Covariant derivative, returning a pointer.
//-------------------------------------------

Tensor* Connection_fspher::p_derive_cov(const Tensor& uu) const {

    // Notations: suffix 0 in name <=> input tensor
    //            suffix 1 in name <=> output tensor

    int valence0 = uu.get_valence() ; 
    int valence1 = valence0 + 1 ; 
    int valence1m1 = valence1 - 1 ; // same as valence0, but introduced for 
                                    // the sake of clarity
	
    // Protections
    // -----------
    if (valence0 >= 1) {
        assert(uu.get_triad() == triad) ; 
    }

    // Creation of the result (pointer)
    // --------------------------------
    Tensor* resu ;

    // If uu is a Scalar, the result is a vector
    if (valence0 == 0) 
        resu = new Vector(*mp, COV, triad) ;
    else {

        // Type of indices of the result :
        Itbl tipe(valence1) ; 
        const Itbl& tipeuu = uu.get_index_type() ;  
        for (int id = 0; id<valence0; id++) {
            tipe.set(id) = tipeuu(id) ;   // First indices = same as uu
        }
        tipe.set(valence1m1) = COV ;  // last index is the derivation index

        // if uu is a Tensor_sym, the result is also a Tensor_sym:
        const Tensor* puu = &uu ; 
        const Tensor_sym* puus = dynamic_cast<const Tensor_sym*>(puu) ; 
        if ( puus != 0x0 ) {    // the input tensor is symmetric
            resu = new Tensor_sym(*mp, valence1, tipe, *triad,
                                  puus->sym_index1(), puus->sym_index2()) ;
        }
        else {  
            resu = new Tensor(*mp, valence1, tipe, *triad) ;  // no symmetry  
        }

    }

    int ncomp1 = resu->get_n_comp() ; 
	
    Itbl ind1(valence1) ; // working Itbl to store the indices of resu
    Itbl ind0(valence0) ; // working Itbl to store the indices of uu
    Itbl ind(valence0) ;  // working Itbl to store the indices of uu
	
    Scalar tmp(*mp) ;	// working scalar

    // Loop on all the components of the output tensor
    // -----------------------------------------------
    for (int ic=0; ic<ncomp1; ic++) {
    
        // indices corresponding to the component no. ic in the output tensor
        ind1 = resu->indices(ic) ; 
    
        // Component no. ic:
        Scalar& cresu = resu->set(ind1) ; 
		
        // Indices of the input tensor
        for (int id = 0; id < valence0; id++) {
            ind0.set(id) = ind1(id) ; 
        }
 
        // dzpuis
        int dz_resu = uu(ind0).get_dzpuis() ;
        bool dz4 = (dz_resu == 4) ;
        dz_resu += dz4 ? 0 : 2 ;

        // Value of last index (derivation index)
        int k = ind1(valence1m1) ; 
        
        switch (k) {
        
        case 1 : {  // Derivation index = r
                    //---------------------
	
            cresu = (uu(ind0)).dsdr(dz_resu) ; 	// d/dr
		
            // all the connection coefficients Gamma^i_{jk} are zero for k=1
            break ; 
        }

        case 2 : {  // Derivation index = theta
                    //-------------------------
			
            cresu = (uu(ind0)).srdsdt(dz_resu) ;  // 1/r d/dtheta 
		
            // Loop on all the indices of uu
            for (int id=0; id<valence0; id++) {
		
                switch ( ind0(id) ) {
				
                case 1 : {	// Gamma^r_{l theta} V^l 
	                        // or -Gamma^l_{r theta} V_l 
	            ind = ind0 ; 
	            ind.set(id) = 2 ;   // l = theta

	            // Division by r in all domains but the CED (where a
	            //  multiplication by r is performed instead)
	            tmp = uu(ind) ; 
	            dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	            cresu -= tmp ; 
	            break ; 
                }
				
                case 2 : {	// Gamma^theta_{l theta} V^l 
	                        // or -Gamma^l_{theta theta} V_l
	            ind = ind0 ; 
	            ind.set(id) = 1 ;   // l = r
	            tmp = uu(ind) ; 
	            dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	            cresu += tmp ; 
	            break ; 
                }
				
                case 3 : {	// Gamma^phi_{l theta} V^l 
	                        // or -Gamma^l_{phi theta} V_l
	        break ; 
                }
				
                default : {
	            cerr << "Connection_fspher::p_derive_cov : index problem ! "
	           << endl ; 
	            abort() ;  
                }
                }

            }
            break ; 
        }


        case 3 : {  // Derivation index = phi
                    //-----------------------
					
            cresu = (uu(ind0)).srstdsdp(dz_resu) ;  // 1/(r sin(theta)) d/dphi 	
		
            // Loop on all the indices of uu
            for (int id=0; id<valence0; id++) {
		
                switch ( ind0(id) ) {
				
                case 1 : {	// Gamma^r_{l phi} V^l 
	                        // or -Gamma^l_{r phi} V_l 
	            ind = ind0 ; 
	            ind.set(id) = 3 ;   // l = phi
	            tmp = uu(ind) ; 
	            dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	            cresu -= tmp ; 
	            break ; 
                }
				
                case 2 : {	// Gamma^theta_{l phi} V^l 
	                        // or -Gamma^l_{theta phi} V_l
	            ind = ind0 ; 
	            ind.set(id) = 3 ;   // l = phi
	            tmp = uu(ind) ; 
	            dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	            tmp.div_tant() ; 	// division by tan(theta)
					
	            cresu -= tmp ; 
	            break ; 
                }
				
                case 3 : {	// Gamma^phi_{l phi} V^l 
	                        // or -Gamma^l_{phi phi} V_l
						
	            ind = ind0 ; 
	            ind.set(id) = 1 ;   // l = r
	            tmp = uu(ind) ; 
	            dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	            cresu += tmp ; 

	            ind.set(id) = 2 ;   // l = theta
	            tmp = uu(ind) ; 
	            dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	            tmp.div_tant() ; 	// division by tan(theta)

	            cresu += tmp ; 
	            break ; 
                }
				
                default : {
	            cerr << "Connection_fspher::p_derive_cov : index problem ! "
	            << endl ; 
	            abort() ;  
                }
                }

            }
            
            break ; 
        }

        default : {
	    cerr << "Connection_fspher::p_derive_cov : index problem ! \n" ;
	    abort() ;  
        }

        } // End of switch on the derivation index


    } // End of loop on all the components of the output tensor

    // C'est fini !
    // -----------
    return resu ; 

}



// Divergence, returning a pointer.
//---------------------------------

Tensor* Connection_fspher::p_divergence(const Tensor& uu) const {

    // Notations: suffix 0 in name <=> input tensor
    //            suffix 1 in name <=> output tensor

    int valence0 = uu.get_valence() ; 
    int valence1 = valence0 - 1 ; 
    int valence0m1 = valence0 - 1 ; // same as valence1 but introduced for 
                                    // the sake of clarity
	
    // Protections
    // -----------
    assert (valence0 >= 1) ;
    assert (uu.get_triad() == triad) ; 

    // Last index must be contravariant:
    assert (uu.get_index_type(valence0-1) == CON) ;


    // Creation of the pointer on the result tensor
    // --------------------------------------------
    Tensor* resu ;

    if (valence0 == 1)      // if u is a Vector, the result is a Scalar
        resu = new Scalar(*mp) ;
    else {
    
        // Type of indices of the result :
        Itbl tipe(valence1) ; 
        const Itbl& tipeuu = uu.get_index_type() ;  
        for (int id = 0; id<valence1; id++) {
            tipe.set(id) = tipeuu(id) ;     // type of remaining indices = 
        }                                   //  same as uu indices

        if (valence0 == 2) {  // if u is a rank 2 tensor, the result is a Vector
            resu = new Vector(*mp, tipe(0), *triad) ;
        }
        else {
            // if uu is a Tensor_sym, the result might be also a Tensor_sym:
            const Tensor* puu = &uu ; 
            const Tensor_sym* puus = dynamic_cast<const Tensor_sym*>(puu) ; 
            if ( puus != 0x0 ) {    // the input tensor is symmetric

                if (puus->sym_index2() != valence0 - 1) {
                 
                    // the symmetry is preserved: 

                    if (valence1 == 2) {
                        resu = new Sym_tensor(*mp, tipe, *triad) ;
                    }
                    else {
                        resu = new Tensor_sym(*mp, valence1, tipe, *triad,
                                  puus->sym_index1(), puus->sym_index2()) ;
                    }
                }
                else { // the symmetry is lost: 
                
	            resu = new Tensor(*mp, valence1, tipe, *triad) ;
                }
            }
            else { // no symmetry in the input tensor:
	        resu = new Tensor(*mp, valence1, tipe, *triad) ;
            }
        }
 
    }

    int ncomp1 = resu->get_n_comp() ;
	
    Itbl ind0(valence0) ; // working Itbl to store the indices of uu
    Itbl ind1(valence1) ; // working Itbl to store the indices of resu
    Itbl ind(valence0) ;  // working Itbl to store the indices of uu
	
    Scalar tmp1(*mp) ;	// working scalar
    Scalar tmp2(*mp) ;	// working scalar

    // Loop on all the components of the output tensor
    for (int ic=0; ic<ncomp1; ic++) {
	
        ind1 = resu->indices(ic) ; 
        Scalar& cresu = resu->set(ind1) ;

        // Derivation index = r
        // --------------------
        int k = 1 ; 	

        // indices (ind1,k) in the input tensor
        for (int id = 0; id < valence1; id++) {
	    ind0.set(id) = ind1(id) ; 
        }
        ind0.set(valence0m1) = k ; 

        // dzpuis of the result
        //---------------------
        int dz_resu = uu(ind0).get_dzpuis() ;
        bool dz4 = (dz_resu == 4) ;
        dz_resu += dz4 ? 0 : 2 ;

        cresu = uu(ind0).dsdr(dz_resu) ; //dT^{l r}/dr

        // Derivation index = theta
        // ------------------------
        k = 2 ; 	

        // indices (ind1,k) in the input tensor
        for (int id = 0; id < valence1; id++) {
	    ind0.set(id) = ind1(id) ; 
        }
        ind0.set(valence0m1) = k ; 
		
        tmp1 = uu(ind0).dsdt() ; //dT^{l theta} /dtheta

        ind = ind0 ;
        ind.set(valence0m1) = 1 ;
        tmp1 += uu(ind) ;//Gamma^theta_{r theta}T^{l r} (div_r is done later)
    

        // Loop on all the (valence0-1) first indices of uu
        for (int id=0; id<valence0m1; id++) {
		
            switch ( ind0(id) ) {
            case 1 : {	// Gamma^r_{l theta} V^l 
	                // or -Gamma^l_{r theta} V_l 
	        ind = ind0 ; 
	        ind.set(id) = 2 ;   // l = theta
	        tmp1 -= uu(ind) ; 
	        break ; 
            }
				
            case 2 : {	// Gamma^theta_{l theta} V^l 
	                // or -Gamma^l_{theta theta} V_l
	        ind = ind0 ; 
	        ind.set(id) = 1 ;   // l = r
	        tmp1 += uu(ind) ; 
	        break ; 
            }
				
            case 3 : {	// Gamma^phi_{l theta} V^l 
	                // or -Gamma^l_{phi theta} V_l
	        break ; 
            }
				
            default : {
	        cout << "Connection_fspher::p_divergence : index problem ! "
	            << endl ; 
	        abort() ;  
            }
            }
            
        }

        // Derivation index = phi
        // ----------------------
        k = 3 ; 			

        // indices (ind1,k) in the input tensor
        for (int id = 0; id < valence1; id++) {
            ind0.set(id) = ind1(id) ; 
        }
        ind0.set(valence0m1) = k ; 
    
        tmp1 += uu(ind0).stdsdp() ; // 1/sin(theta) dT^phi / dphi
    
        ind = ind0 ;
        ind.set(valence0m1) = 1 ;
        tmp1 += uu(ind) ;//Gamma^phi_{r phi} T^{l r} (div_r is done later)
        ind.set(valence0m1) = 2 ;
        tmp2 = uu(ind) ;//Gamma^phi_{theta phi} T^{l theta} (div_r is done later)

        // Loop on all the (valence0-1) first indices of uu
        for (int id=0; id<valence0-1; id++) {
      
            switch ( ind0(id) ) {
            case 1 : {	// Gamma^r_{l phi} V^l 
	                // or -Gamma^l_{r phi} V_l 
	        ind = ind0 ; 
	        ind.set(id) = 3 ;   // l = phi
	        tmp1 -= uu(ind) ; 
	        break ; 
            }
				
            case 2 : {	// Gamma^theta_{l phi} V^l 
	                // or -Gamma^l_{theta phi} V_l
	        ind = ind0 ; 
	        ind.set(id) = 3 ;   // l = phi
	        tmp2 -= uu(ind) ; 
	        break ; 
            }
				
            case 3 : {	// Gamma^phi_{l phi} V^l 
	                // or -Gamma^l_{phi phi} V_l
	        ind = ind0 ; 

	        ind.set(id) = 1 ;   // l = r
	        tmp1 += uu(ind) ; 

	        ind.set(id) = 2 ;   // l = theta
	        tmp2 += uu(ind) ; 
	        break ; 
            }
				
            default : {
	        cout << "Connection_fspher::p_divergence : index problem ! "
	            << endl ; 
	        abort() ;  
            }
            }
        }
        // There remains a division by tan(theta) and r:
        //----------------------------------------------
        tmp2.div_tant() ;
        tmp1 += tmp2 ;
        dz4 ? tmp1.div_r() : tmp1.div_r_inc2() ;

        cresu += tmp1 ; // the d/dr term...

    }

    // C'est fini !
    // -----------
    return resu ; 

}







