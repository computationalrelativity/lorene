/*
 *  Methods of class Connection.
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

char connection_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.14  2004/01/22 16:15:53  e_gourgoulhon
 * First operational version of ricci().
 *
 * Revision 1.13  2004/01/19 16:57:44  e_gourgoulhon
 * First implementation of method ricci().
 * Not tested yet.
 *
 * Revision 1.12  2004/01/13 21:33:33  e_gourgoulhon
 * Corrected a bug in method p_derive_cov: inverted case CON and case COV.
 *
 * Revision 1.11  2004/01/04 20:57:51  e_gourgoulhon
 * -- Data member delta is now of type Tensor_sym (and no longer
 *    Tensor_delta).
 * -- Better handling of tensor symmetries in method p_derive_cov().
 *
 * Revision 1.10  2004/01/01 11:24:04  e_gourgoulhon
 * Full reorganization of method p_derive_cov: the main loop is now
 * on the indices of the *output* tensor (to take into account
 * symmetries in the input and output tensors).
 *
 * Revision 1.9  2003/12/30 22:58:27  e_gourgoulhon
 * -- Replaced member flat_conn (flat connection) by flat_met (flat metric)
 * -- Added argument flat_met to the constructors of Connection.
 * -- Suppressed method fait_ricci() (the computation of the Ricci is
 *    now devoted to the virtual method ricci()).
 * -- Implementation of methods fait_delta() and derive_cov().
 *
 * Revision 1.8  2003/12/27 14:59:05  e_gourgoulhon
 * Method derive_cov() suppressed.
 *
 * Revision 1.7  2003/10/16 14:21:36  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.6  2003/10/11 14:39:49  e_gourgoulhon
 * Suppressed declaration of unusued arguments in some methods.
 *
 * Revision 1.5  2003/10/06 13:58:46  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.4  2003/10/03 14:16:04  e_gourgoulhon
 * Added set_der_0x0 in some constructors.
 *
 * Revision 1.3  2003/10/02 21:32:06  e_gourgoulhon
 * Added constructor from Metric.
 * Added functions fait_delta and update.
 *
 * Revision 1.2  2003/10/01 15:42:49  e_gourgoulhon
 * still ongoing...
 *
 * Revision 1.1  2003/09/29 21:13:08  e_gourgoulhon
 * First version --- not ready yet.
 *
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
#include "metric.h"


                        //-----------------------//
                        //      Constructors     //
                        //-----------------------//


// Constructor ab initio

Connection::Connection(const Tensor_sym& delta_i, 
                       const Metric_flat& flat_met_i) 
                      : mp(&(delta_i.get_mp())),
		        triad(delta_i.get_triad()),
		        delta(delta_i), 
		        assoc_metric(false),
		        flat_met(&flat_met_i) {
                        
    assert( delta_i.get_valence() == 3 ) ; 
    assert( delta_i.sym_index1() == 1 ) ; 
    assert( delta_i.sym_index2() == 2 ) ;
    assert( delta_i.get_index_type(0) == CON ) ; 
    assert( delta_i.get_index_type(1) == COV ) ; 
    assert( delta_i.get_index_type(2) == COV ) ; 
		
    set_der_0x0() ; 
}		


// Standard constructor from a metric. 

Connection::Connection(const Metric& met,
                       const Metric_flat& flat_met_i) 
                      : mp(&(met.get_mp())),
		        triad(met.cov().get_triad()),
		        delta(*mp, CON, COV, COV, *triad, 1, 2),
		        assoc_metric(true),
		        flat_met(&flat_met_i) {
		
    fait_delta(met) ; 	// Computes delta

    set_der_0x0() ; 
}


// Copy constructor

Connection::Connection(const Connection& conn_i) : mp(conn_i.mp),
		triad(conn_i.triad),
		delta(conn_i.delta), 
		assoc_metric(conn_i.assoc_metric),
		flat_met(conn_i.flat_met) {
			
    set_der_0x0() ; 
	
}		


// Constructor for derived classes

Connection::Connection(const Map& mpi, const Base_vect& bi) : mp(&mpi),
		triad(&bi),
		delta(mpi, CON, COV, COV, bi, 1, 2),
		assoc_metric(false),
		flat_met(0x0){
		
    set_der_0x0() ; 
	
}		


	
                        //-----------------------//
                        //      Destructor       //
                        //-----------------------//

Connection::~Connection(){

    del_deriv() ; 
	
}

			//-----------------------------//
			//        Memory management    //
			//-----------------------------//

void Connection::del_deriv() const {

	if (p_ricci != 0x0) delete p_ricci ; 
	
	set_der_0x0() ; 
	
}

void Connection::set_der_0x0() const {

	p_ricci = 0x0 ; 
	
}


            //-----------------------------//
            //     Mutators / assignment   //
            //-----------------------------//


void Connection::operator=(const Connection& ci) {
	
	assert( triad == ci.triad ) ; 
	delta = ci.delta ; 
    flat_met = ci.flat_met ; 
	
	del_deriv() ; 

}	

void Connection::update(const Tensor_sym& delta_i) {

	assert(assoc_metric == false) ;
	
	assert(flat_met != 0x0) ; // to guarantee we are not in a derived class
	
        assert( delta_i.get_valence() == 3 ) ; 
        assert( delta_i.sym_index1() == 1 ) ; 
        assert( delta_i.sym_index2() == 2 ) ;
        assert( delta_i.get_index_type(0) == CON ) ; 
        assert( delta_i.get_index_type(1) == COV ) ; 
        assert( delta_i.get_index_type(2) == COV ) ; 
		
	delta = delta_i ; 
	
	del_deriv() ; 
	
}


void Connection::update(const Metric& met) {

	assert(assoc_metric == true) ;
	
	assert(flat_met != 0x0) ; // to guarantee we are not in a derived class
	
	fait_delta(met) ; 
	
	del_deriv() ; 
	
}



			//-----------------------------//
			//    Computational methods    //
			//-----------------------------//

// Covariant derivative
//---------------------

Tensor* Connection::p_derive_cov(const Tensor& uu) const {

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
    assert(flat_met != 0x0) ; 

    // Creation of the result (pointer)
    // --------------------------------
    Tensor* resu ;

    // If uu is a Scalar, the result is a Vector
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
	
    *resu = uu.derive_cov(*flat_met) ;   // Initialisation to the flat derivative 
	
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
 
        // Value of last index (derivation index)
        int k = ind1(valence1m1) ; 
        
        // Loop on the number of indices of uu 
        for (int id=0; id<valence0; id++) {
            
            ind = ind0 ;
                
            switch( uu.get_index_type(id) ) {
                
                case CON : {
                    for (int l=1; l<=3; l++) {
                        ind.set(id) = l ; 
                        cresu += delta(ind0(id), k, l) * uu(ind) ;
                    }
                    break ; 
                }
                
                case COV : {
                    for (int l=1; l<=3; l++) {
                        ind.set(id) = l ; 
                        cresu -= delta(l, k, ind0(id)) * uu(ind) ;
                    }
                    break ; 
                }
                
                default : {
                    cerr << 
                    "Connection::p_derive_cov : unexpected type of index !\n" ;
                    abort() ; 
                    break ; 
                }
                
            }   // end of switch on index type 
                
        }   // end of loop on the number of indices of uu               
        
    }   // end of loop on all the components of the output tensor

    // C'est fini !
    // ------------
  
    return resu ; 
  
} 


Tensor* Connection::p_divergence(const Tensor& ) const {

	cout << "Connection::p_divergence : not implemented yet !" << endl ; 
	abort() ; 
	return 0x0 ;
} 


const Tensor& Connection::ricci() const {

    if (p_ricci == 0x0) {  // a new computation is necessary
    
        if (assoc_metric) {     // The Ricci tensor is symmetric if the
                                // connection is associated with some metric
            p_ricci = new Sym_tensor(*mp, COV, *triad) ; 
        }
        else {
            p_ricci = new Tensor(*mp, 2, COV, *triad) ; 
        }
        
        const Tensor& d_delta = delta.derive_cov(*flat_met) ; 
                
        for (int i=1; i<=3; i++) {
        
            int jmax = assoc_metric ? i : 3 ; 
            
            for (int j=1; j<=jmax; j++) {

                Scalar tmp1(*mp) ;
                tmp1.set_etat_zero() ; 
                for (int k=1; k<=3; k++) {
                    tmp1 += d_delta(k,i,j,k) ; 
                } 
                
                Scalar tmp2(*mp) ;
                tmp2.set_etat_zero() ; 
                for (int k=1; k<=3; k++) {
                    tmp2 += d_delta(k,i,k,j) ; 
                } 
                
                Scalar tmp3(*mp) ;
                tmp3.set_etat_zero() ; 
                for (int k=1; k<=3; k++) {
                    for (int m=1; m<=3; m++) {
                        tmp3 += delta(k,k,m) * delta(m,i,j) ; 
                    }
                } 
                tmp3.dec_dzpuis() ;  // dzpuis 4 -> 3
                
                Scalar tmp4(*mp) ;
                tmp4.set_etat_zero() ; 
                for (int k=1; k<=3; k++) {
                    for (int m=1; m<=3; m++) {
                        tmp4 += delta(k,j,m) * delta(m,i,k) ; 
                    }
                } 
                tmp4.dec_dzpuis() ;  // dzpuis 4 -> 3
                
                p_ricci->set(i,j) = tmp1 - tmp2 + tmp3 - tmp4 ; 
                
            }
        }

    }
	
    return *p_ricci ; 
	
}



void Connection::fait_delta(const Metric& gam) {

    assert(flat_met != 0x0) ; 
        
    const Tensor& dgam = gam.cov().derive_cov(*flat_met) ; 
    
    for (int k=1; k<=3; k++) {
        for (int i=1; i<=3; i++) {
            for (int j=1; j<=i; j++) {
                Scalar& cc = delta.set(k,i,j) ; 
                cc = 0 ; 
                for (int l=1; l<=3; l++) {
                    cc += gam.con()(k,l) * ( 
                        dgam(l,j,i) + dgam(i,l,j) - dgam(i,j,l) ) ; 
                        
                }
                cc = 0.5 * cc ; 
            }
        }
    }


}  












