/*
 *  Methods of template class Evolution
 *
 *    (see file evolution.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon & Jerome Novak
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

char evolution_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/02/13 15:53:20  e_gourgoulhon
 * New (template) class for time evolution.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h" 

// C headers
#include <stdlib.h>
#include <assert.h>



                    //-------------------------//
                    //      Constructors       //
                    //-------------------------//

                    
template<typename TyT> 
Evolution<TyT>::Evolution(const TyT& initial_value, double initial_time)
      : size(100),
        fixed_size(false),
        jlast(0) {

    val = new TyT*[size] ; 
    
    val[0] = new TyT(initial_value) ; 

    for (int j=1; j<size; j++) {
        val[j] = 0x0 ; 
    }
    
    the_time = new double[size] ; 
    
    the_time[0] = initial_time ; 

    for (int j=1; j<size; j++) {
        the_time[j] = -1e20 ; 
    }
    
    
}                    
                    

template<typename TyT> 
Evolution<TyT>::Evolution(const TyT& initial_value, double initial_time, 
    int nstored)
      : size(nstored),
        fixed_size(true),
        jlast(0) {

    assert(nstored > 0) ; 
    
    val = new TyT*[size] ; 
    
    val[0] = new TyT(initial_value) ; 

    for (int j=1; j<size; j++) {
        val[j] = 0x0 ; 
    }
    
    the_time = new double[size] ; 
    
    the_time[0] = initial_time ; 

    for (int j=1; j<size; j++) {
        the_time[j] = -1e20 ; 
    }
    
    
}                    
                    

template<typename TyT> 
Evolution<TyT>::Evolution(const Evolution<TyT>& a_in)
      : size(a_in.size),
        fixed_size(a_in.fixed_size),
        jlast(a_in.jlast) {

    val = new TyT*[size] ; 
    
    for (int j=0; j<size; j++) {
        if (a_in.val[j] != 0x0) {
            val[j] = new TyT( *(a_in.val[j]) ) ; 
        }
        else {
            val[j] = 0x0 ; 
        }
    }
    
    the_time = new double[size] ; 
    
    for (int j=0; j<size; j++) {
        the_time[j] = a_in.the_time[j] ; 
    }
    
    
}                    
                    

                    
                    
                    //-----------------------//
                    //      Destructor       //
                    //-----------------------//

                    
template<typename TyT> 
Evolution<TyT>::~Evolution(){

    for (int j=0; j<size; j++) {
        if (val[j] != 0x0) delete val[j] ; 
    }
    
    delete [] val ;
    delete [] the_time ; 
    
}
                    
                    
                    //-----------------------//
                    //      Mutators         //
                    //-----------------------//

                    
template<typename TyT> 
void Evolution<TyT>::operator=(const Evolution<TyT>& a_in) {

    assert(a_in.fixed_size == fixed_size) ; 
    
    cerr << "void Evolution<TyT>::operator= : not implemented yet ! \n" ; 
    abort() ; 

}


                    
template<typename TyT> 
void Evolution<TyT>::update(const TyT& new_value, double new_time) {

    jlast++ ; 
         
    if (jlast == size) {  // re-organization of arrays val and the_time is necessary
    
        if (fixed_size) {   // Permutation of values inside array val 
        
            delete val[0] ; 
            
            for (int j=0; j<size-1; j++) {
                val[j] = val[j+1] ; 
            }
        
            for (int j=0; j<size-1; j++) {
                the_time[j] = the_time[j+1] ; 
            }
            
            jlast-- ;
        
        }
        else {  // Enlargement of array val : doubling the size

            int size_new = 2 * size ; 
            TyT** val_new = new TyT*[size_new] ; 
            for (int j=0; j<size; j++) {
                val_new[j] = val[j] ; 
            }
            
            double* the_time_new = new double[size_new] ;
            for (int j=0; j<size; j++) {
                the_time_new[j] = the_time[j] ; 
            }
            
            size = size_new ;
            delete [] val ; 
            val = val_new ;  
            delete [] the_time ; 
            the_time = the_time_new ; 
            
        }
    }
    else {
        assert( jlast < size ) ; 
    }
    
    val[jlast] = new TyT( new_value ) ; 
    the_time[jlast] = new_time ; 

}                   
                    
                    
                    //-----------------------//
                    //      Accessors        //
                    //-----------------------//

                 
template<typename TyT> 
const TyT& Evolution<TyT>::operator[](int pos) const {

    assert(pos >= 0) ;
    assert(pos < size) ; 
    
    TyT* pval = val[pos] ; 
    assert(pval != 0x0) ; 
    
    return *pval ; 

}                  
                    
template<typename TyT> 
TyT Evolution<TyT>::operator()(double ) const {

    cerr << 
    "Evolution<TyT>::operator()(double ) const : not implemented yet !"
        << endl ; 
    abort() ; 

}                  
                    
