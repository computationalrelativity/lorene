/*
 *  Methods of template class Evolution_full
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

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/02/17 22:13:34  e_gourgoulhon
 * Suppressed declaration of global char[] evolution_C = ...
 *
 * Revision 1.2  2004/02/16 12:37:01  e_gourgoulhon
 * Method update: initialization of extended part of arrays val and the_time.
 *
 * Revision 1.1  2004/02/15 21:55:33  e_gourgoulhon
 * Introduced derived classes Evolution_full and Evolution_std.
 * Evolution is now an abstract base class.
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
#include <assert.h>



                    //-------------------------//
                    //      Constructors       //
                    //-------------------------//

                    
template<typename TyT> 
Evolution_full<TyT>::Evolution_full(const TyT& initial_value, 
                            double initial_time, int fact_resize_i)
      : Evolution<TyT>(initial_value, initial_time, 100), 
        fact_resize(fact_resize_i) 
{ }    
                                        

template<typename TyT> 
Evolution_full<TyT>::Evolution_full(const Evolution_full<TyT>& evo)
      : Evolution<TyT>(evo), 
        fact_resize(evo.fact_resize) 
{ }

                    
                    
                    //-----------------------//
                    //      Destructor       //
                    //-----------------------//

                    
template<typename TyT> 
Evolution_full<TyT>::~Evolution_full(){ }
                    
                    
                    //-----------------------//
                    //      Mutators         //
                    //-----------------------//

                    
template<typename TyT> 
void Evolution_full<TyT>::operator=(const Evolution_full<TyT>& ) {

    cerr << "void Evolution_full<TyT>::operator= : not implemented yet ! \n" ; 
    abort() ; 

}

template<typename TyT> 
void Evolution_full<TyT>::operator=(const Evolution<TyT>& ) {

    cerr << "void Evolution_full<TyT>::operator= : not implemented yet ! \n" ; 
    abort() ; 

}


                    
template<typename TyT> 
void Evolution_full<TyT>::update(const TyT& new_value, double new_time) {

    jlast++ ; 
         
    if (jlast == size) {  // re-organization of arrays val and the_time is necessary
    
        int size_new = fact_resize * size ; 

        TyT** val_new = new TyT*[size_new] ; 
        for (int j=0; j<size; j++) {
            val_new[j] = val[j] ; 
        }
        for (int j=size; j<size_new; j++) {
            val_new[j] = 0x0 ; 
        }
            
        double* the_time_new = new double[size_new] ;
        for (int j=0; j<size; j++) {
            the_time_new[j] = the_time[j] ; 
        }
        for (int j=size; j<size_new; j++) {
            the_time_new[j] = -1e20 ; 
        }
            
        size = size_new ;
        delete [] val ; 
        val = val_new ;  
        delete [] the_time ; 
        the_time = the_time_new ; 
            
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
const TyT& Evolution_full<TyT>::operator[](int j) const {

    assert(j >= 0) ;
    assert(j < size) ; 
    
    TyT* pval = val[j] ; 
    assert(pval != 0x0) ; 
    
    return *pval ; 

}                  
                    
template<typename TyT> 
double Evolution_full<TyT>::get_time(int j) const {

    assert(j >= 0) ;
    assert(j < size) ; 
        
    return the_time[j] ; 

}                  
                    
