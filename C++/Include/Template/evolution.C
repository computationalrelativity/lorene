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
 * Revision 1.2  2004/02/15 21:55:33  e_gourgoulhon
 * Introduced derived classes Evolution_full and Evolution_std.
 * Evolution is now an abstract base class.
 *
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
Evolution<TyT>::Evolution(const TyT& initial_value, double initial_time,
    int size_i)
      : size(size_i),
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
Evolution<TyT>::Evolution(const Evolution<TyT>& evo)
      : size(evo.size),
        jlast(evo.jlast) {

    val = new TyT*[size] ; 
    
    for (int j=0; j<size; j++) {
        if (evo.val[j] != 0x0) {
            val[j] = new TyT( *(evo.val[j]) ) ; 
        }
        else {
            val[j] = 0x0 ; 
        }
    }
    
    the_time = new double[size] ; 
    
    for (int j=0; j<size; j++) {
        the_time[j] = evo.the_time[j] ; 
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
void Evolution<TyT>::operator=(const Evolution<TyT>& ) {

    cerr << "void Evolution<TyT>::operator= : not implemented yet ! \n" ; 
    abort() ; 

}


                    
                    
template<typename TyT> 
TyT Evolution<TyT>::operator()(double ) const {

    cerr << 
    "Evolution<TyT>::operator()(double ) const : not implemented yet !"
        << endl ; 
    abort() ; 

}                  
                    
